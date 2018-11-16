#######
# Clare Bycroft, May 2016
#
# This script contains functions auxiliary to the main R scripts, and which are useful for multiple other R scripts involved in sample QC. It also encodes some directory names that may need to be changed in future versions. It does call functions involved in other files, e.g batch2sample.R, but those might change, so we don't hard code the links to those files in here.
#
######

library(stringr)
library(data.table)
library(RColorBrewer)
library(grid)
library(gridBase)  # install into own 
library(sp)
library(png)
library(dplyr)
library(rhdf5)
library(abind)
library(qqman)
library(quadprog)
library(hexbin)
library(latticeExtra)
library(viridis)
library(ggplot2)

plink = '/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink'
qctool = '/apps/well/qctool/1.5/qctool'
bgenix = '/apps/well/bgenix/20160708/bgenix'
ldstore = '/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/ldstore_v1.1_x86_64/ldstore'
#opar=par()
# set some global options
options(stringsAsFactors = FALSE,check.names=FALSE) # these apply to reading and writing data


# set some global variables
baseDataDir="/well/ukbiobank/expt/V2_QCed.SNP-QC/data"
#baseSampleQCDir="/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing"
baseSampleQCDir="/well/ukbiobank/qcoutput.V2_QCed.sample-QC"
GenotypesForRelease="/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr"

# downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes on 10/8/2018
chromLengthsFile = paste0(baseSampleQCDir,"/data/hg19.chrom.sizes")

sexChroms = c(23,24,25,26,25)
names(sexChroms) = c("X","Y","XY","MT","PAR")
sexChroms2 = c(23,24,25,25,26)
names(sexChroms2) = c("X","Y","XPAR1","XPAR2","MT")


# set par boundaries (from ucsc browser, hg19)
par1 = 2699520
par2 = 154931044



get.sample.table.name <- function(Batch){    
    paste0(Batch,"/V2_QCed_",Batch,"_Sample_Table_Pheno.csv")    
}

list.references <- function(Batch="Batch_b001"){    
    print(paste0("Listing reference sample names based on ",Batch))
    csvFile = get.sample.table.name(Batch)
    info = get.batch.info(Batch,paste0(baseDataDir,"/",csvFile))    
    references = unique(info$PIID[get.batch.references(info)])
    return(references)
}

list.all.references <- function(){
    batches = all.batches()
    refList = sapply( batches,function(b) list.references(b) )
}

find.my.batch <- function(PIID,batchFile){
    #rand=round(runif(1,100,199),0)
    #write.table(PIID,file=paste0("tmp.",rand),quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    #out = system( paste0("grep -wf ",paste0("tmp.",rand)," ",batchFile),intern=TRUE)
    #bs = t(sapply(out,function(x) str_split(x," ")[[1]],USE.NAMES=F))
    bs = read.table(batchFile,stringsAsFactors=FALSE,header=FALSE)
    ids = bs[,1]    
    index = match(PIID,ids)
    b = bs[index,2]
    
    if(sum(bs[index,1]!=PIID) > 0 )  print("WARNING: problem with function find.my.batch in auxFunctions.R")
#    system(paste0("rm tmp.",rand))
    
    return(b)
}

batchInfoFields = colnames(get.batch.info("Batch_b001",paste0(baseDataDir,"/",get.sample.table.name("Batch_b001"))))
quantInfoFields = batchInfoFields[c(7:8,11,17,19,20,28,29,32,34:42,44,46:52,54)]

read.multiple.batch.info <- function(fields,theseBatches=NULL,exclRef=TRUE,printBatches=TRUE){
    if(is.null(theseBatches)) theseBatches=all.batches()
    if(!"PIID"%in%fields) fields = c("PIID",fields)

    if(sum(!fields%in%batchInfoFields) > 0 ){
        print( paste0("WARNING: ",fields[!fields%in%batchInfoFields]," not an available field!" ))
        print("Choose from: ")
        print(batchInfoFields)
        
        fields = fields[fields%in%batchInfoFields]
    }
    
    info = sapply(theseBatches,function(b){
        if(printBatches) print(b)
        csvFile = get.sample.table.name(b)        
        out = get.batch.info(b,paste0(baseDataDir,"/",csvFile))
        
        if(length(intersect(colnames(out),c("Well","Plate.Name"))) ==2 ) {
            colnames(out)[colnames(out)=="Well"] = "Submitted.Well"
            colnames(out)[colnames(out)=="Plate.Name"] = "Submitted.Plate.Name"            
        }

        if(exclRef) out = out[out$Sample.Source=="Customer",]
        
       # print(fields[!fields%in%colnames(out)])
        out = out[,fields]
    },simplify=F)

    info = rbindlist(info)
    return(info)
}


## functions to help with plotting

# get log scale tick marks
getLogTicks <- function(powers=-c(10:0),nticks=5,base=10,lims=NULL){
    seq1 = base^powers
    smallTicks = unlist(sapply(1:(length(seq1)-1),function(i) seq(seq1[i],seq1[i+1],by=diff(seq1[i:(i+1)])/nticks)[-c(1,nticks+1)],simplify=FALSE))
    bigTicks = seq1
    allTicks = c(bigTicks,smallTicks)
    if(!is.null(lims)){
        lineRange = c(min(allTicks[allTicks>=lims[1]]),max(allTicks[allTicks<=lims[2]]))
        return(list(bigTicks,smallTicks,lineRange))
    } else {
        return(list(bigTicks,smallTicks))                            
    }
}

# get logit scale tick marks
getLogitTicks <- function(range=c(0.0001,0.9999),nticks=5,base=10,lims=NULL,stepSize=0.3){
    edg = range
    seq1 = inv.logit(seq(logit(edg[1]),logit(edg[2]),by=-logit(stepSize)))
    smallTicks = unlist(sapply(1:(length(seq1)-1),function(i) seq(seq1[i],seq1[i+1],by=diff(seq1[i:(i+1)])/nticks)[-c(1,nticks+1)],simplify=FALSE))
    bigTicks = seq1
    allTicks = c(bigTicks,smallTicks)
    if(!is.null(lims)){
        lineRange = c(min(allTicks[allTicks>=lims[1]]),max(allTicks[allTicks<=lims[2]]))
        return(list(bigTicks,smallTicks,lineRange))
    } else {
        return(list(bigTicks,smallTicks))                            
    }
}

                                        # List of colourblind-friendly colours from http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
getColBlindCol <- function(want=NULL){
    o=rgb2hsv(230,159,0) # orange
    s=rgb2hsv(86,180,233) # sky blue
    g=rgb2hsv(0,158,115) # bluegreen
    y=rgb2hsv(240,228,66) # yellow
    b=rgb2hsv(0,114,178) # blue
    v=rgb2hsv(213,94,0) # vermillion
    r=rgb2hsv(204,121,167) # red/pink
    if(is.null(want)) want = c("y","o","v","r","b","s","g")
    cols = sapply(want,function(x) {
        m=get(x)
        return(hsv(m[1],m[2],m[3]))
        })
    return(cols)
}

giantColour = "#0072B2" # From colourblind palette!
genoColour = "#009E73" # From colourblind palette!
impColour = getColBlindCol("r") # From colourblind palette!

sexCols=getColBlindCol(c("r","g","v","b"))
#sexCols = c("purple","darkgreen","blue","red")
names(sexCols) = c("M","F","MtoF","FtoM")


# get hsv properties of standard colours
colorsHsv=rgb2hsv(col2rgb(colors()))
primaryColors = c("red","blue","green","darkgreen","yellow","cyan","orange","purple","maroon4","darkgreen","darkblue","skyblue","orange3")
colorsSat = colors()[colorsHsv[2,]>0.8]

# get set of distinct colours
getColoursDistant3 <- function(n,huFrac=0.6,startHue=0,minLum=20,maxLum=80,chroma=10){
  if(huFrac<1/n) huFrac=1/n
  h=seq(360/n+startHue,360+startHue,by=360/n)
  h=h%%360
  l=seq((maxLum-minLum)/(n*huFrac),maxLum,by=(maxLum-minLum)/(n*huFrac))
  print(l)
  myCols = sapply(1:n,function(x) {
        hcl(h=h[x],l=l[x%%length(l)+1],c=chroma,fixup=TRUE)
      })
  return(myCols)
  #pizza(myCols)
}



# get colour scales from numeric vector
colour.scale <- function(y,colourSet= c("green","yellow","red"),nBreaks=200,fixedLims=NULL,colSpace="Lab",...){
    x = y[!is.na(y)] # remove NAs
    Colors <- colorRampPalette(colourSet,space=colSpace,...)(nBreaks)
    if (!is.null(fixedLims)) scale <- seq(fixedLims[1],fixedLims[2]+0.0001,length.out=nBreaks) else scale <- seq(min(x),max(x)+0.0001,length.out=nBreaks)
    cols = rep(NA,length(y))
    cols[!is.na(y)] = sapply(x,function(g) Colors[which(scale>g)[1]])
    
return(cols)
}

## function to order samples by number of occurances
order.by.number.occurrences <- function(group) {
    x = as.numeric(factor(group))
    x[is.na(x)] = max(x,na.rm=TRUE)+1
    order = sort(table(x)[x],decreasing=TRUE,index.return=TRUE)
    order = order$ix
    return(order)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

My.add.alpha2 <- function(col, alpha=1, value=1,saturation=1){
  if(missing(col))
    # saturation and value are absolute values for the output (rather than fold changes)
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb), 2, 
                     function(x) {
                       Hsv=rgb2hsv(x)
                       hsv(Hsv[1],saturation,value,alpha=alpha)
                     })
}


get.place.of.birth <- function(Info){
    c = colnames(Info)
    if((!"Ethnic.background"%in%c) | (!"Country.of.birth..UK.elsewhere."%in%c) | (!"Country.of.Birth..non.UK.origin."%in%c) ) stop("WARNING: not all necessray columns exist in input data. ")
        
    Info$Eth2 = transform.ethnicity(Info$Ethnic.background)
    Info$Country.of.birth = Info$Country.of.birth..UK.elsewhere.
    Info$Country.of.birth[Info$Country.of.birth..UK.elsewhere.%in%c("Elsewhere","",NA)] = Info$Country.of.Birth..non.UK.origin.[Info$Country.of.birth..UK.elsewhere.%in%c("Elsewhere","",NA)]
    Info$Country.of.birth[is.na(Info$Country.of.birth)] = "Other/Unknown"
    
    return(Info$Country.of.birth)
}


generateColours <- function(x){
    # input a list of unique countries (or any other string vector)
    countryColours = sample(unique(ethnicity2col),length(x)+1,replace=T)
    names(countryColours) = c(x,"Other/Unknown")
    symb = sample(0:6,length(countryColours),replace=T)
    dupes = table(countryColours)
    for(d in names(dupes) ){
        these = countryColours==d
        if(dupes[d]==1) next
        if((dupes[d]>1)&(length(unique(symb[these]))==sum(these))) next
        
        while((dupes[d]>1)&(length(unique(symb[these])) < sum(these))){
            #print(d)
            symb[these] = sample(c(0:6),sum(these),replace=T)
        }
    }
    names(symb) = names(countryColours)
    cols = countryColours
    cols["Other/Unknown"]="gray"
    chars = symb
    return(list(cols,chars))
}


ethnicityOrder = c(10,15,20,6,14,18,8,12,7,3,2,9,11,4,21,22,23,16,5,1,17,13,19) # order based on "ethnicities" object in ukbcolors.R


makeScaleWithHist <- function(values,colourSet,scalenum=100,nbreaks=200,fixedLims=NULL,binWidth=NULL,xlims=NULL,bgColour="transparent",xlab=NULL){
    opar=par()
    par(mar=c(2,0,0,0),mgp=c(1,1,0))
    if(is.null(fixedLims)) colindex<-matrix(seq(min(values),max(values),length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
    if(!is.null(fixedLims)) colindex<-matrix(seq(fixedLims[1],fixedLims[2],length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
   
    if(!is.null(binWidth)) getBreaks <- seq(min(values),max(values)+binWidth,by=binWidth) else getBreaks <- seq(min(values),max(values),length.out=nbreaks)
   
    h <- hist(values,breaks=getBreaks,plot=F)
    histMax <- max(h$density) + max(h$density)/20
    if(is.null(xlims)) xlims=c(min(values),max(values))
    
    hist(values,breaks=getBreaks,xlim=xlims,axes=F,xlab=xlab,ylab=NULL,main=NULL,probability=T,border="gray")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bgColour,border=NA)
    #axis(1,cex.axis=0.5,tck=-0.1)
    
    colours <- colour.scale(colindex[,1],colourSet)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=length(colindex))
    #scalephysicalpos <- seq(1,nbreaks,length.out=length(colindex))
    scalephysicalpos <- seq(1,nbreaks,length.out=length(colindex))
   # image(1:nbreaks,1,colindex,xaxt="n",xlab="",yaxt="n",ylab="",col=colours,zlim=range(colindex))
     image(scalelocs,c(-max(h$density)/20,histMax),colindex,xaxt="n",xlab="",yaxt="n",ylab="",col=colours,zlim=range(colindex),axes=F,add=T)     
    # axis(1,las=1,cex.axis=0.6,tcl=-0.2,padj=-3)
    # axis(1,las=1,cex.axis=0.6)
    
    
    par(mgp=c(0,0,0))
    if(is.null(fixedLims)) hist(values,breaks=getBreaks,col=add.alpha("white",0.7),axes=F,xlab=NULL,ylab=NULL,main=NULL,border="transparent",probability=T,add=T,xlim=xlims)
    if(!is.null(fixedLims)) hist(values,breaks=getBreaks,col=add.alpha("white",0.8),
                                 axes=F,xlab=NULL,ylab=NULL,main=NULL,border="transparent",xlim=c(fixedLims[1],fixedLims[2]),probability=T,add=T)
    axis(1,cex.axis=0.8,tck=-0.03,at=signif(seq(xlims[1],xlims[2],length.out=7),2))
     #box(lwd=0.5)    
    par(mar=opar$mar,mgp=opar$mgp)
}



## function for plotting PCA in grid

plot.grid.PCs <- function(PCs,n=NULL,maxPC=NULL,marg=0.05,gridMar=0.08,marLab=1.2,shapes=1,colours="black",...){

  plot.new()
  grid.newpage()
  pushViewport(viewport(marg,marg,width=1-2*marg,height=1-2*marg,just=c("left","bottom")))
  #grid.rect()
  #grid.show.layout(gl) 
  if(is.null(n)) n = ceiling(ncol(PCs)/2)
  if(is.null(maxPC)) maxPC = n*2
  print(n)
  
  gl <- grid.layout(n, n)
  pushViewport(viewport(layout=gl))     
  #
  for(i in 1:n) { 
      for(j in 1:n) {
          
          if(i==j) next
              
          print(c(i,j))

          if(i < j){
              pc1 = i+n
              pc2 = j+n
          } else {
              pc1 = i
              pc2 = j
          }
          print(c(pc1,pc2))

          if((pc1>maxPC)|(pc2>maxPC)) next

          pushViewport(viewport(layout.pos.row=i, layout.pos.col=j, name="foo")) 
                                        #  grid.rect() 
          par(fig=gridFIG(),new=TRUE,mgp=c(2,0.3,0),mar=gridMar*c(1,1,1,1),...)
          
          ylabel=NA
          xlabel=NA
          xaxty="n"
          yaxty="n"        
          
             # if pc1 or pc2 is bigger than the number of columns, move on
          if((pc2 > ncol(PCs))|(pc1 > ncol(PCs))) { upViewport(); next }

          plot(range(PCs[,pc2]),range(PCs[,pc1]),col=NA,xaxt=xaxty,yaxt=yaxty,xlab=xlabel,ylab=ylabel)    
          points(PCs[,pc2],PCs[,pc1],cex=2/n,xpd=NA,col=colours,pch=shapes,xpd=NA,lwd=0.5)    

          if(pc2==1) { axis(2,tcl=-0.2); mtext(colnames(PCs)[pc1],2,xpd=NA,line=marLab) } 
          if(pc1==n) { axis(1,tcl=-0.2,padj=-1); mtext(colnames(PCs)[pc2],1,xpd=NA,line=marLab) }

          if(pc1==(n+1)) { axis(3,tcl=-0.2,padj=-1); mtext(colnames(PCs)[pc2],3,xpd=NA,line=marLab) } 
          if(pc2==2*n) { axis(4,tcl=-0.2); mtext(colnames(PCs)[pc1],4,xpd=NA,line=marLab) }

          upViewport()
      }
  }
    print(colnames(PCs))

} 

#png(paste0(OutDir,'/plots/',OutputFile,'-grid-%02d.png'),width=1500,height=1500,res=150)

#for(l in seq(1,nPCs,10)){
#    print(l)
#    plot.grid.PCs(as.data.frame(PCs[Order,paste0("PC",c(l:(l+9)))]),shapes=Chars[Order],colours=Colors[Order])
#}

#dev.off()


# Plot a heat map!
plotMixtureHeat <- function(toPlot,filename,title="Linear coefficients",cap=T,capValue=NULL,colourSet=NULL,Cex=1,ylab="PC",xlab="Ethnic background",negativeScale=T,nCols=100,fixedRange=NULL,plotPoints=F,rowLabels=NULL,colLabels=NULL,Legend=NULL,scaleTitle=NA,x.cex=1,y.cex=1,labYoff=1,ylabMar=17,scaleWidth=0.08,...){

    toPlot2 <- toPlot

    if((negativeScale==F)&(is.null(colourSet))) colourSet=c("white","yellow","red")
    if(is.null(colourSet)) colourSet = c("red","white","blue")

    colourScale=colorRampPalette(colourSet,interpolate="linear")(nCols+(nCols+1)%%2)# make sure there's always an odd number of colours so 0 can be in the centre.
    
    
    if(cap){
        if(is.null(capValue)){
            cap=mean(abs(toPlot)) + 2*sd(abs(toPlot))
        } else {
            cap = capValue
        }
        print('capping matrix')
        toPlot2[toPlot2>cap] <- cap
        toPlot2[toPlot2<(-cap)] <- -cap
        
        colourScale=c(My.add.alpha2(colourScale[1],value=0.2),colourScale,My.add.alpha2(colourScale[length(colourScale)],value=0.5))
        print( head(colourScale) )
        print(range(toPlot2))
    }
                               #cap=5
    

    Range <- max(abs(toPlot2),na.rm=T)
    if(negativeScale==F) Range=c(min(toPlot2,na.rm=T),max(toPlot2,na.rm=T))
    if(negativeScale==T) Range=c(-Range,Range)
    if(!is.null(fixedRange)) Range=fixedRange

    if(is.null(rowLabels)) rowLabels <- rownames(toPlot2)
    if(is.null(colLabels)) colLabels <- colnames(toPlot2)


    if(plotPoints==T){      
        rowShapes <- as.numeric(Legend[rownames(toPlot2),"shape"])
        rowColours <- Legend[rownames(toPlot2),"colour"]
        
        colShapes <- as.numeric(Legend[colnames(toPlot2),"shape"])
        colColours <- Legend[colnames(toPlot2),"colour"]
    } 


        
    if(grepl("png",filename)) png(filename,...)
    if(grepl("pdf",filename)) pdf(filename,...)
    
    #png(paste0(OutDir,"/plots/test.png"),height=1200,width=2000,res=150)
    layout(mat=matrix(c(1,2),1,2),widths=c(1-scaleWidth,scaleWidth))
                                        #par(mar=c(5,11,2,0.5),mgp=c(10,1,0),cex=Cex)
    par(mar=c(5,ylabMar,2,0.5),cex=Cex)

    image(as.matrix(toPlot2),col=colourScale,zlim=Range,axes=F,main=title,
          xlab=xlab,ylab=ylab)
    axis(1,at=seq(0,1,length.out=nrow(toPlot2)),labels=rowLabels,las=1,tick=T,cex.axis=x.cex)
    axis(2,at=seq(0,1,length.out=ncol(toPlot2)),labels=colLabels,las=2,tick=F,cex.axis=y.cex,line=labYoff)
    
    if(plotPoints==T){
        points(x=seq(0,1,length.out=nrow(toPlot2)),y=rep(line2user(1,1)-(line2user(2,1)-line2user(1,1))/2,nrow(toPlot2)),col=rowColours,pch=rowShapes,xpd=NA,bg=rowColours,lwd=1.5)
        points(y=seq(0,1,length.out=ncol(toPlot2)),x=rep(line2user(1,2)-(line2user(labYoff,2)-line2user(1,2))/2,nrow(toPlot2),ncol(toPlot2)),col=colColours,pch=colShapes,xpd=NA,bg=colColours,lwd=1.5)
    }
                                        #1300,1933
    # The scale!
    par(mar=c(5,0,2,4),c(3, -1, 0))
    colindex<-t(matrix(seq(Range[1],Range[2],length.out=nCols),ncol=1,nrow=nCols))
    image(1,1:nCols,colindex,xaxt="n",yaxt="n",xlab="",ylab="",col=colourScale,zlim=range(colindex),main=scaleTitle,cex.main=0.8,xpd=NA)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=10)
    axis(2,at=seq(1,nCols,length.out=10),labels=signif(scalelocs,2),
         las=2,cex.axis=0.8,side=4,tck=-0.2,hadj=0.5)
#    axis(2,at=(nCols+1)/2,labels=0,
#         las=2,cex.axis=0.8,side=4,tck=-0.2,hadj=0.5)

    dev.off()
}



# get x and y coordinates of current plot line-heights (used in the above function) 
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}





## Infer kinship classes from kinship estimates ===> From KING.
degrees = c(1/(2^(3/2)),1/(2^(5/2)),1/(2^(7/2)),1/(2^(9/2)))
                                        #colDef = c("sibs"="red","dupe/twins"="pink","2nd degree"="green","3rd degree"="blue","parent/child"="purple")
colDef = c("sibs"="blue","dupe/twins"="red","2nd degree"="green","3rd degree"="yellow","parent/child"="purple")
shapeDef = c("sibs"=1,"dupe/twins"=1,"2nd degree"=1,"3rd degree"=1,"parent/child"=2)
kinNames = c("sibs"="Full siblings","dupe/twins"="Monozygotic twins","2nd degree"="2nd degree","3rd degree"="3rd degree","parent/child"="Parent-offspring")
    
kinWeights = c(6,10,4,2,8); names(kinWeights)=names(colDef) # weights for graphing 

kinshipLegendOrder = rev(names(kinWeights)[order(kinWeights)])

    
get.kin.classes <- function(kin,ibs0Threshold=0.002 ){
    
    ## get indexes of kinship classes
    classes = rep(NA,dim(kin)[1])
    p = kin$Kinship
    q = kin$IBS0
    classes[ p > degrees[1]  ] = "dupe/twins"
    classes[(p <= degrees[1] ) & ( q >= ibs0Threshold ) & (p > degrees[2] )] = "sibs"    
    classes[( p <= degrees[2] ) & ( p > degrees[3] )] = "2nd degree"
    classes[( p  <= degrees[3] ) & (p  > degrees[4] ) ] = "3rd degree"
    classes[(p <= degrees[1] ) & ( q < ibs0Threshold )] = "parent/child"

    if( sum(is.na(classes)) > 0 ) print("Some pairs unable to be classified...")
    return(classes)

}


# Get labels for family ids in network graphs
getFamLabels <- function(inds,myFams,myGraph){
    famIDs = myFams$membership[inds]
    labelInds = sapply(unique(famIDs[!is.na(famIDs)]),function(fid){
        is = inds[famIDs==fid]; is=is[!is.na(is)]
        degs = degree(myGraph)[is]
        keep = is[degs==min(degs)][1]
        print(length(is))        
        return(keep)
    })
    famIDs[!inds%in%labelInds]=NA
    names(famIDs) = inds
    famIDs = famIDs[V(myGraph)$name] # order by the vertices in the graph
    return(famIDs)
}


# function to parse command input strings with numbers separated by commas and colons, such as "1:2,5:10"

parse.range.string <- function(s){
    s = str_split(s,",")[[1]]
    colon = grepl(":|-",s)
    if(sum(colon)>0) {
        o = unlist(sapply(grep(":|-",s,value=TRUE),function(x) {
            t = as.numeric(str_split(x,":|-")[[1]])
            out = t[1]:t[2]
        },USE.NAMES=FALSE))
        s = unique(c(as.numeric(s[!colon]),o))
    }
    s=as.numeric(s)
    return(s)
}

read.recomb <- function(sex=FALSE,chrom=NULL){
    print('Reading recombination rate files from /well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr*')

    recombrates = list()
    
    if(is.null(chrom) ) recombChroms = c(1:22,"X_nonPAR","X_PAR1","X_PAR2") else recombChroms=as.character(parse.range.string(chrom))
    if(length(intersect(recombChroms,sexChroms))>0) sex = TRUE
        
    if(sex) recombChroms = c("X_nonPAR","X_PAR1","X_PAR2")
    
    for( i in recombChroms ){    
        print(i)
        r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
        recombrates[[i]] = r
    }
    recombrates[["23"]] = recombrates[["X_nonPAR"]]
    recombrates[["25"]] = Reduce(rbind,list(recombrates[["X_PAR1"]],recombrates[["X_PAR2"]]))
    
    return(recombrates)
}


read.SNPQC.files <- function(otherFiles=NULL,hweThreshold=10^-100,type="autosome",justSNPs=TRUE){

    arrayFile = paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".array_test.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    batchHweFile = paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".hwe.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    mfFile =  paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".malefemale.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    plateFile =  paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".plate_effect.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    batchFile =  paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".batch_effect.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    concordanceFile = paste0(baseDataDir,"/SNP-QC_summary/V2_QCed.",type,".controls_discordance.0.05.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt")
    imageArtefactFile = paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/SelectSNPs/probesets_in_subimage_37555.txt")
    hweFile = paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-",type,"-oxfordqc.hwe")
        
    arraySNPs = read.table(arrayFile,header=FALSE,stringsAsFactors=FALSE)
    batchHweSNPs = read.table(batchHweFile,header=FALSE,stringsAsFactors=FALSE)
    mfSNPs = read.table(mfFile,header=FALSE,stringsAsFactors=FALSE)
    plateSNPs = read.delim(plateFile,header=FALSE,stringsAsFactors=FALSE,sep=" ",fill = TRUE,col.names=c("Batch","SNP",paste0("Plate",1:100)))
    batchSNPs = read.table(batchFile,header=FALSE,stringsAsFactors=FALSE)

    concordanceSNPs = read.table(concordanceFile,header=FALSE,stringsAsFactors=FALSE)

    imageSNPs = read.table(imageArtefactFile,header=TRUE,stringsAsFactors=FALSE)[,"affy_snp_id"]
    
    hwe = read.table(hweFile,stringsAsFactors=FALSE,header=TRUE)
    hwe = tbl_df(hwe)
    hwe = filter(hwe,TEST=="ALL")
    hweSNPs = hwe[(hwe$P < hweThreshold) & (!is.na(hwe$P)),]

    
    otherSNPs=c()
    for(otherFile in otherFiles){
        other = read.table(otherFile,stringsAsFactors=FALSE,header=FALSE)[,1]
        otherSNPs = unique(c(other,otherSNPs))
    }

    if(justSNPs) snpList = list("arraySNPs"=arraySNPs[,1],"batchHweSNPs"=batchHweSNPs[,2],"mfSNPs"=mfSNPs[,2],"concordanceSNPs"=concordanceSNPs[,1],"hweSNPs"=hweSNPs[["SNP"]],"plateSNPs"=plateSNPs[,2],"batchSNPs"=batchSNPs[,2],"imageSNPs"=imageSNPs,"otherSNPs"=otherSNPs)

    if(!justSNPs) snpList = list("arraySNPs"=arraySNPs,"batchHweSNPs"=batchHweSNPs,"mfSNPs"=mfSNPs,"concordanceSNPs"=concordanceSNPs,"hweSNPs"=hweSNPs,"imageSNPs"=imageSNPs,"plateSNPs"=plateSNPs[,1:2],"batchSNPs"=batchSNPs,"otherSNPs"=otherSNPs)

    return(snpList)
}


# quantile normalise function
quantile.norm <- function(y){
    
    r = rank(y,na.last="keep",ties.method="random")
    R = ( r/max(r,na.rm=TRUE) ) - 0.5/max(r,na.rm=TRUE)
    qn = qnorm(R,0,1)
    return(qn)
    
}

# inverse logit function
inv.logit <- function(x){
    y = exp(x) / (exp(x) + 1)
    return(y)
}

# logit function
logit <- function(x){
    y = log( x/(1-x) )
    return(y)
}


# load UK map files. Downloaded 7/06/2015 using:
## wget http://biogeo.ucdavis.edu/data/gadm2.8/rds/GBR_adm0.rds -O QC-Scripts/R/GBR_adm0.rds
## wget http://biogeo.ucdavis.edu/data/gadm2.8/rds/GBR_adm1.rds -O QC-Scripts/R/GBR_adm1.rds
## wget http://biogeo.ucdavis.edu/data/gadm2.8/rds/GBR_adm2.rds -O QC-Scripts/R/GBR_adm2.rds

## conversion to British National Grid was done on Clare's laptop (needed rgdal libraries) using the script:  /Users/clare/Documents/ClareDPhil/DPhil/UKBiobank_V2/testMirror/QC-Scripts/GWAS/otherGWAS/convertCoordinates.R
                                        #mapFileUK0 = paste0(baseSampleQCDir,"/QC-Scripts/R/GBR_adm0-converted.RData") # most detailed
mapFileUK0 = paste0(baseSampleQCDir,"/QC-Scripts/R/GBR_IRL_IMN_sp_maps.RData")
mapFileUK1 = paste0(baseSampleQCDir,"/QC-Scripts/R/GBR_adm1-converted.RData")
mapFileUK2 = paste0(baseSampleQCDir,"/QC-Scripts/R/GBR_adm2-converted.RData") # just national borders

wgs84 = '+proj=longlat +datum=WGS84'  # coordinate system of the rds polygons
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'  # coordinate system of the east-west coordinates


# create points object for ssplot
make.points <- function(x,y,coord.system=bng,...){
    # other options are: cex, pch, col, bg, etc.
    Points <- SpatialPoints(cbind(x,y),proj4string=CRS(coord.system))
    mappoints = list("sp.points",Points,...)
}

hwe.test <- function(HOM1,HET,HOM0){ # freq is frequency of A1
    freqs = (2*HOM1 + HET)/(2*(HOM1 + HET + HOM0))
    Totals = HOM1 + HET + HOM0
    e1 = Totals*freqs^2
    eh = Totals*2*freqs*(1-freqs)
    e0 = Totals*(1-freqs)^2
    
    qs = (HOM1 - e1)^2/e1 + (HET - eh)^2/eh + (HOM0 - e0)^2/e0
    p = dchisq(qs,df=1)
    
    return(p)
}


# read in GWAS output and apply some QC
signifLevel = 5e-8

read.gwas <- function(dataFile,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=NULL,QCexclusions=c(),extraTitle="",useLmmInf=TRUE,bayesFactorFile=NULL,...){

    print( paste0("reading in GWAS output file for ",chrom) )

    if(is.null(minmaf)) minmaf=0
    if(is.null(mininfo)) minmaf=0
    if(is.null(maxmiss)) minmaf=1

    if(grepl("%%",dataFile)) dataFile = gsub("%%",chrom,dataFile)
    print(dataFile)


    if(grepl("GIANT",dataFile)){
        print("Reading GIANT data...")
        
        if(grepl(".gz$",dataFile)) DFraw = tryCatch(read.table(gzfile(dataFile),sep="",header=TRUE,stringsAsFactors=FALSE,...), error=function(e) NULL) else DFraw = tryCatch(read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE,...), error=function(e) NULL)        
        
        if(is.null(DFraw)) {
            "No data to read!"
            return(NULL)
        }

        if(!is.null(bayesFactorFile)){
            print(bayesFactorFile)
                    
            bfRaw = tryCatch(read.table(bayesFactorFile,sep="",header=TRUE,stringsAsFactors=FALSE,...), error=function(e) NULL)
            if(nrow(DFraw)!=nrow(bfRaw)) print("Cannot combine bayes factors, as files don't have the same number of rows!") else DFraw = cbind(DFraw,bfRaw)            
        }
        
        DFraw = DFraw[!is.na(DFraw$BP),]
        DFraw$BETA = DFraw$b
        
        DF = dplyr::tbl_df(DFraw)

        # add in alternative (unique) ID
        DF$SNP2 = altSNPID(DF$CHR,DF$BP,DF$Allele1,DF$Allele2)
        DF$SNP2.swap = altSNPID(DF$CHR,DF$BP,DF$Allele2,DF$Allele1)

        if(is.null(Ymax)) Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10
        DF$P2=DF$P
        DF$P2[DF$P<(10^-Ymax)] = 10^-Ymax

        Pvalset=paste0("GIANT.chr",chrom,extraTitle)

        return(list("DF"=DF,"Pvalset"=Pvalset))
    }

    # Otherwise for bolt results
    print("Reading BOLT data...")
        
    DFraw = tryCatch(read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE,...), error=function(e) NULL)

    if(is.null(DFraw)) {
        "No data to read!"
        return(NULL)
    }

    if(!is.null(bayesFactorFile)){
        print(bayesFactorFile)
        bfRaw = tryCatch(read.table(bayesFactorFile,sep="",header=TRUE,stringsAsFactors=FALSE,...), error=function(e) NULL)
        if(nrow(DFraw)!=nrow(bfRaw)) print("Cannot combine bayes factors, as files don't have the same number of rows!") else DFraw = cbind(DFraw,bfRaw)            
    }
    
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[(DFraw$A1FREQ > 0.5) & (!is.na(DFraw$A1FREQ))] = 1-DFraw$A1FREQ[(DFraw$A1FREQ > 0.5) & (!is.na(DFraw$A1FREQ))]

    DF = dplyr::tbl_df(DFraw)

    print(head(DF))
    
    if( length(unique(DF$CHR))==1 ) chrom=unique(DF$CHR)
    
    if("INFO"%in%colnames(DFraw)){
        Pvalset = paste(basename(dataFile),".chr",chrom,".maf",minmaf,".info",mininfo,".pruned",extraTitle,sep="")
        if(minmaf!=0) DF = dplyr::filter(DF, MAF > minmaf & INFO > mininfo)
        if(minmaf==0) DF = dplyr::filter(DF, INFO > mininfo) # don't filter out maf==0 snps
        
    } else {
        Pvalset = paste(basename(dataFile),".chr",chrom,".maf",minmaf,".miss",maxmiss,".pruned",extraTitle,sep="")
        if(minmaf!=0) DF = dplyr::filter(DF, MAF > minmaf & F_MISS < maxmiss)
        if(minmaf==0) DF = dplyr::filter(DF, F_MISS < maxmiss) # don't filter out maf==0 snps

    }

    DF = dplyr::filter(DF,!SNP %in% QCexclusions)

    if(("P_BOLT_LMM_INF" %in% colnames(DF)) & (useLmmInf)) {
        print("using BOLT_LMM_INF")
        DF = dplyr::rename(DF, P = P_BOLT_LMM_INF)
    } else {
        print("using P_LINREG")            
        DF = dplyr::rename(DF, P = P_LINREG)
    }

    if(chrom%in%c("X","XY","Y","MT","PAR","PAR1","PAR2")) {

        if(chrom=="PAR") DF$CHR = 25 else DF$CHR = sexChroms[chrom]
        
    } else {
        DF$CHR = as.numeric(DF$CHR)
    }

    
    if(is.null(Ymax)) Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10

    DF$P2=DF$P
    DF$P2[DF$P<(10^-Ymax)] = 10^-Ymax

    # add in alternative (unique) ID
    DF$SNP2 = altSNPID(DF$CHR,DF$BP,DF$ALLELE1,DF$ALLELE0)
    DFraw$SNP2 = altSNPID(DFraw$CHR,DFraw$BP,DFraw$ALLELE1,DFraw$ALLELE0)
    
    return(list("DFraw"=DFraw,"DF"=DF,"Pvalset"=Pvalset))
    
}

get.colClasses <- function(file,columnNames,colTypes=NA){
    header = read.table(file,nrow=1,stringsAsFactors=FALSE)[1,]
    these = match(columnNames,header)
    out = rep("NULL",length(header))
    out[these[!is.na(these)]] = colTypes[!is.na(these)]
    if(sum(is.na(these))>0) print(paste0("No column found with name: ",columnNames[is.na(these)]))
    return(out)
}

# find hit regions using just bp distance
find.hits <- function(DF,chrom,minGap=(10^6)/2){

    hits = (DF$P < 5e-8) & (DF$CHR == chrom)

    DFthese = DF[hits,]
    gaps = diff(DFthese$BP)
    regions = c(which(gaps > minGap)) # treat 1MB gaps as hits
    print(paste0(length(regions)," gaps"))
    starts = c(min(DFthese$BP),DFthese$BP[regions+1])
    ends = c(DFthese$BP[regions],max(DFthese$BP))
    regs = cbind(starts,ends)
                                        # find the top hits in each region
    topVariants = apply(regs,1,function(x){
        p = (DFthese$BP <= x[2]) & (DFthese$BP >= x[1])
        check = DFthese[p,]
        check$SNP[order(check$P)][1]
    })

    return(list("DFthese"=DFthese,"regs"=regs,"topVariants"=topVariants))
}

# find hit regions using recombination rates
find.hits.recomb <- function(DF,chrom,minGap=0.125,buffer=25,recombRates,aveRate = 1.2e-6){
    
    # minGap is maximum cM distance from 'top' snp
    # buffer is extra distance in KB around top snp

    DFchrom = DF[(DF$CHR == chrom),]

    if( !is.null(recombRates) ) {

                                        # linear interpolation of recombrates at snps
        lin = approx(x=recombRates$position, y = recombRates$Genetic_Map.cM., xout=DFchrom$BP)
        DFchrom$Genetic_Map.cM. = lin$y

    } else {
        # assume genetic distance of 0.865 million pb per cM. 
        DFchrom$Genetic_Map.cM. = c(0,diff( DFchrom$BP )*aveRate)  # average rate taken as 3615/3.02GB
    }
    
    # anything with NA get
    NArec = which(is.na(DFchrom$Genetic_Map.cM.))
    if( length(NArec) > 0 ) {
        print(paste0(length(NArec)," snps outside recombination map"))
        for(n in NArec){
            bp = DFchrom$BP[n]
            nearest = recombRates[which(abs(recombRates$position-bp)==min(abs(recombRates$position-bp)))[1],]
            if(bp > nearest$position) rec = nearest$Genetic_Map.cM. + (bp - nearest$position)*nearest$COMBINED_rate.cM.Mb.
            if(bp < nearest$position) rec = 0 # just assume anything beyond the recombination map is 0 cM away
            DFchrom$Genetic_Map.cM.[n] = rec            
        }
    }
    
    hits = (DFchrom$P < signifLevel) & (DFchrom$CHR == chrom)
    
    DFthese = DFchrom[hits,]

    sortedHits = DFthese$SNP[order(DFthese$P)]
    sortedHits2 = DFthese$SNP2[order(DFthese$P)]

    regs = as.data.frame(matrix(NA,nrow=0,ncol=2))
    topVariants = c()
    topVariants2 = c()
    print( paste0( "Thinning ",length(sortedHits)," hits..." ) )

    if(length(sortedHits)>0){
        
        for( i in 1:length(sortedHits) ) {
            s = sortedHits[i]
            s2 = sortedHits2[i]

            bp = DFthese$BP[DFthese$SNP==s][1]
            inRegs = apply(regs,1,function(y) (bp >= y[1]) & (bp <= y[2]) )
            if( sum(inRegs,na.rm=TRUE) > 0 ) next        
            rec = DFthese$Genetic_Map.cM.[DFthese$SNP==s][1]
            minR = rec - minGap
            maxR = rec + minGap
            region = DFchrom$BP[(DFchrom$Genetic_Map.cM. >= minR) & (DFchrom$Genetic_Map.cM. <= maxR)]
            reg = c(min(region,na.rm=TRUE) - buffer*1000,max(region,na.rm=TRUE) + buffer*1000)
            regs = rbind(regs,reg)
            topVariants = c(topVariants,s)
            topVariants2 = c(topVariants2,s2)
        }
        colnames(regs) = c("start","end")
        print( paste0("Found ",nrow(regs)," hit regions after thinning with ",minGap,"cM + ",buffer,"KB") )
        return(list("DFthese"=DFthese,"regs"=regs,"topVariants"=topVariants,"topVariants2"=topVariants2))
        
    } else {
        return(NULL)
    }        
}


# plot region of GWAS. NOTE: needs some objects to exist, e.g catFile, and BB.ps2snp.

plot.gwas.region <- function(DF,topSNP,catFile=NULL,ld=NULL,plotDir,width=NULL,region=NULL,recombRates=NULL,extra="",moreToPlot=FALSE,Ymax=NULL,Pvalset="plot",Ymin=NULL,...){

    print(topSNP)
    #print(head(ld,40))
    #print(nrow(ld))
    #print(head(ld[order(ld$R2,decreasing=TRUE),]))
    
    if((!topSNP%in%BB.ps2snp$AffySNPID)&(!topSNP%in%BL.ps2snp$AffySNPID)) rsid = topSNP else rsid = BB.ps2snp$dbSNPRSID[BB.ps2snp$AffySNPID==topSNP][1]
    
    if((!topSNP%in%BB.ps2snp$AffySNPID)&(topSNP%in%BL.ps2snp$AffySNPID)) rsid = BL.ps2snp$dbSNPRSID[BL.ps2snp$AffySNPID==topSNP][1]
    
    maf = DF$MAF[DF$SNP==topSNP]
        
    pos = DF$BP[DF$SNP==topSNP]
    thisChrom = DF$CHR[DF$SNP==topSNP]
    if(is.null(region)) region = c(pos - width/2,pos + width/2)
            # if region is smaller than specified width then plot width/2 around the SNP
    if(!is.null(width)) region = c( min(region[1],(pos-width/2)),max(region[2],(pos + width/2)) )
    
    these = (DF$BP >= region[1]) & (DF$BP <= region[2]) & (DF$CHR == thisChrom)
        
    if(!is.null(ld)){
        LD = ld[ld$SNP_A==topSNP,]
        ldIndex = match(DF$SNP[these],LD$SNP_B)
        r2 = LD$R2[ldIndex]
        colors = rep("lightgray",sum(these))
        colors[(r2 <= 0.1)&(!is.na(r2))] = "darkgray"
        colors[(r2 > 0.1)&(!is.na(r2))] = colour.scale(r2[(r2 > 0.1)&(!is.na(r2))],colourSet=c("yellow","orange","red"),fixedLims=c(0.1,1))
    } else {
        colors = rep("darkgray",sum(these))
        colors[DF$SNP[these]==topSNP] = "red"
    }
    if(is.null(ld)) print('no ld information for this snp')
        
    shapes = rep(16,sum(these))
    shapes[DF$SNP[these]==topSNP] = 17

    if(is.null(Ymax)) {
        Ymax = min(100,max(-log10(DF$P)))
    }
    DF$P2 = DF$P
    DF$P2[-log10(DF$P)>Ymax] = 10^-Ymax
    
    # is there are known hit in this region??
    if(!is.null(catFile)) known = catFile[(catFile$BP >= region[1]) & (catFile$BP <= region[2]) & (catFile$CHR == thisChrom),] else known=c()

    # if topSNP uses rsid, then change to affy id, as this is what's used to create cluster plots
    topSNPaffy = topSNP
    if(topSNPaffy %in% BL.ps2snp$dbSNPRSID)  topSNPaffy = BL.ps2snp$AffySNPID[BL.ps2snp$dbSNPRSID==topSNPaffy]
    if(topSNPaffy %in% BB.ps2snp$dbSNPRSID)  topSNPaffy = BB.ps2snp$AffySNPID[BB.ps2snp$dbSNPRSID==topSNPaffy]
        
    png(paste(plotDir,"/",Pvalset,"-region-",pos,"-",topSNPaffy,extra,".png",sep=""),width=41,height=12,units="in",res=150)
    par(mar=c(5,5,4.1,5),cex.axis=1.5,cex.lab=1.5)
    plot(NULL,xlim=region,
         xlab="Position",ylab = expression(-log[10](italic(p))),
         ylim=c(0,Ymax),...)
    
    if(!is.null(recombRates)){
        recombMax = max(recombRates$COMBINED_rate.cM.Mb.) # maximum recombination rate on chromsome
        recomPos = (recombRates$position >= region[1]) & (recombRates$position <= region[2])
        xrecom=recombRates$position[recomPos]
        yrecom=recombRates$COMBINED_rate.cM.Mb.[recomPos]
        frac = Ymax/recombMax
        lines(xrecom,yrecom*frac,col="blue")
        axis(4,at = frac*seq(0,recombMax,10),labels=seq(0,recombMax,10))
        mtext("Recombination rate (cM/Mb)", side=4, line=3,cex=par("cex.lab"))
    }
    
    points(DF$BP[these],-log10(DF$P2[these]),col=colors,pch=shapes,cex=2)
    
    if( ( length(known) > 0 ) & (!is.null(catFile)) ) {
        points(known$BP,-log10(known$Pvalue),pch=17,cex=2,col="blue")
        text(known$BP,-log10(known$Pvalue),known$SNP,pos=2,offset=0.5,cex=2)
        # are there any SNPs in UKBiobank that are also in catalogue?
        borderShapes =rep(1,sum(DF$BP%in%known$BP)); borderShapes[DF$SNP[DF$BP%in%known$BP]==topSNP] = 2
        points(DF$BP[DF$BP%in%known$BP],-log10(DF$P2[DF$BP%in%known$BP]),col="black",pch=borderShapes,cex=2,lwd=2)
    }
    
    realPval = ""
    if( DF$P2[DF$SNP==topSNP] != DF$P[DF$SNP==topSNP]) realPval = paste0(" (",DF$P[DF$SNP==topSNP],")")
    text(DF$BP[DF$SNP==topSNP],-log10(DF$P2[DF$SNP==topSNP]),labels=paste0(rsid,realPval),pos=4,offset=0.5,cex=3)
    # allele frequency of top SNP
    text(DF$BP[DF$SNP==topSNP],-log10(DF$P2[DF$SNP==topSNP]),labels=round(maf,4),pos=2,offset=0.5,cex=3)

    abline(h=-log10(5e-8),col="red")
    
    if(!moreToPlot) dev.off()

    return(NULL)
}


 

plot.gwas.region.v2 <- function(DF,topSNP,catFile=NULL,ld=NULL,ldIMP=NULL,plotDir,width=NULL,region=NULL,recombRates=NULL,extra="",moreToPlot=FALSE,Ymax=NULL,Ymin=NULL,Pvalset="plot",showRegion=NULL,addGenes=list(),spaceForGenes=0.3,imputedSNPs=NULL,basicCol="darkgray",imputeCol="lightgray",imputeCex=1.8,extraPoints=NULL,extraPointsCol="blue",extraPointsPch=16,extraPointsCex=1,Legend=NULL,thisChrom=NULL,ldIMPvariant=NULL,newPlot=TRUE,printTopSNP=TRUE,ldColours=c("yellow","orange","red"),pdf=FALSE,resScale=1,axisScale=NULL,xlabExtra="",...){

    print(topSNP)
    #print(head(ld,40))
    #print(nrow(ld))
    #print(head(ld[order(ld$R2,decreasing=TRUE),]))
    
    if((!topSNP%in%BB.ps2snp$AffySNPID)&(!topSNP%in%BL.ps2snp$AffySNPID)) rsid = topSNP else rsid = BB.ps2snp$dbSNPRSID[BB.ps2snp$AffySNPID==topSNP][1]
    
    if((!topSNP%in%BB.ps2snp$AffySNPID)&(topSNP%in%BL.ps2snp$AffySNPID)) rsid = BL.ps2snp$dbSNPRSID[BL.ps2snp$AffySNPID==topSNP][1]

#    rsid = topSNP
    maf = DF$MAF[DF$SNP==topSNP]    
        
    pos = DF$BP[DF$SNP==topSNP]
    if( is.null(thisChrom) ) thisChrom = DF$CHR[DF$SNP==topSNP]    
    
    if(is.null(region)) region = c(pos - width/2,pos + width/2)
            # if region is smaller than specified width then plot width/2 around the SNP
    if(!is.null(width)) region = c( min(region[1],(pos-width/2)),max(region[2],(pos + width/2)) )

    if( region[1]<0 ) region[1] = 0
    
    these = (DF$BP >= region[1]) & (DF$BP <= region[2]) & (DF$CHR == thisChrom)
    
    if(!is.null(imputedSNPs)) {
        impSnps = (imputedSNPs$BP <= region[2]) & (imputedSNPs$BP >= region[1]) & (imputedSNPs$CHR == thisChrom)
        
        if(sum(impSnps) == 0) {
            print("oh no, there are no imputed snps here! ")
            print(thisChrom)
            print(paste0(region[1]," to ",region[2]))
            imputedSNPs=NULL
        } else {
            theseImp = imputedSNPs[impSnps,]
        }

        if(!is.null(ldIMP)){
            print("Plotting imputation data LD")
        }
    }
    
    if(!is.null(ld)){
        LD = ld[ld$SNP_A==topSNP,]
        ldIndex = match(DF$SNP[these],LD$SNP_B)
        r2 = LD$R2[ldIndex]
        colors = rep("gray48",sum(these))
        colors[(r2 <= 0.1)&(!is.na(r2))] = basicCol
        colors[(r2 > 0.1)&(!is.na(r2))] = colour.scale(r2[(r2 > 0.1)&(!is.na(r2))],colourSet=ldColours,fixedLims=c(0.1,1))
    } else {
        colors = rep(basicCol,sum(these))
        colors[DF$SNP[these]==topSNP] = "red"
    }
    if(is.null(ld)) print('no ld information for this snp')
        
    shapes = rep(23,sum(these)) # shape is a diamond for the genotyped snps (or whatever is put in the first argument.)
#    shapes[DF$SNP[these]==topSNP] = 17

    if(is.null(Ymax)) {
        Ymax = min(100,max(-log10(DF$P[these])))
        if(!is.null(imputedSNPs)) Ymax = min(100,max(c(-log10(DF$P[these]),-log10(theseImp$P))))
    }
    if(is.null(Ymin)) {
        Ymin = 0
    }
    
    DF$P2 = DF$P
    DF$P2[-log10(DF$P)>Ymax] = 10^-Ymax
    
    # is there are known hit in this region??
    if(!is.null(catFile)) {
        print(region)
        overlapCat = (catFile$BP >= region[1]) & (catFile$BP <= region[2]) & (catFile$CHR == thisChrom) & (!is.na(catFile$CHR))
        
        known = catFile[overlapCat,]
        if(sum(overlapCat) == 0 ) known = c()
    } else {
        known=c()
    }

    # if topSNP uses rsid, then change to affy id, as this is what's used to create cluster plots
    topSNPaffy = topSNP
    #if(topSNPaffy %in% BL.ps2snp$dbSNPRSID)  topSNPaffy = BL.ps2snp$AffySNPID[BL.ps2snp$dbSNPRSID==topSNPaffy]
    #if(topSNPaffy %in% BB.ps2snp$dbSNPRSID)  topSNPaffy = BB.ps2snp$AffySNPID[BB.ps2snp$dbSNPRSID==topSNPaffy]
        
    if(newPlot) {
        if(pdf) pdf(paste(plotDir,"/",Pvalset,"-region-",pos,"-",topSNPaffy,extra,".pdf",sep=""),width=30,height=12,bg="transparent") else png(paste(plotDir,"/",Pvalset,"-region-",pos,"-",topSNPaffy,extra,".png",sep=""),width=30,height=12,units="in",res=resScale*150,bg="transparent")
    
        par(mar=c(5,5,4.1,5),cex.axis=2.2,cex.lab=2.2)
    }
    
    if( length( addGenes ) > 0 ) layout(mat=c(1,2),heights=c( 1-spaceForGenes,spaceForGenes ) ) # two rows - one for genes
        
#    plot(NULL,xlim=region,xaxt="n",
#         xlab=paste0( "Position on chromosome ",thisChrom,xlabExtra),ylab = expression(-log[10](italic(p))),
#         ylim=c(Ymin,Ymax),...)

    plot(NULL,xlim=region,xaxt="n",
         xlab="",ylab = expression(-log[10](italic(p))),
         ylim=c(Ymin,Ymax),...)
    
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col="lightgray", border=NA)

    if(is.null(axisScale)) axisScale = 1
    axis(1,at = axTicks(1),labels=axTicks(1)/axisScale)

    # plot recombination rates
    if(!is.null(recombRates)){
        recombMax = max(recombRates$COMBINED_rate.cM.Mb.) # maximum recombination rate on chromsome
        recomPos = (recombRates$position >= region[1]) & (recombRates$position <= region[2])
        xrecom=recombRates$position[recomPos]
        yrecom=recombRates$COMBINED_rate.cM.Mb.[recomPos]
        frac = Ymax/recombMax
        lines(xrecom,yrecom*frac,col="blue")
        axis(4,at = frac*seq(0,recombMax,10),labels=seq(0,recombMax,10))
        mtext("Recombination rate (cM/Mb)", side=4, line=3,cex=par("cex.lab"))
    }    

    # significance line
    abline(h=-log10(5e-8),col="red")

    # imputed snps
    if(!is.null(imputedSNPs)) {
        theseImp$P2 = theseImp$P  # fix ymax to match the genotype data
        theseImp$P2[-log10(theseImp$P)>Ymax] = 10^-Ymax
        
        points(theseImp$BP,-log10(theseImp$P2),col=imputeCol,pch=16,cex=imputeCex)

        topImp = theseImp$P==min(theseImp$P)
        
        if( !is.null(ldIMPvariant) ) {
           if(!is.na(ldIMPvariant)) topImp = theseImp$SNP==ldIMPvariant else print("WARNING: the snp given as the top imputed snp is NA.")
        }
        
        if( sum(topImp)>1 ) print("Actually two equal-sized top SNPs")
        if(( !is.null(ldIMP) ) & ( (length(intersect(theseImp$SNP[topImp],ldIMP$SNP_A))>0) ) ) {

                                        # plot LD with top SNP.
            LD = ldIMP[ldIMP$SNP_A%in%theseImp$SNP[topImp],]                
            ldIndex = match(theseImp$SNP,LD$SNP_B)
                
            r2 = LD$R2[ldIndex];
            inLD = !is.na(r2);
            r2 = r2[inLD]; 
            colorsImp = rep("gray48",length(r2))
            colorsImp[(r2 <= 0.1)&(!is.na(r2))] = basicCol
            colorsImp[(r2 > 0.1)&(!is.na(r2))] = colour.scale(r2[(r2 > 0.1)&(!is.na(r2))],colourSet=ldColours,fixedLims=c(0.1,1))
            print("Imputed LD colours:")
            print(table(colorsImp))
            ORDER = order.by.number.occurrences(colorsImp)
            points(theseImp$BP[inLD][ORDER],-log10(theseImp$P2[inLD][ORDER]),col=colorsImp[ORDER],pch=16,cex=imputeCex)
            
        } else {
            print( paste0("No ld information for this imputed snp: ",theseImp$SNP[topImp]))
        }
        
        if(printTopSNP) points(theseImp$BP[topImp],
               -log10(theseImp$P2[topImp]),col="black",pch=21,bg=ldColours[length(ldColours)],cex=3,lwd=2)

    }

    # plot all genotyped snps (as diamonds)
    print(table(colors))
    ORDER = order.by.number.occurrences(colors)

    points(DF$BP[these][ORDER],-log10(DF$P2[these][ORDER]),col=colors[ORDER],bg=colors[ORDER],pch=shapes[ORDER],cex=2)
    
    if(printTopSNP) points(DF$BP[DF$SNP==topSNP],-log10(DF$P2[DF$SNP==topSNP]),col="black",pch=23,bg=ldColours[length(ldColours)],cex=3,lwd=2) # top genotyped snp

    # any extra points
    if(!is.null(extraPoints)){
        
        if(sum(-log10(extraPoints$P2) > Ymax)>0) print("Warning: some of your extraPoints are going to go outside the current axis limit. Check Ymax.")
        points(extraPoints$BP,-log10(extraPoints$P2),col=extraPointsCol,pch=extraPointsPch,cex=extraPointsCex)

    }
    
    if( ( length(known) > 0 ) & (!is.null(catFile)) ) {
        
        points(known$BP,-log10(known$Pvalue),pch=17,cex=2,col="blue")
        text(known$BP,-log10(known$Pvalue),known$SNP,pos=2,offset=0.5,cex=2)

        if(!is.null(imputedSNPs)){
                                        # any known snps in imputed data? Match on BP
            borderShapes =rep(1,sum(theseImp$BP%in%known$BP)); borderCex=rep(2,sum(theseImp$BP%in%known$BP));
            borderCex[theseImp$SNP[theseImp$BP%in%known$BP]==topSNP] = 3
            
            points(theseImp$BP[theseImp$BP%in%known$BP],-log10(theseImp$P2[theseImp$BP%in%known$BP]),col="blue",pch=borderShapes,cex=borderCex,lwd=2)

        }
        
        # are there any SNPs in UKBiobank genotype data that are also in catalogue?
        borderShapes =rep(5,sum(DF$BP%in%known$BP)); borderCex=rep(2,sum(DF$BP%in%known$BP));
        borderCex[DF$SNP[DF$BP%in%known$BP]==topSNP] = 3
        
        points(DF$BP[DF$BP%in%known$BP],-log10(DF$P2[DF$BP%in%known$BP]),col="blue",pch=borderShapes,cex=borderCex,lwd=2)
                                        

    }

                                        # print names of top snps at the top
    if(printTopSNP) {
        
        textY=par("usr")[4]
        textX=par("usr")[2] - (par("usr")[2]-par("usr")[1])/2
                                        #if( DF$P2[DF$SNP==topSNP] != DF$P[DF$SNP==topSNP]) realPval = paste0(" (",DF$P[DF$SNP==topSNP],")")
       # print1 = paste0(rsid,"  (p-value=",format(DF$P[DF$SNP==topSNP],style="scientific",digits=4),";  MAF=",round(maf,4),")")
        print1 = paste0(rsid," ( MAF=",round(maf,4)," )")

        if( !is.null(imputedSNPs) ) {
                                        #   print(topImp)
            realPval2 = ""
                                        #if( theseImp$P[topImp] != theseImp$P2[topImp] ) realPval2 = paste0(" (",theseImp$P[topImp],")")
            
            maf2 = theseImp$MAF[topImp]
            rsid2 = theseImp$SNP[topImp]
            #print2 = paste0(rsid2,"  (p-value=",format(theseImp$P[topImp],style="scientific",digits=4),";  MAF=",round(maf2,4),")")
            print2 = paste0(rsid2," ( MAF=",round(maf2,4)," )")

            PCH = c(23,rep(21,length(rsid2))) 
            
            legend(textX,textY,horiz=TRUE,xjust=0.5,yjust=0.25,bty="n",
                   legend=c(paste0(print1,"   "),print2),col="black",pch=PCH,pt.bg=ldColours[length(ldColours)],cex=2.5,pt.lwd=2,xpd=NA,x.intersp = 0.7)

        } else {
            
            legend(textX,textY,horiz=TRUE,xjust=0.5,yjust=0.25,bty="n",
                   legend=print1,col="black",pch=23,pt.bg=ldColours[length(ldColours)],cex=2.5,pt.lwd=2,xpd=NA,x.intersp = 0.7)
            
        }
    }
    
    ## Add in a legend
    if(!is.null(Legend)) {
        do.call(legend,Legend)
    }

            
    if( !is.null( showRegion )  ) {
        outWidth = diff(showRegion)
        abline(v=showRegion,col="darkgray",lty=3,lwd=3)
        text(mean(showRegion),y=par("usr")[4],labels=paste0(round(outWidth/1000000,1),"Mb region"),xpd=NA,cex=2,pos=1)
    }

    if( length( addGenes ) > 0 ) {
        par( mar=c(2,5,0,5), bty="n") 
                                        # add genes

        plot.genes(chromosome = formatChrom(thisChrom),region=region,local.genes=addGenes$genes,exons=addGenes$exons,label.cex=2,height_in_inches = 1,xaxt="n")        
    }

    if(!moreToPlot) dev.off()

    return(region)
}


##ld.scale <- function(colourSet=c("yellow","orange","red"),fixedLims=c(0.1,1)){
    #colour.scale(r2[(r2 > 0.1)&(!is.na(r2))],colourSet=c("yellow","orange","red"),fixedLims=c(0.1,1))



#}

compute.ld.around.snp <- function(variants,genotypeFile,chrom,samplesFile,bgenSampleFile=NULL,DF){
# NOTE: this calls plink
    rNumber = round(runif(1,1,100000),0)

    if( sum(is.na(variants))>0 ) print("some snps are NA. reducing ld calculations to non-na snps.")
    variants = unique(variants[!is.na(variants)])

    if(length(variants)>0){
        
        write.table(variants,file = paste0("topVariants.",rNumber,".tmp"),quote=FALSE,col.names=FALSE,row.names=FALSE)

                                        # for Y chromosome and MT, just compute ld with everything
        if (chrom %in% c(24,26,"Y","MT","M")){
            
            forSystem = paste0(plink," --bfile ",genotypeFile," --keep-allele-order --keep ",samplesFile," --chr ",chrom," --r2 with-freqs dprime --ld-snp-list topVariants.",rNumber,".tmp --ld-window 100000000 --ld-window-r2 0 --out taggedSNPs.",rNumber,".tmp")
            
        } else {
            
            forSystem = paste0(plink," --bfile ",genotypeFile," --keep-allele-order --keep ",samplesFile," --chr ",chrom," --r2 with-freqs dprime --ld-snp-list topVariants.",rNumber,".tmp --ld-window 1000 --ld-window-r2 0 --out taggedSNPs.",rNumber,".tmp")
            
        }
        
        print(forSystem)
        
        system(forSystem)
        
        ld = read.table(paste0("taggedSNPs.",rNumber,".tmp.ld"),header=TRUE)
        
        system(paste0("rm topVariants.",rNumber,".tmp"))
        system(paste0("rm taggedSNPs.",rNumber,".tmp.ld"))

        return(ld)
    } else {
        print("No non-missing unique snps found for this set of regions.")
        return(NULL)
    }

}


altSNPID <- function(chr,bp,a1,a2,strpad=TRUE){
    if(strpad) out = paste0(str_pad(chr,2,side="left","0"),":",bp,"_",a1,"_",a2)
    if(!strpad) out = paste0(chr,":",bp,"_",a1,"_",a2)
    return(out)
}

compute.ld.around.snp.bgen <- function(variants,bgenFile,bgenSampleFile,chrom,samplesFile,regions=NULL,DF,rNumber=NULL,callthreshold=0.1,sex=NULL){

    print(head(variants))
    
    if(is.null(bgenSampleFile)) {
        print("WARNING: must specify bgen sample file with argument -sample. Otherwise can't compute LD properly.")
        return(NULL)
    }

                                        # first subset bgen file using bgenx. make sure to subset based on the set of samples in samplesFile (e.g whatever set was used in the GWAS)
    convertBgen=FALSE   
    if(is.null(rNumber)) {
        rNumber = round(runif(1,1,100000))
        convertBgen=TRUE # actually get the chunk from the full bgen file and convert to plink. Otherwise just use the existing one!        
    }
    
    if(length(variants)!=length(unique(variants))) print( paste0("Warning: duplicated SNPs in top variants. Not a big problem, but shrinking analysis to ",length(unique(variants))," unique variants."))
    variants2 = unique(variants)
                                        # NOTE: these should be in the form for SNP2 created by altSNPID function.

    if( sum(is.na(variants2))>0 ) print("some snps are NA. reducing ld calculations to non-na snps.")
    variants2 = variants2[!is.na(variants2)]

    if( length(variants2)>0 ){
        
        write.table(variants2,file = paste0("topVariants.",rNumber,".tmp"),quote=FALSE,col.names=FALSE,row.names=FALSE)

        
        DFchrom = DF[DF$CHR==chrom,]
        
        positions = sapply(1:length(variants2),function(s){

            v = variants2[s]
            print(v)
                                        # This won't necessarily match the plink '-ldwindow' argument, so just compute ld within the specified hit region (might contain a heap of SNPs!). If no region argument is supplied, just extract snps within 1MB of the centre SNP.
            sr = which(variants==v)[1]
            
            if(!is.null(regions)){
                reg = regions[sr,]
                starts = as.numeric(reg[1])
                ends = as.numeric(reg[2])
                                        # if region is smaller than 1MB either side of SNP, then take that region.

                                        # if( starts > ( DF$BP[DF$SNP==v] - 20e3 ) ) starts = DF$BP[DF$SNP==v] - 20e3
                                        # if( ends < ( DF$BP[DF$SNP==v] + 20e3 ) ) ends = DF$BP[DF$SNP==v] + 20e3
                
            } else {
                starts = DF$BP[DF$SNP2==v] - 20e3; starts[starts<0]=0 
                ends = DF$BP[DF$SNP2==v] + 20e3
            }

            
#######  # TEMP FOR TESTING
                                        #starts = DF$BP[DF$SNP==v] - 10e3; starts[starts<0]=0 
                                        #ends = DF$BP[DF$SNP==v] + 10e3
########
            print(starts)
            print(ends)
            if(starts<0) starts=0
            
            nSNPs = sum( (DFchrom$BP<=ends)&(DFchrom$BP>=starts) )
            if( nSNPs > 10000 ) print( paste0( nSNPs ," SNPs in chunk. May take a long time to process...") )
            
            chrs = sprintf("%02d",DF$CHR[DF$SNP2==v])
            if(chrs%in%c(23,24,26)) chrs=names(sexChroms)[sexChroms==chrs]
            if(chrs==25) chrs="PAR"

############
            ldThisSNP = ld.bgen(snp,chrs,starts,ends,rNumber,convertBgen=convertBgen,callthreshold,sex)
            
############
            
        },simplify=FALSE)
        
        #system(paste0("rm topVariants.",rNumber,".tmp"))
        system(paste0("rm taggedSNPs.",rNumber,".tmp.ld"))
                                        #system(paste0("rm taggedSNPs.",rNumber,".tmp.ld.log"))
        
        ld = rbind_all(positions)

                                        # revert back to original SNP IDs?
        linkSNP = DF$SNP; names(linkSNP)=DF$SNP2
        ld$SNP2_A = ld$SNP_A
        ld$SNP2_B = ld$SNP_B
        
        ld$SNP_A = linkSNP[ld$SNP2_A]
        ld$SNP_B = linkSNP[ld$SNP2_B]
        
        return(ld)
    } else {
        print("No non-missing unique snps found for this set of regions.")
        return(NULL)
    }
}

##### This function uses plink to do the ld calculation
ld.bgen <- function(snp,chrs,starts,ends,rNumber,convertBgen=TRUE,callthreshold,sex=NULL,extraOption=""){

    if(convertBgen){
        bgenFile1 = bgenFile
        
        if(chrs[1]=="PAR"){
            if(starts < 2699509) {
                chrs="PAR1"
                bgenFile1 = gsub("PAR2","PAR1",bgenFile)
            }
            if(starts > 2699507) {
                chrs="PAR2"
                bgenFile1 = gsub("PAR1","PAR2",bgenFile)
            }
        }
        
        positions = paste0(chrs,":",starts,"-",ends)
        print(chrs)
                                        #run bgenx to get subset bgen file.
        forSystem1 = paste0(bgenix," ",bgenFile1," ",positions," > genotypes.",rNumber,".tmp.bgen")
        print(forSystem1)
        system(forSystem1)

                                        # The above creates a v1.1 bgen file with sample names in the header. Need to convert to something that plink can read!
        
                                        #forSystem2 = paste0(qctool," -g genotypes.",rNumber,".tmp.bgen -og genotypes.",rNumber,".tmp2.bgen")
                                        #print(forSystem2)
                                        #system(forSystem2)

                                        # ~1,000 variants in a typical region for imputed data!? Takes about 20secs to chunk up.
                                        #tempBgen = paste0("genotypes.",rNumber,".tmp2.bgen")
        tempBgen = paste0("genotypes.",rNumber,".tmp.bgen")
        
                                        # first get list of duplicate variants (if any!)
        
                                        #    forSystem3 = paste0( plink," --bgen ",tempBgen," --sample ",bgenSampleFile," --list-duplicate-vars suppress-first --keep-allele-order --out duplicates.",rNumber,".tmp")
                                        #    system(forSystem3)

                                        #    forSystem4 = paste0("awk '{print $4}' duplicates.",rNumber,".tmp.dupvar | tail -n +2 > duplicates.",rNumber,".tmp.dupvar.pos")
                                        #    system(forSystem4)

                                        # update bim file so rsids are all unique, and match the top.variants list.

        print(paste0("Converting bgen probabilities to plink format using hard-call-threshold of ",callthreshold))
        if(!is.null(sex)) {
            if(sex=="M") extraOption = " --filter-males"
            if(sex=="F") extraOption = " --filter-females"
        }
            
        forSystem4 = paste0( plink," --bgen ",tempBgen," --hard-call-threshold ",callthreshold," --sample ",bgenSampleFile," --allow-extra-chr --keep-allele-order",extraOption," --make-bed --out genotypes.",rNumber,".tmp")
        print(forSystem4)
        system(forSystem4)
        
                                        # make dummy IDs for snp names based on position etc. THIS SHOULD MATCH THE FUNCTION altSNPID()
        system(paste0("mv genotypes.",rNumber,".tmp.bim genotypes.",rNumber,".tmp.bim1"))
        chromName = chrs[1]
        if(chrs[1]%in%c("X")) chromName=23
        if(chrs[1]%in%c("PAR","PAR1","PAR2")) chromName=25
        
        system(paste0("awk 'OFS=\"\t\" {print \"",chromName,"\",\"",chromName,":\"$4\"_\"$5\"_\"$6,$3,$4,$5,$6}' genotypes.",rNumber,".tmp.bim1 > genotypes.",rNumber,".tmp.bim" ))
        
                                        # NOTE: with genotype data we used --ld-window 1000 which means compute LD with only 1000 markers either side of main SNP. SNPs are much denser on imputed data, so will use a kb distance instead. 1000 markers in genotype data ~20000 BP = 20KB.
                                        # test = read.table(/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.bim)
                                        # markers = test[,4]
                                        # diffs1K = sapply(1:(length(markers)-1000), function(i) markers[i+1000]-markers[i])
                                        # mean(diffs1K)
                                        # call plink to compute LD

                                        #    forSystem5 = paste0(plink," --bgen ",tempBgen," --sample ",bgenSampleFile," --keep-allele-order --exclude duplicates.",rNumber,".tmp.dupvar.pos --keep ",samplesFile," --chr ",chrom," --r2 with-freqs dprime --ld-snp-list topVariants.",rNumber,".tmp --ld-window 100000 --ld-window-kb 500 --ld-window-r2 0.1 --out taggedSNPs.",rNumber,".tmp")
    }
    
                                        #forSystem5 = paste0(plink," --bfile genotypes.",rNumber,".tmp --keep-allele-order --keep ",samplesFile," --r2 with-freqs dprime --ld-snp-list topVariants.",rNumber,".tmp --ld-window 100000 --ld-window-kb 500 --ld-window-r2 0.1 --out taggedSNPs.",rNumber,".tmp")
# NEW: 11/10/2017.  No dprime; removed window restriction.
    forSystem5 = paste0(plink," --bfile genotypes.",rNumber,".tmp --keep-allele-order --keep ",samplesFile," --r2 --ld-window 10000000 --ld-window-kb 2000 --ld-snp-list topVariants.",rNumber,".tmp --ld-window-r2 0.1 --out taggedSNPs.",rNumber,".tmp")
    system(forSystem5)           
                                        # read back the LD results for this SNP
    
    ldThisSNP = read.table(paste0("taggedSNPs.",rNumber,".tmp.ld"),header=TRUE)

    return(ldThisSNP)    
    
}

##### This function uses ldscore to do the ld calculation http://www.christianbenner.com/
ld.bgen2 <- function(snp,chrs,starts,ends,rNumber,convertBgen=TRUE,callthreshold){
    
    if(chrs%in%c(23:26)) chrs = names(sexChroms)[sexChroms==chrs]
                                        # get only this chromosome.
    positions = paste0(chrs,":",starts,"-",ends)
                                            #run bgenx to get subset bgen file.
    forSystem1 = paste0(bgenix," ",bgenFile," ",positions," > genotypes.",rNumber,".tmp.bgen")
    print(forSystem1)
    system(forSystem1)

                                        #run bgenx to get subset bgen file.
    forSystem2 = paste0(ldstore," --bgen genotypes.",rNumber,".tmp.bgen --bcor taggedSNPs.",rNumber,".tmp.bcor --n-threads 8 --ld-thold 0.1")
    print(forSystem2)
    system(forSystem2)

    forSystem2 = paste0(ldstore," --bcor taggedSNPs.",rNumber,".tmp.bcor --merge 8 --table taggedSNPs.",rNumber,".tmp.ld")
    print(forSystem2)    
    system(forSystem2)
    
}

    
                                        # read in imputed snps list
inputedDataPrefix='/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr%%.hrc+uk10k_sorted_8bit_rsids'

imputedSNPLists=paste0(inputedDataPrefix,'.list')
read.imputed.snps <- function(chroms,imputedSNPFile=imputedSNPLists){
    imputedSNPs = list()
    print(paste0(length(chroms)," chroms requested. Might take a while. ~2mins each for large chromosomes."))
    for(i in chroms){
        r = read.table(gsub("%%",i,imputedSNPFile),header=TRUE,stringsAsFactors=FALSE)
        r$SNP2 = altSNPID(r$chromosome,r$position,r$first_allele,r$alternative_alleles,strpad=TRUE)
        imputedSNPs[[as.character(i)]] = r
    }
    return(imputedSNPs)
}


# read in phased file?
phaseFile = "/well/ukbiobank/expt/V2_QCed.export/data/imputation_pipeline_input/v4/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.phasing_QC.chr"

read.phased.snps <- function(chroms,phasedSNPFile=phaseFile){
    phasedBims = list()
    if(chroms=="genome") chroms=c(1:22,"X","XY")
    
    for(i in chroms){
                                        # NOTE: sex chromosomes not there yet!
        print(i)
        r = read.table(paste0(phasedSNPFile,i,".bim"),header=FALSE,stringsAsFactors=FALSE)
        r$SNP2 = altSNPID(r$V1,r$V4,r$V5,r$V6)  # reconstruct ID in imputed data. Assume that alleles are in the same order!
        if(i%in%names(sexChroms)) i =  as.character(sexChroms[i])
        phasedBims[[as.character(i)]] = r

    }
    return(phasedBims)    
}


imputedMAFLists='/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/chr%%.maf+info'
read.imputed.maf <- function(chroms,imputedMAFFile=imputedMAFLists){
    
    imputedMAF = list()
    print(paste0(length(chroms)," chroms requested. Might take a while. ~1min each for large chromosomes."))
    for(i in chroms){
        r = read.table(gsub("%%",i,imputedMAFFile),colClasses=c("character","numeric","character","character","numeric","numeric"),header=FALSE,stringsAsFactors=FALSE,comment.char = "#")
        r$SNP2 = altSNPID(i,r$V2,r$V3,r$V4)
        imputedMAF[[as.character(i)]] = r
    }
    return(imputedMAF)
}

getMafBins <- function(mafs,mafBins,incZero=FALSE){
    bins=rep(NA,length(mafs))
    for(b in 1:(length(mafBins)-1)){
        bins[(!is.na(mafs))&(mafs >= mafBins[b])&(mafs < mafBins[b+1])]=b
    }
    if(incZero) bins[mafs==0] = 0

    return(bins)
}

snpFrequencyFiles= "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr"
# frequency calculations done using /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/InterimComparison/compute-frequencies.sh
read.genotyped.maf <- function(snpFrequencyData=snpFrequencyFiles,mafBins = c(0,1/1000,1/100,5/100,Inf) ,incZero=FALSE){

    snpFreqFiles = list.files(path=dirname(snpFrequencyData),pattern=paste0(basename(snpFrequencyData),".*frqx"),full.names=TRUE)

    snpFreqs= sapply(snpFreqFiles,function(x) read.delim(x,header=TRUE,stringsAsFactors=FALSE,sep="\t"),simplify=FALSE)

    freqs=rbind_all(snpFreqs)

     # compute frequency of A1 allele.

                                        # valid for all snps
    freqs$acA1 = (2*freqs$C.HOM.A1. + freqs$C.HET. + freqs$C.HAP.A1.)
    freqs$acA2 = (2*freqs$C.HOM.A2. + freqs$C.HET. + freqs$C.HAP.A2.)
    
    freqs$freq = (2*freqs$C.HOM.A1. + freqs$C.HET. + freqs$C.HAP.A1.)/(2*(freqs$C.HOM.A1. + freqs$C.HET. + freqs$C.HOM.A2.) + freqs$C.HAP.A1. + freqs$C.HAP.A2.)

    print( paste0( sum(is.na(freqs$freq) )," markers with missing freq." ) )
    
    freqs$MAF = freqs$freq
    freqs$MAF[(freqs$freq>0.5)&(!is.na(freqs$freq))] = 1 - freqs$freq[(freqs$freq>0.5)&(!is.na(freqs$freq))]

    
    freqs$MinorAllele = 1
    freqs$MinorAllele[freqs$freq>0.5] = 2
    
    # minor allele counts
    freqs$MAC = freqs$acA1
    freqs$MAC[freqs$MinorAllele==2] = freqs$acA2[freqs$MinorAllele==2]
    
    freqs$bin=NA
    
    for(b in 1:(length(mafBins)-1)){
        freqs$bin[(freqs$MAF >= mafBins[b])&(freqs$MAF < mafBins[b+1])]=b
    }
    if(incZero) freqs$bin[freqs$MAF==0] = 0

    return(freqs)
}

read.genotyped.snps <- function(chroms,genotypedSNPFile=GenotypesForRelease){
    genoBims = list()
    if(chroms=="genome") chroms=c(1:22,names(sexChroms))
    
    for(i in chroms){
                                        # NOTE: sex chromosomes not there yet!
        print(i)
        r = read.table(paste0(genotypedSNPFile,i,".bim"),header=FALSE,stringsAsFactors=FALSE)
        r$SNP2 = altSNPID(r$V1,r$V4,r$V5,r$V6)  # reconstruct ID in imputed data. Assume that alleles are in the same order!
        colnames(r) = c("CHR","SNP","V3","BP","A1","A2","SNP2")
        genoBims[[as.character(i)]] = r
    }
    return(genoBims)    
}


###############
## Plot manhattans!!
source('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/myManhattan.R')

plot.BOLT.pvalues <- function(Data,chrom,plotOutDir,plotQQ=TRUE,catFile = NULL,Ymax=NULL,extraData=NULL,extraPointsCol="lightgray",moreToPlot=FALSE,extraFilename="",qqylim=NULL,cexManhattan=1,...) {

    if(is.null(Data)){
        print("No data for this chromosome, I'm afraid.")
        return(NULL)
    }
    
    # Data is output from read.gwas() function
    DF = Data$DF
    if(chrom!="genome") DF = filter(DF,CHR %in% chrom)

    if(nrow(DF)==0){
        print("No data for this chromosome, I'm afraid.")
        return(NULL)
    }
        
    Pvalset =  gsub("genome",chrom,Data$Pvalset)
    print(Pvalset)
        
    maxP = round(max(-log10(DF$P),na.rm=T))
    print(paste('max -log(pval) = ',maxP))

    if(is.null(Ymax)) Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10



    
    
########### qq plot p-values
    if(plotQQ){

                                        # compute Lambda
        chisq1 = qchisq(DF$P,1,lower.tail=FALSE)
        lambda1 = median(chisq1)/qchisq(0.5,1)    # careful with LMM...
        
        png(paste(plotOutDir,"/",Pvalset,extraFilename,"-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
        DF$P2=DF$P
        DF$P2[DF$P<(10^-Ymax)] = 10^-Ymax
        qqman::qq(DF$P2)
        mtext(3,text = paste0(expression(lambda)," = ",round(lambda1,2)))
        dev.off()
    }


    
########### plot p-values manhattan

    png(paste(plotOutDir,"/",Pvalset,extraFilename,"-manhattan%02d.png",sep=""),width=61,height=12,units="in",res=150,bg="transparent")
    par(las=1,font.main=1,cex.axis=2,cex.lab=2,mar=c(7 ,7, 5 ,2))
    myManhattan(DF,ymax=Ymax,suggestiveline = FALSE,xpd=NA,cex=cexManhattan,...)

    # are we plotting any extra points?
    if(!is.null(extraData)){
        extraData$P2 = extraData$P
        extraData$P2[-log10(extraData$P) > Ymax] = 10^-Ymax
        print( sum(-log10(extraData$P) > Ymax) ) 
        
        points(extraData$BP,-log10(extraData$P2),col=extraPointsCol)
    }
#    myManhattan(DF,ymax=Ymax,suggestiveline = FALSE,xpd=NA,cex=1,col="transparent")

    # add extra catalogue hits?
    if(!is.null(catFile) & (chrom!="genome")) {
        print( "Printing catalogue hits..." )

        catFileSub = catFile[catFile$CHR %in% chrom,]
        print( paste0( sum(-log10(catFileSub$Pvalue)>Ymax)," catalogue hits above ",Ymax) )
        
        catFileSub$Pvalue[-log10(catFileSub$Pvalue)>Ymax] = 10^(-Ymax)
        # plot non-European hits differently
        colors = rep("red",dim(catFileSub)[1])
#        colors[!grepl("European",catFileSub$Ancestry)] = "blue"
        print( table(catFileSub$Ancestry[!grepl("European",catFileSub$Ancestry)]) )

        # do we have these SNPs in UKBiobank? match on chrom and position
        #inHere = (catFileSub$BP %in% DF$BP)&(catFileSub$CHR == DF$CHR)
        #catFileSub$Pvalue[inHere] = DF$Pvalue[]
        
        points(catFileSub$BP/1000000,-log10(catFileSub$Pvalue),pch=8,col=colors,cex=4,lwd=1,xpd=NA)
        points(catFileSub$BP/1000000,-log10(catFileSub$Pvalue),pch=16,col=colors,cex=2,xpd=NA)
        
    }
    
    if(( !moreToPlot ) ) dev.off()


    
########## plot effect sizes
    DF$index = DF$BP

    if( length(unique(DF$CHR)) > 1 ){
        for(i in unique(DF$CHR)){        
            if(i>1) DF$index[DF$CHR==i] = DF$index[DF$CHR==i] + max(DF$BP[DF$CHR==(i - 1)])
        }
    }        
    snps = which((DF$P < 5e-8)&(!is.na(DF$P)))
    
    if("BETA"%in%colnames(DF) & ( length(snps) > 0 )){
        
        beta = DF
        beta$BETA[DF$A1FREQ > 0.5] = -beta$BETA[DF$A1FREQ > 0.5]
        beta = beta[snps,]               
        
#        png(paste(plotOutDir,"/",Pvalset,"-EffectSizes.png",sep=""),height=1000,width=1000,res=150)
#        myManhattan(beta,p="BETA",logtransform=FALSE,genomewideline=0,suggestiveline=FALSE)
        
#        dev.off()
        
    }
}

plot.NULL.qqplot <- function(Data,chrom,plotOutDir,Ymax=NULL,extraData=NULL,moreToPlot=FALSE,extraFilename="",...){

    
    if(is.null(Data)){
        print("No data for this chromosome, I'm afraid.")
        return(NULL)
    }
    
    # Data is output from read.gwas() function
    DF = Data$DF
    if(chrom!="genome") DF = filter(DF,CHR %in% chrom)

    if(nrow(DF)==0){
        print("No data for this chromosome, I'm afraid.")
        return(NULL)
    }
        
    Pvalset =  gsub("genome",chrom,Data$Pvalset)
    print(Pvalset)
        
    maxP = round(max(-log10(DF$P),na.rm=T))
    print(paste('max -log(pval) = ',maxP))

    if(is.null(Ymax)) Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10



}


myQQplot <- function(pvector,...){
    # THIS function is just the same as qqman:qq, but just outputs the x/y values.
    if (!is.numeric(pvector)) 
        stop("Input must be numeric.")
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
        is.finite(pvector) & pvector < 1 & pvector > 0]
    o = -log10(sort(pvector, decreasing = FALSE))
    e = -log10(ppoints(length(pvector)))
    return(list("o"=o,"e"=e))
}




combine.cluster.regions.plots <- function(snp,prefix,filename=NULL,clusterPlotDir="./"){
    
    print(paste0(prefix,"-region.*",snp,".png"))
    plots1 = list.files(path="./",pattern=paste0(prefix,"-region.*",snp,".png"),recursive=TRUE)
    plots2 = list.files(path=clusterPlotDir,pattern=paste0("clusterplot.*",snp),recursive=TRUE)

    continue = TRUE
    plotRegionPlot = TRUE
    
    if(length(plots1)==0) {
        print(snp)   
        print("Could not find a region plot for this snp!")
                                        # some region plots might have rsids instead of affy ids. then search on position and chromosome
        rsid=snp
        if((!snp%in%BB.ps2snp$AffySNPID)&(!snp%in%BL.ps2snp$AffySNPID)) rsid = snp else rsid = BB.ps2snp$dbSNPRSID[BB.ps2snp$AffySNPID==snp][1]
        
        if((!snp%in%BB.ps2snp$AffySNPID)&(snp%in%BL.ps2snp$AffySNPID)) rsid = BL.ps2snp$dbSNPRSID[BL.ps2snp$AffySNPID==snp][1]

        if(rsid==snp){
                                        # still nothing, look in sex chromosomes
            if((!snp%in%BB.ps2snpSex$AffySNPID)&(!snp%in%BL.ps2snpSex$AffySNPID)) rsid = snp else rsid = BB.ps2snpSex$dbSNPRSID[BB.ps2snpSex$AffySNPID==snp][1]
            
            if((!snp%in%BB.ps2snpSex$AffySNPID)&(snp%in%BL.ps2snpSex$AffySNPID)) rsid = BL.ps2snpSex$dbSNPRSID[BL.ps2snpSex$AffySNPID==snp][1]
        }

        print(rsid)
        
        plots1 = list.files(path="./",pattern=paste0(prefix,"-region.*",rsid,".png"),recursive=TRUE)
        
        if(length(plots1)>0){
            print("Yay, we've now found one with rsid!")
            continue=TRUE
        } else {
            print("still could not find a region plot even after looking for rsid...")
            if(length(plots2)==0) continue=FALSE else plotRegionPlot = FALSE
        }
    }
    if(length(plots2)==0) {
        print(snp)   
        print("Could not find a cluster plot for this snp!")
        continue=FALSE
    }
        
    if(length(plots1) > 1) {
        print(paste0("Two versions of this region plot available. Using just the latest one. "))
        print(plots1)
        plots1=plots1[2]
    }

    if(continue){
        if(is.null(filename)){
            if(length(plots1)!=0) filename = gsub(".png","-withClusterPlots.pdf",plots1) else filename = gsub(".png","-withClusterPlots.pdf",plots2[1])
            
        }

        print(filename)
        if(length(plots2) < 12) print(paste0(length(plots2)," cluster plots found for this snp."))
        
        pdf(filename,width=10,height=10)
        
        par(mar=c(0,0,0,0))
        plot.new()
        
        if(plotRegionPlot) {
            img <- readPNG(plots1)
                                        # plot top row, then 3 x 4 matrix below
            rasterImage(img,0,3/4,1,1)
        }
        
        p = 0
        for(i in plots2[1:min(length(plots2),12)]) {  # print just the first 12 files
            img <- readPNG(paste0(clusterPlotDir,"/",i))
            p = p+1
            r = ceiling(p/4); c = p%%4; if(c==0) c = 4
            rasterImage(img,(c-1)/4,(3-r)/4,c/4,(4-r)/4)
        }

        dev.off()
    }
}




##################
# Functions for cnv extraction.

# read data for each batch in sex chromosomes. Must have run submit-extract-sex-CNV.sh for this to work.

extract.cnv <- function(inds=NULL,batch,type,snpsToInclude=NULL,ps2snp=ps2snp){
    
# most of the arguments are actually redundant...
    cnvTextFile = paste0( baseSampleQCDir,"/data/CNV/",batch,"/AxiomGT1.cnv.",type,".final.sexchrom.txt" )
    
    out = read.table(cnvTextFile,header=TRUE,stringsAsFactors=FALSE)

    # fix best array names
    oldNames = colnames(out)[-c(1:3)]
    newNames = gsub("\\.","-",toupper(gsub("\\.CEL","",oldNames)))
    
    return(out)

}

# extract CNV data for a certain set of individuals in a batch (or all of them if inds=NULL)

read.cnv.hdf5 <- function(inds,snps,batch=NULL,type,otherInfo=NULL,batchCheck=TRUE){

    if( (is.null(otherInfo)) & (is.null(batch)) ) {
        print("Can't extract. Must specify otherInfo or batch")
        return(NULL)
    }
                                        # snps must be in probesetID format
    if( (is.null(batch)) & (!is.null(otherInfo)) ){

        if(is.null(inds)) {
            print("Are you sure you want to extract ALL the cnv data??")
            theseBatches = unique(otherInfo$Batch)
        } else {
            theseBatches = unique(otherInfo$Batch[otherInfo$PIID %in% inds])
        }
        
    }
    
    if(!is.null(batch)) theseBatches=batch
    
    h5 = paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095.V2_All.",type,".sexchrom.h5")
    print("Reading hdf5 file")    
    batchesInH5 = unique(sapply(h5ls(h5)$name,function(s) str_split(s,pattern="\\.")[[1]][1]))
    missingBatches = theseBatches[!theseBatches%in%batchesInH5]
    if(length(missingBatches) > 0 ) {
        print( paste0("WARNING: the following batches are not in the hdf5 file. ",paste(missingBatches,collapse=",")) )
        theseBatches = theseBatches[theseBatches%in%batchesInH5]
    }
        
    dataBatches = sapply(theseBatches,function(batch){
        print(batch)
        snpOrder = h5read(h5,paste0(batch,".snpOrder"))
        sampleOrder = h5read(h5,paste0(batch,".sampleOrder"))
        
        if(!is.null(snps)) indexSnp = which(snpOrder$probeset_id%in%snps) else indexSnp=c(1:length(snpOrder))
        if(!is.null(inds)) indexSamples = which(sampleOrder%in%inds) else indexSamples=c(1:length(sampleOrder))

        if((!is.null(inds)) & ( sum(inds%in%sampleOrder) == 0 ) ) {
            print("No relevant samples in this batch.")
            dataSubset=NULL
        } else {
            dataSubset = h5read(h5,paste0(batch,".data"),index=list(indexSnp,indexSamples))
            rownames(dataSubset) = snpOrder$probeset_id[indexSnp]
            colnames(dataSubset) = sampleOrder[indexSamples]
        }
        return(dataSubset)
    },simplify=FALSE)
    
#    out = Reduce(cbind,dataBatches)     # this doesn't work with all batches. Error in f(init, x[[i]]) : long vectors not supported yet: bind.c:1304! Just output dataBatches
    out = dataBatches     # this doesn't work with all batches! Just output dataBatches
    return(out)
}


is.in.ellipse <- function(x1,y1,x0,vx,y0,vy,cov,c=2){  
  theta = 0.5 * atan2(2*cov, vx-vy)
  sint = sin(theta)
  cost = cos(theta)
  a = c*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
  b = c*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
  sint = sin(theta)
  cost = cos(theta)
  test = ((cost*(x1 - x0) + sint*(y1 - y0))^2)/a^2 + ((sint*(x1 - x0) - cost*(y1 - y0))^2)/b^2
  return(test <= 1)
}


#########
# functions for plotting cytobands

#cytoPlot <- function(data,chrom)

#rect()




#################
# Function to apply linear model for each PC with country of birth as covariates
# Used in numbersForPaper.R, and based on the script in pca-ethnic-background-correlation.R

cobByPClm <- function(PCs,cob2,eth2){
    
    nPCs = sum(grepl("^PC",colnames(PCs)))

    print(paste0("Running linear regression on ",nPCs," PCs for country of birth..."))

    cobLMnorm = sapply(1:nPCs,function(pc){
        print(pc)
        y = PCs[[paste0("PC",pc)]]
        yNorm = y/sd(y)  # normalise by standard deviation
                                        #    LMfit = lm(y~0+eth2+cob)
        LMfit = lm(yNorm~0+cob2)   
        s = summary(LMfit)
        return(s)
    },simplify=FALSE)


    cobLMraw = sapply(1:nPCs,function(pc){
        print(pc)
        y = PCs[[paste0("PC",pc)]]
                                        #yNorm = y/sd(y)  # normalise by standard deviation and mean
                                        #    LMfit = lm(y~0+eth2+cob)
        LMfit = lm(y~0+cob2)   
        s = summary(LMfit)
        return(s)
    },simplify=FALSE)

    uCob2 = unique(cob2);uCob2 = uCob2[!is.na(uCob2)] 
    cobLMnormMean1 = sapply(1:nPCs,function(pc){
        print(pc)
        y = PCs[[paste0("PC",pc)]]
        yNorm = y/sd(y)  # normalise by standard deviation
        o = sapply(uCob2,function(u) mean(yNorm[cob2==u],na.rm=TRUE))
        return(o)
    },simplify=TRUE)


    coeffMat = sapply(cobLMnorm,function(s) s$coefficients[,1])
    pvals = sapply(cobLMnorm,function(s) s$coefficients[,4])

    signif = 0.05/prod(dim(coeffMat))
    coeffMat2 = coeffMat
    coeffMat2[pvals>=signif] = 0

                                        # raw coefficients
    coeffMatraw = sapply(cobLMraw,function(s) s$coefficients[,1])
    pvalsRaw = sapply(cobLMraw,function(s) s$coefficients[,4])
    coeffMatraw2 = coeffMatraw
    coeffMatraw2[pvalsRaw>=signif] = 0

                                        # just the means
    rownames(cobLMnormMean1) = uCob2
    cobLMnormMean = cobLMnormMean1[match(gsub("cob2","",rownames(coeffMat)),rownames(cobLMnormMean1)),]
    cobLMnormMean2 = cobLMnormMean
    cobLMnormMean2[pvals>=signif] = 0


    rownames(coeffMat) = rownames(coeffMat2) = rownames(pvals) = rownames(coeffMatraw) = rownames(coeffMatraw2) = rownames(pvalsRaw) = gsub("cob2","",rownames(coeffMat))
    colnames(coeffMat) = colnames(coeffMat2) = colnames(pvals) = colnames(coeffMatraw) = colnames(pvalsRaw) = paste0("PC",1:nPCs)


#################
                                        # Plot the matrix of coefficients
    cobCols = getColoursDistant3(length(table(cob2)))
    names(cobCols) = names(sort(table(cob2)))
    cobChars = rep(c(1:5),100)[1:length(cobCols)]
    names(cobChars) = names(cobCols)

                                        # Make colours based on most common ethnicity value
    maxEthnic = t(sapply(names(sort(table(cob2))), function(country){
        if(country=="Other/Unknown") return(list("Other/Unknown",1))
        tab = table(eth2[cob2==country])
        tab = tab[!is.na(tab)]
        #tab = tab[names(tab)!="Other/Unknown"]
        if(length(tab)==0) print(country)
        if(length(tab)>1) dist = (max(tab) - max(tab[tab!=max(tab)]))/sum(tab) else dist = 1
                                        #if(dist < 0.5) print(country)
        maxEth = names(tab)[tab==max(tab)][1]
                                        # }
        return(list(maxEth,dist))
    },simplify=TRUE))

    cobCols = ethnicity2col[unlist(maxEthnic[,1])]
    cobChars = ethnicity2char[unlist(maxEthnic[,1])]
    names(cobCols) = names(cobChars) = rownames(maxEthnic)


    Legend = cbind.data.frame(unlist(maxEthnic[,1]),cobCols,cobChars[names(cobCols)])  # make legend
    colnames(Legend) = c("Factor","colour","shape")
    rownames(Legend) = names(cobCols)


    treeOrder = hclust(dist(coeffMat)) # order by normalised coefficients
    treeOrder = hclust(dist(coeffMat2)) # order by normalised coefficients after exluding non-significant elements

    return(list(cobLMnormMean2,cobLMnormMean,coeffMat2,Legend,maxEthnic,signif,pvals,treeOrder))
    
}




###########
# Some functions for plotting tables

italic <- function(x){
paste0('{\\emph{', x, '}}')
}
large <- function(x){
paste0('{\\Large ', x, '}')
}
bold <- function(x){
paste0('{\\bfseries ', x, '}')
}
white <- function(x){
    paste0('\\textcolor{White}{', x,'}')    # NOTE: this won't work unless you've defined "White" in your latex preamble.    
}

#########
# plot rectangle on hexbin plot

myHexbin.rect <- function(hvp,x0,y0,x1,y1,
    col = "black", lty = 1, lwd = 2, ...){
  
    pushHexport(hvp, clip = "off")
    xx <- current.viewport()$xscale
    yy <- current.viewport()$yscale
    grid.lines(x=c(x0,x1),y = c(y0,y0), default.units = "native", gp = gpar(col = col, lty = lty, lwd = lwd))
    grid.lines(x=c(x1,x1),y = c(y0,y1), default.units = "native", gp = gpar(col = col, lty = lty, lwd = lwd))
    grid.lines(x=c(x1,x0),y = c(y1,y1), default.units = "native", gp = gpar(col = col, lty = lty, lwd = lwd))
    grid.lines(x=c(x0,x0),y = c(y1,y0), default.units = "native", gp = gpar(col = col, lty = lty, lwd = lwd))
    
    popViewport()
}               

#########
                                        # plot x and y labels on hexbin

myHexbin.xylab <- function(hvp,xlab,ylab,lineHeight=(-2),gp = gpar(fontsize = 16),...){
    pushHexport(hvp, clip = "off")
    grid.text(xlab, y = unit(lineHeight, "lines"), gp = gp,...)
    grid.text(ylab, x = unit(lineHeight, "lines"), gp = gp,rot = 90,...)
    popViewport()
}



#############
# Turn scientific notation into x10^y

get10power <- function(X,minPower=0,extra=NULL,...){
    Y = format(X,scientific=TRUE,...)
    parts=str_split(Y,"e")
    dec = sapply(parts,function(i) i[1])
    ten = sapply(parts,function(i) as.numeric(i[2]))
    if(!is.null(extra)) Ynew = sapply(1:length(Y),function(i) as.expression(bquote(.(extra)~.(dec[i])~x10^.(ten[i])))) else Ynew = sapply(1:length(Y),function(i) as.expression(bquote(.(dec[i])~x10^.(ten[i]))))
    Ynew[is.na(X)]=NA
    if(minPower>0) Ynew[which(X<(10^minPower))] = as.character(X[which(X<(10^minPower))])
    if(minPower<0) Ynew[which(X>(10^minPower))] = as.character(X[which(X>(10^minPower))])

    return(Ynew)
}


get10power2 <- function(X,minPower=0){
    ten=log10(X)
    Ynew = sapply(1:length(X),function(i) as.expression(bquote(10^.(ten[i]))))
    Ynew[is.na(X)]=NA
    if(minPower>0) Ynew[which(X<(10^minPower))] = as.character(X[which(X<(10^minPower))])
    if(minPower<0) Ynew[which(X>(10^minPower))] = as.character(X[which(X>(10^minPower))])

    return(Ynew)
}



#get10power(c(0,NA,10,1000,100200))

####################################################
# EVERYTHING BELOW IS FOR POSTERIOR STUFF
####################################################

# Compute (approximate) Bayes Factors from betas and ses.

                                        # Formlar, for a given SNP:
# B = effect size estimate
# SE = standard error of the effect size estimate
# g = sqrt( prior on variance of effect sizes under the null )
# BF = N(0,SE^2)/N(0,SE^2+g^2) = sqrt(SE^2/(SE^2+g^2))*exp( g^2 * B^2 / ( 2* SE^2 ( SE^2 + g^2) ))


#bf2 = sqrt( (SE^2+g^2)/SE^2 ) * exp( -(g^2 * (B/SE)^2 )/ ( 2* ( SE^2 + g^2 ) ) )

calculate.BF <- function(B,SE,g=0.2){
    
    bf = sqrt( SE^2/(SE^2+g^2) ) * exp( g^2 * B^2 / ( 2* SE^2 * ( SE^2 + g^2 ) ) )

    # do the calculation in log space.
    logbf = log( sqrt( SE^2/(SE^2+g^2) ) ) +  g^2 * B^2 / ( 2* SE^2 * ( SE^2 + g^2 ) )
    
    #check = dnorm(B,0,sd=sqrt(SE^2+g^2))/dnorm(B,0,sd=SE)

    #logcheck = log( dnorm(B,0,sd=sqrt(SE^2+g^2)) ) - log( dnorm(B,0,sd=SE) )
        
    #return(cbind(bf,check,logbf,logcheck))
    return(cbind(bf,logbf))
}

compute.BF.posteriors <- function(BF){
    
    return(BF/sum(BF,na.rm=TRUE))
}

compute.BF.posteriors2 <- function(BF,logBF,region.indexes=NULL){


    # Use the logs of BF to avoid extremely high values
    hasInfinite = (Inf %in% BF) | (-Inf %in% BF)
    
    if( hasInfinite ){
        
        myShrink = max(abs(logBF),na.rm=TRUE) - 600

        BF = exp(logBF - myShrink) # subtract something from the log values and take the exponential
        print(paste0(myShrink," shrunk logs."))
    }
    
    if(!is.null( region.indexes)){
        
        totals = aggregate(BF,sum,na.rm=TRUE,by=list(region.indexes))
                                        #sapply(dataThis$region.index)
        if( ("Inf"%in%totals[,2])| (NA%in%totals[,2]) ) print("ERROR: You have infinite or NA values in your totals. Consider increasing the shrinkage factor.")
        out = BF/totals[match(region.indexes,totals$Group.1),2]
        
    } else {
        out = BF/sum(BF,na.rm=TRUE)
    }
    return(out)
}





getQCsubset <- function(MergedData,minmaf,mininfo,maxmiss){
    
    MergedDataQC = lapply(MergedData,function(x){

        print(paste0('chromosome: ',x$CHR[1]))
        
        toExclude = sapply(types,function(t){
                                        #print(t)
            mafcol = paste0(t,"..MAF")
            infocol = paste0(t,"..INFO")
            misscol = paste0(t,"..F_MISS")
                                        # get a list of exclusions
            
            excl1 = excl2 = excl3 = c()
            if(mafcol %in% colnames(x)) excl1 = which(x[,mafcol] <= minmaf)
            if(infocol %in% colnames(x)) excl2 = which(x[,infocol] <= mininfo)
            if(misscol %in% colnames(x)) excl3 = which(x[,misscol] >= maxmiss)

            excl = 1:nrow(x) %in% c(excl1,excl2,excl3)
            
           # print(paste0("excluding ",sum(excl)," SNPs from ",t," because of QC."))

            if(sum(excl)>0) {
                print( paste0("MAF: ",length(excl1)) )
                print( paste0("INFO: ",length(excl2)) )
                print( paste0("MISS: ",length(excl3)) )       
#                print(head(x[excl,grep(t,colnames(x))]))
            }
                                        # spit out an T/F list of exclusions
            return(excl)
        })
        
        return(toExclude)
    })

    totalsExcluded = sapply(types,function(t) {
        sapply(MergedDataQC,function(x) sum( x[,t] ) )
    })
        
    return(list( MergedDataQC , totalsExcluded )) 
}


extractInfoForRegions <- function(chrs,region.indices,prior=NULL,columns=NULL,MergedData,MergedDataQC,scaleFactor=1,allOverlapping=FALSE,types=c("IMP","GENO","GIANT"),mustBeSignificant=FALSE) {

    print(head(chrs))
    print(head(region.indices))
    
    theChroms = lapply(unique(chrs),function(chr){
        
        print(chr)
        
        chromIndex = which(MergedDataChromOrder==chr)
        data = MergedData[[chromIndex]]
        qc = MergedDataQC[[chromIndex]]

        region.indexes = unique(region.indices[chrs==chr])
        
        if(allOverlapping){
            these = (data$region.index %in% region.indexes) & (rowSums(qc[,types])==0) & ( !is.na(rowSums(data[,paste0(types,c("..Posterior"))]) ) )
        } else {            
            these = data$region.index%in%region.indexes
        }
        print(sum(these))
        data = data[these,]
        qc = qc[these,]
        
        
        dat = sapply(types,function(type){
            
            rawPost = data[,paste0(type,c("..Posterior"))]
            rawPost[qc[,type]] = NA
            print(type)
            
            if(!is.null(prior)){
            #   print(prior)
                seCol = grep(paste0(type,c("..SE")),colnames(data),value=TRUE)
                bCol = grep(paste0(type,c("..BETA")),colnames(data),value=TRUE)
                
                print(paste0("Re-computing bayes factors using prior ",prior))

                if(type!="GIANT") {
                    print(paste0("Scaling BETAs and SEs for ",type," by ",scaleFactor))
                } else {
                    scaleFactor=1
                }
                betas = data[,bCol]/scaleFactor
                ses = data[,seCol]/scaleFactor
                
                bfNew = calculate.BF(betas,ses,g=prior)
                
                bf = bfNew[,"bf"]
                logbf = bfNew[,"logbf"]
                
            } else {
                bfcol = grep(paste0(type,c("..bf")),colnames(data),value=TRUE)
                logbfcol = grep(paste0(type,c("..logbf")),colnames(data),value=TRUE)

                print(bfcol)
                print(logbfcol)
                
                bf = data[,bfcol]
                logbf = data[,logbfcol]
            }

            # APPLY THE QC HERE, BEFORE COMPUTING POSTERIORS
            bf[qc[,type]] = NA
            logbf[qc[,type]] = NA

            if( mustBeSignificant ){
                
                notSignif = which(data[,paste0(type,c("..P"))] >= signifLevel)
                print(paste0(length(notSignif)," snps not significant."))
                print(head(notSignif))
                print(head(bf))
                print(length(logbf))
                print(length(bf))

                bf[notSignif] = NA
                logbf[notSignif] = NA
            }

            print(length(logbf))
            print(length(bf))

            out = compute.BF.posteriors2(bf,logbf,data$region.index)

            return(out)
        })

        colnames(dat) = paste0(types,"..MyNewPosteriors")
        
        dat = cbind(data[,c("region.index","CHR","BP")],dat)
        
        
        if(!is.null(columns)) {
            print(columns)
            foundColumns = unique(unlist(sapply(columns,function(x) grep(paste0(x,"$"),colnames(data),value=TRUE))))
            dat = cbind(dat,data[,colnames(data)%in%foundColumns])
        }
        
        return(dat)
    })

    print("hello")
    OUT = abind(theChroms,along=1,force.array=FALSE)
    
    overlappingSNPs = is.na(OUT[,paste0(types,"..MyNewPosteriors")])
    
    
    if(allOverlapping) {
        keep = rowSums(overlappingSNPs)==0 # nothing is missing
        print(paste0(sum(keep)," SNPs that are in all studies."))
       
    } else {
        keep = rowSums(!overlappingSNPs)>=1 # at least one is not missing
        print(paste0(sum(keep)," SNPs that is still in at least one study."))
    }

    
    OUT = OUT[keep,]

    OUT$region_ID = paste0(OUT$CHR,"_",OUT$region.index)

    return(OUT)
}





getTopHits <- function(type="GIANT",MergedData,MergedDataQC){

    topHitsGIANT = lapply(1:nChroms,function(i){
        
        raw = MergedData[[i]]
        qc = MergedDataQC[[i]]
        x = raw[!qc[,type],]

        region.indexes = x[,"region.index"]
        ranks = order(x[,paste0(type,"..P")]) # Ordered by p-value (might not be the same as posterior?)
        o = region.indexes[ranks]
        first = which(!duplicated(o))
        these = ranks[first]
        print(length(these))
        #print(head(x))
        #print(grep("..P",colnames(x)))
        topHits = x[these,c("CHR","region.index","SNP2","GIANT..SNP2.swap",grep("\\.\\.P$",colnames(x),value=TRUE))]
        
        return(topHits)
    })

    topHitsGIANT = abind(topHitsGIANT,along=1,force.array=FALSE)
    topHitsGIANT = topHitsGIANT[order(topHitsGIANT[,paste0(type,"..P")]),]
    topHitsGIANT$region_ID = paste0(topHitsGIANT$CHR,"_",topHitsGIANT$region.index)
    
    if( sum(duplicated(topHitsGIANT$region_ID))>0 ) print("WARNING: some top snps counted more than once.")
    
    return(topHitsGIANT)
}

getRegionInfo <- function(chr,reg,AllRegions){
    chromOrder = sapply(AllRegions,function(x) x$chr[1])
    d = AllRegions[[which(chromOrder==chr)]]
    o = d[reg,]
    return(o)
}


getNonOverlappingHits <- function(topHitsGIANT,AllRegions){
                                        # NOTE: Some regions might overlap. This is because of the region definition.
    overlappingRegions = lapply(AllRegions,function(regions){
        
        check = lapply(1:nrow(regions),function(i) {
            x = regions[i,]
            start = x["start"][1,]
            end = x["end"][1,]
            o = which(! ( (regions[,"end"] < start) | (regions[,"start"] > end) ))
            return(o)
        })
    })

                                        # Go down the list, and skip anything that is overlapping with an existing one.

    # first make sure that topHitsGIANT are ordered by the original order!!
    topHitsGIANT = topHitsGIANT[order(topHitsGIANT$CHR,topHitsGIANT$region.index),]
    
    k = 1
    toKeep = c()
    toKeepIDs = c()
    toSkip = c()
    nHits=0
    maxHits=nrow(topHitsGIANT)

    while( (nHits < maxHits) & (k <= nrow(topHitsGIANT)) ){
        #print(k)
        me = topHitsGIANT[k,]
        ch = overlappingRegions[[which(MergedDataChromOrder==me$CHR)]]
        id = me$region_ID
        hasPair = ch[[me$region.index]]
        pairID = paste0(me$CHR,"_",hasPair)
        
        if(length(intersect(pairID,toKeepIDs))==0){
            toKeepIDs = c(toKeepIDs,id)
            toKeep = c(toKeep,k)
            nHits = nHits + 1
        } else {
            toSkip = c(toSkip,k)
            print(paste0("Skipping the ",k,"th highest hit."))
        }    
        k = k + 1
    }

                                        #print(paste0(length(overlaps)," regions actually overlap one-another!"))

    nRegions = length(toKeep)
    print(paste0(nRegions," out of ",nrow(topHitsGIANT)," GIANT regions kept after pruning for overlaps."))

    regionsToCheck = topHitsGIANT[toKeep,]
    
    return(regionsToCheck)
}


plotRegionPosteriors <- function(X,chr,reg,theseTypes=types,cols=c(impColour,genoColour,giantColour),title="",shapes=c(21,23,24),alpha=0.5,tol=0.00001,maxY=150,buffer=25000,...){
    
    Xsub = X[X$CHR==chr & X$region.index == reg,]
    N = nrow(Xsub)
        
    BPrange = getRegionInfo(chr,reg,AllRegions)
    print(BPrange)
    xLims = c(BPrange[1,"start"],BPrange[1,"end"])
    
    if(!is.null(tol)){
        PostRange = rowSums(Xsub[,paste0(theseTypes,"..MyNewPosteriors")],na.rm=TRUE)
        xLims = range(Xsub$BP[PostRange>tol])+c(-1000,1000)
    }
    print(xLims)
    
    scale = "MB"
    scaleF = 1e6
    if(diff(xLims) < 10000) {
        scale = "KB"
        scaleF = 1e3
    }

    xLims = xLims/scaleF
    print(xLims)

    plot(NULL,xlim=xLims,ylim=c(0,1.05),col=NA,xlab=paste0("Position on chromosome ",chr," (",scale,")"),ylab="Posterior probability",main=title)

    Ns = c()
    for(type in theseTypes){

        if(is.null(tol)) abline(v=xLims+c(buffer,-buffer)/scaleF,lty=3,col="gray") # 25KB from the ends.
        
        y = Xsub[,paste0(type,"..MyNewPosteriors")]
        x = Xsub[,"BP"]/scaleF
        minP = min(Xsub[,paste0(type,"..P")],na.rm=TRUE)
        topHit = which(Xsub[,paste0(type,"..P")]==minP)
        Ns = c(Ns,sum(!is.na(y)))
        print(head(x))
        print(head(y))
        colIndex = which(theseTypes==type)
        points(x,y,xlim=BPrange,col=cols[colIndex],pch=shapes[colIndex],bg=add.alpha(cols[colIndex],alpha),...)
        points(x[topHit],y[topHit],col="black",pch=shapes[colIndex],bg="transparent",...)

    }
    legend("topleft",legend=paste0(theseTypes,": ",Ns),col=cols,pch=shapes,pt.bg=add.alpha(cols,alpha),bty="n",horiz=TRUE,...)

                                        # Plot the p-values

    PostRange = rowSums(Xsub[,paste0(theseTypes,"..MyNewPosteriors")],na.rm=TRUE)
    P2 = -log10(Xsub[,paste0(theseTypes,"..P")])
    Pmax = max( P2[which( PostRange>0 & P2!=Inf),] ,na.rm=TRUE)
    maxY = Pmax*1.5
    
    print(Pmax)
    print(maxY)
    
    plot(NULL,xlim=xLims,ylim=c(0,maxY*1.05),col=NA,xlab=paste0("Position on chromosome ",chr," (",scale,")"),ylab="-log10(p)",main=title)
    abline(h=maxY,col="gray",lty=3)

    Ns = c()
    for(type in theseTypes){
        
        y = P2[,paste0(type,"..P")]
        y[which(y>Pmax)] = maxY
        
        x = Xsub[,"BP"]/scaleF
        minP = min(Xsub[,paste0(type,"..P")],na.rm=TRUE)
        topHit = which(Xsub[,paste0(type,"..P")]==minP)
        Ns = c(Ns,sum(!is.na(y)))
        print(head(x))
        print(head(y))
        colIndex = which(theseTypes==type)
        points(x,y,xlim=BPrange,col=cols[colIndex],pch=shapes[colIndex],bg=add.alpha(cols[colIndex],alpha),...)
        points(x[topHit],y[topHit],col="black",pch=shapes[colIndex],bg="transparent",...)

    }
    legend("topleft",legend=paste0(theseTypes,": ",Ns),col=cols,pch=shapes,pt.bg=add.alpha(cols,alpha),bty="n",horiz=TRUE,...)


    
}



getRegressionInfo <- function(r){
    s = summary(r)
    
    o = bquote( "slope = "* .(round(coef(r)[2],3)) *" ("* .( round(s$coefficients[2,"Std. Error"],4) )*")   "* r^2 * " = " * .(round(s$r.squared,3)) )
    
    return(o)
}


PlotPosteriorsAndEffectSizes <- function(X,plotName,theseRegions,nRegions=1:20,plotdir,plotHistograms=FALSE,types=c("IMP","GENO","GIANT"),effectSizesOnly=FALSE,...){

    # subset X
    ind = paste0(X$CHR,"_",X$region.index)
    keep = ind %in% theseRegions$region_ID
    X = X[which(keep),]
    print(nrow(X))

    if(!effectSizesOnly){

        pdf(paste0(plotdir,"/plots/region-posteriors-",plotName,".pdf"),bg="transparent",width=10,height=6)

        for(i in  nRegions ){
            snp = getRegionInfo(theseRegions$CHR[i],theseRegions$region.index[i],AllRegions)
            snpInfo = paste0(plotName,"; top hit in GIANT: ",snp[,"HITS.topVariants2"],",  ",snp[,"HITS.topVariants"],", -log10(p) = ",round(-log10(theseRegions$GIANT..P[i]),3))
            
            plotRegionPosteriors(X,chr=theseRegions$CHR[i],reg=theseRegions$region.index[i],title=snpInfo,theseTypes=types,...)
            
        }

        dev.off()
    }


    pdf(paste0(plotdir,"/plots/effect-sizes-",plotName,".pdf"),bg="transparent",width=10,height=10)

    par(...)
    
    giantBFlipped = X$GIANT..BETA
    giantBFlipped[X$GIANT.swapped] = - X$GIANT..BETA[X$GIANT.swapped]

                                        # exclude anything that doesn't overlap
                                        #    snp = getRegionInfo(theseRegions$CHR[i],theseRegions$region.index[i],AllRegions)

    for(type in types[types!="GIANT"]){

        print(type)
        
        x = giantBFlipped
        y = X[,paste0(type,"..BETA")]

        xse = X$GIANT..SE
        yse = X[,paste0(type,"..SE")]
            
        keep = !is.na(x) & !is.na(y) & !duplicated(X$SNP2) & X$GIANT..P < signifLevel
        print(sum(keep))
        
        x = x[keep]
        y = y[keep]
        xse = xse[keep]
        yse = yse[keep]

        yScaled = y/scaleFactor
        
        regressionBraw = lm(y~x)
        sBraw = summary(regressionBraw)

        regressionB = lm(yScaled~x)
        sB = summary(regressionB)

        xz = x/xse
        yz = y/yse

        regressionZ = lm(yz~xz)
        sZ = summary(regressionZ)

        Waldx = mean(xz,na.rm=TRUE)-xz
        Waldy = mean(yz,na.rm=TRUE)-yz

        regressionW = lm(Waldy~Waldx)
        sW = summary(regressionW)


        topSnps = which(X$SNP2[keep]%in%theseRegions$SNP2)

        pointCol= add.alpha("black",0.5)
        
        plot(x,y,main=paste0(plotName,"; ",sum(keep)," SNPs."),xlab="Effect size in GIANT (2014)",ylab=paste0("Raw effect size in ",sourceNames[type]),pch=16,col=pointCol)
        points(x[topSnps],y[topSnps],col="red")
        abline(0,1,col="red",lty=3)
        abline(coef(regressionBraw),col="blue")
        legend("top",legend=getRegressionInfo(regressionBraw),col=NA,bty="n")


        plot(x,yScaled,main=paste0(plotName,"; ",sum(keep)," SNPs."),xlab="Effect size in GIANT (2014)",ylab=paste0("Scaled effect size in ",sourceNames[type] ),pch=16,col=pointCol)
        points(x[topSnps],yScaled[topSnps],col="red")

        abline(0,1,col="red",lty=3)
        abline(coef(regressionB),col="blue")
        legend("top",legend=getRegressionInfo(regressionB),col=NA,bty="n")

        l = max(abs(c(xz,yz))); limits = c(-l,l)
        
        plot(xz,yz,xlim=limits,ylim=limits,main=paste0(plotName,"; ",sum(keep)," SNPs."),xlab="Z-score in GIANT (2014)",ylab=paste0("Z-score in ",sourceNames[type] ),pch=16,col=pointCol)
        points(xz[topSnps],yz[topSnps],col=giantColour)

        abline(0,1,col="black",lty=3,lwd=1.5)
        abline(coef(regressionZ),col="red",lwd=1.5)
        legend("top",legend=getRegressionInfo(regressionZ),col=NA,bty="n")

                                        #plot(Waldx,Waldy,main=paste0(plotName,"; ",sum(keep)," SNPs."),xlab="Wald stat in GIANT (2014)",ylab="Wald stat in UK Biobank (imputed data)")
                                        #abline(0,1,col="red",lty=3)
                                        #abline(coef(regressionW),col="blue")
                                        #legend("top",legend=getRegressionInfo(regressionW),col=NA,bty="n")

        if(plotHistograms){
            hist(xz,xlab="Z-score in GIANT (2014)",freq=FALSE)
            curve(dnorm(x, mean=0, sd=1),col="darkblue",add=TRUE, lwd=2,yaxt="n")

            hist(yz,xlab="Z-score in UK Biobank (imputed data)",freq=FALSE)
            curve(dnorm(x, mean=0, sd=1),col="darkblue",add=TRUE, yaxt="n")

            hist(Waldx,xlab="Wald stat in GIANT (2014)",freq=FALSE)
            curve(dnorm(x, mean=0, sd=1),col="darkblue",add=TRUE, yaxt="n",freq=FALSE)

            hist(Waldy,xlab="Wald stat in UK Biobank (imputed data)",freq=FALSE)
            curve(dnorm(x, mean=0, sd=1),col="darkblue",add=TRUE, yaxt="n")
        }

    }
    
    dev.off()

       
}

## credible sets analysis

summarisePosteriorComparison <- function(X,theseRegions,theseTypes,bins=seq(0.8,1,by=0.01)){
                                        # subset X
    ind = paste0(X$CHR,"_",X$region.index)
    
    keep = ( ind %in% theseRegions$region_ID ) & ( !( (is.na(X[,paste0(theseTypes[1],"..MyNewPosteriors")])) & (is.na(X[,paste0(theseTypes[2],"..MyNewPosteriors")])) ) ) 
    
    Xsub = X[which(keep),]
    print(nrow(Xsub))
    indSub = ind[keep]
        
    PosType1 = Xsub[,paste0(theseTypes[1],"..MyNewPosteriors")]
    PosType2 = Xsub[,paste0(theseTypes[2],"..MyNewPosteriors")]
    PValType1 = Xsub[,paste0(theseTypes[1],"..P")]
    PValType2 = Xsub[,paste0(theseTypes[2],"..P")]
    
    uniqueRegions = unique(ind)
    nRegions = length(uniqueRegions)
    
    if(nRegions!=nrow(theseRegions)) print(paste0("You asked for ",nrow(theseRegions)," but some are missing in your X data. Only plotting ",nRegions," regions."))

    regSummaries = lapply(1:nrow(theseRegions),function(k){
        
        i = theseRegions$region_ID[k]
        if(!i %in% uniqueRegions) return(NULL)
        
        if(k%%100==0 ) print(k)
        
        these = indSub==i
            
        pA = PosType1[these]
        pB = PosType2[these]

        minPvalA = min(PValType1[these],na.rm=TRUE)
        minPvalB = min(PValType2[these],na.rm=TRUE)

        if(minPvalA >= signifLevel) print(paste0(theseTypes[1]," has no significant snp! For region ",k," ",i))

        if(minPvalB >= signifLevel) print(paste0(theseTypes[2]," has no significant snp! For region ",k," ",i))
            
        orderA = order(pA,decreasing=TRUE)
        orderB = order(pB,decreasing=TRUE)
        a = cumsum(pA[orderA])
        b = cumsum(pB[orderB])

        if( length(a) != length(b) ) print("ERROR !!")
        
        summaryOfSets = sapply(bins, function(x) {
            
            inSetA = unique(orderA[ 1: (sum(a < x,na.rm=TRUE) + 1) ])
            inSetB = unique(orderB[ 1: (sum(b < x,na.rm=TRUE) + 1) ])

            if(x==1) {
                inSetA = which(!is.na(pA))
                inSetB = which(!is.na(pB))
            }

            sizeA = length(inSetA)
            sizeB = length(inSetB)
            overlapAB = sum( pB[inSetA] ,na.rm=TRUE) # how much does x% of the posterior in A cover in B?
            overlapBA = sum( pA[inSetB] ,na.rm=TRUE) # how much does x% of the posterior in B cover in A?
            overlapBnotA = sum( pB[ inSetB[!inSetB%in%inSetA] ] ,na.rm=TRUE) # how much of B's posterior is covered by snps not in the credible set of A?

            overlapAnotB = sum( pA[ inSetA[!inSetA%in%inSetB] ] ,na.rm=TRUE) # how much of A's posterior is covered by snps not in the credible set of B?
            
            overlapSize = length(intersect(inSetA,inSetB))
            unionSize = length(unique(c(inSetA,inSetB)))
            
            o = cbind(sizeA,sizeB,overlapSize,unionSize,overlapAB,overlapBA,overlapBnotA,overlapAnotB)
            
            return(o)
        },simplify=FALSE)

        summaryOfSets = abind(summaryOfSets,along=1,force.array=FALSE)
        rownames(summaryOfSets) = bins
        
        return(summaryOfSets)
    })

    return(regSummaries)
}



getIn95CredibleSet <- function(X,theseRegions,theseTypes,BIN=0.95){
    
    ind = paste0(X$CHR,"_",X$region.index)
    
    keep = ( ind %in% theseRegions$region_ID ) & ( !( (is.na(X[,paste0(theseTypes[1],"..MyNewPosteriors")])) & (is.na(X[,paste0(theseTypes[2],"..MyNewPosteriors")])) ) ) 
    
    Xsub = X[which(keep),]
    print(nrow(Xsub))
    indSub = ind[keep]
    
    PosType1 = Xsub[,paste0(theseTypes[1],"..MyNewPosteriors")]
    PosType2 = Xsub[,paste0(theseTypes[2],"..MyNewPosteriors")]
    PValType1 = Xsub[,paste0(theseTypes[1],"..P")]
    PValType2 = Xsub[,paste0(theseTypes[2],"..P")]
    
    uniqueRegions = unique(ind)
    nRegions = length(uniqueRegions)
    
    if(nRegions!=nrow(theseRegions)) print(paste0("You asked for ",nrow(theseRegions)," but some are missing in your X data. Only plotting ",nRegions," regions."))

    regSummariesBIN = lapply(1:nrow(theseRegions),function(k){
        
        i = theseRegions$region_ID[k]
        if(!i %in% uniqueRegions) return(NULL)
        
        if(k%%100==0 ) print(k)
        
        these = indSub==i
        
        pA = PosType1[these]
        pB = PosType2[these]

        minPvalA = min(PValType1[these],na.rm=TRUE)
        minPvalB = min(PValType2[these],na.rm=TRUE)

        if(minPvalA >= signifLevel) print(paste0(theseTypes[1]," has no significant snp! For region ",k," ",i))

        if(minPvalB >= signifLevel) print(paste0(theseTypes[2]," has no significant snp! For region ",k," ",i))
        
        orderA = order(pA,decreasing=TRUE)
        orderB = order(pB,decreasing=TRUE)
        a = cumsum(pA[orderA])
        b = cumsum(pB[orderB])

        if( length(a) != length(b) ) print("ERROR !!")
        
        inSetA = unique(orderA[ 1: (sum(a < BIN,na.rm=TRUE) + 1) ])
        inSetB = unique(orderB[ 1: (sum(b < BIN,na.rm=TRUE) + 1) ])

        inABoolian = 1:length(a)%in%inSetA
        inBBoolian = 1:length(b)%in%inSetB

        D = cbind(inABoolian,inBBoolian)
        colnames(D) = c(paste0(theseTypes[1],"..in95credibleSet"),paste0(theseTypes[2],"..in95credibleSet"))

        o = cbind(Xsub[these,],D)
        return(o)            
    })
    
    return(regSummariesBIN)
    
}



extractCredibleData <- function(CredibleSets,theseRegions){
    
        #o = cbind(sizeA,sizeB,overlapSize,overlapAB,overlapBA)
    NonNull = which(sapply(CredibleSets,function(x) !is.null(x)))
    nRegs = length(NonNull)
    
    regs = theseRegions[NonNull,]
    CredibleSets = CredibleSets[NonNull]

    bins = rownames(CredibleSets[[1]])

    # extract sizes
    sizes = lapply(bins,function(y) {
        p = sapply(CredibleSets,function(x) x[y,c("sizeA","sizeB")])
        return(p)
    })
    
    totals = sapply(CredibleSets,function(x) x["1",c("sizeA","sizeB")])
        
    # extract overlap sizes
    overlapsizes = lapply(bins,function(y) {
        p = sapply(CredibleSets,function(x) x[y,c("overlapSize")])
        return(p)
    })
    unionsizes = lapply(bins,function(y) {
        p = sapply(CredibleSets,function(x) x[y,c("unionSize")])
        return(p)
    })

    # First one is always summing the posteriors for type B
    overlapABBA = lapply(bins,function(y) {
        p = sapply(CredibleSets,function(x) x[y,c("overlapAB","overlapBA")])
        return(p)
    })

    overlapABnotBA = lapply(bins,function(y) {
        p = sapply(CredibleSets,function(x) x[y,c("overlapBnotA","overlapAnotB")])
        return(p)
    })


    
    binsN = as.numeric(bins)

    bin95 = which(binsN==0.95)
    
    print("Size-1 numbers")
    print(sum(sizes[[bin95]][1,]==1))
    print(sum(sizes[[bin95]][2,]==1))

    print("median size")
    print(median(sizes[[bin95]][1,]))
    print(median(sizes[[bin95]][2,]))

    print("median size per snp")
    print(median(sizes[[bin95]][1,]/totals[1,]))
    print(median(sizes[[bin95]][2,]/totals[2,]))

    print("total number markers")
    print(sum(totals[1,]))
    print(sum(totals[2,]))

    
    return(list("sizes"=sizes,"totals"=totals,"overlapsizes"=overlapsizes,"unionsizes"=unionsizes,"overlapABBA"=overlapABBA,"overlapABnotBA"=overlapABnotBA,"binsN"=binsN,"nRegs"=nRegs))
    
}


PlotCredibleSets <- function(CredibleSets,theseRegions,type1,type2,plotbins=c(0.95),...){
                                        # Output from summarisePosteriorComparison()

    colours = c(impColour,impColour,genoColour,giantColour)[match(c(type1,type2),names(sourceNames))]

    print(type1)
    print(type2)
    DATA = extractCredibleData(CredibleSets)

    sizes = DATA$sizes
    totals = DATA$totals
    overlapsizes = DATA$overlapsizes
    unionsizes = DATA$unionsizes
    overlapABBA = DATA$overlapABBA
    overlapABnotBA = DATA$overlapABnotBA
    binsN = DATA$binsN
    nRegs = DATA$nRegs
    

    # sizes
    for(b in plotbins){
        
        print("sizes")
        
        x = sizes[[which(binsN==b)]]
        fracType2Smaller = sum(x[2,]<x[1,])/ncol(x)
        fracType1Smaller = sum(x[1,]<x[2,])/ncol(x)

        #print(x[2,]<x[1,])
        
        plot(x=x[1,],y=x[2,],ylim=c(0,max(x)),xlim=c(0,max(x)),xlab=sourceNames[type1],ylab=sourceNames[type2],col=add.alpha("blue",0.5),main=paste0("Credible set size of ",nRegs," regions at level ",b),...)
        legend("topright",bty="n",legend=c(paste0(round(fracType1Smaller*100,2)," % of regions are smaller in ",type1),paste0(round(fracType2Smaller*100,2)," % of regions are smaller in ",type2)),col=NA)
        abline(0,1,col="red")


#### Hex-bin plot

        x1 = x[1,]
        x2 = x[2,] 

        myTheme = theme_bw() + theme(aspect.ratio=1)
        maxBreaks = max(c(x1[(x1<20)&(x2<20)],x2[(x1<20)&(x2<20)]))
        my_breaks = 2^c(0:ceiling(log2(maxBreaks)))
        data = data.frame(cbind(x1,x2))

        d <- ggplot(data[(x1<20)&(x2<20),], aes(x=x1, y=x2)) + stat_binhex(bins=20) +  scale_fill_gradientn(name="Count" ,breaks=rev(my_breaks),trans="log2",colors=viridis_pal()(10)) +  scale_x_continuous( expand = c(0.1, 0.1),name=sourceNames[type1]) +
            scale_y_continuous(expand = c(0.1, 0.1),name=sourceNames[type2]) + geom_abline(slope=1, intercept=0,color="red",linetype="longdash") + myTheme + ggtitle(paste0(b*100,"% credible set sizes for ",nRegs," (zoomed in)"))

        print(d)
        
####
        plot(x=x[1,],y=x[2,],ylim=c(0,20),xlim=c(0,20),xlab=sourceNames[type1],ylab=sourceNames[type2],col=add.alpha("blue",0.1),main=paste0("Credible set size of ",nRegs," regions at level ",b," (zoomed in)"),cex=3,...)
        abline(0,1,col="red")

        oo = abind(sapply(1:length(my_breaks),function(w) rep(w,my_breaks[w])),along=1)
        plot(rep(1,length(oo)),oo,xlab=NA,axes=FALSE,ylab=NA,col=add.alpha("blue",0.1),main="Count",pch=16,cex=8,xpd=NA)
        text(x=rep(1,length(my_breaks)),y=1:length(my_breaks),labels=my_breaks,cex=2,xpd=NA,adj=0.5)
            
        # As fraction of SNPs in region.
        plot(x=x[1,]/totals[1,],y=x[2,]/totals[2,],xlab=sourceNames[type1],ylab=sourceNames[type2],col=add.alpha("blue",0.3),main=paste0("Proportion of markers in ",b*100,"% credible set for ",nRegs," regions"),...)
        abline(0,1,col="red")

        x1 = x[1,]/totals[1,]
        x2 = x[2,]/totals[2,]

        lim = 0.05
        data = data.frame(cbind(x1,x2))
        #maxBreaks = max(data)
        #my_breaks = 2^c(0:ceiling(log2(maxBreaks)))

        d <- ggplot(data, aes(x=x1, y=x2)) + stat_binhex(bins=100) +  scale_fill_gradientn(name="Count" ,trans="log2",colors=viridis_pal()(10)) +  scale_x_continuous( expand = c(0.01, 0.01),name=sourceNames[type1]) +
            scale_y_continuous(expand = c(0.01, 0.01),name=sourceNames[type2]) + geom_abline(slope=1,color="red",linetype="longdash") + myTheme + ggtitle(paste0("Proportion of markers in ",b*100,"% credible set for ",nRegs," regions (zoomed in)"))

        print(d)
        


         # Histgram of the differences
        hist(x[2,]-x[1,],xlab=paste0("Size difference (",sourceNames[type2], " minus ",sourceNames[type1],")"),breaks=50,main = paste0("Difference in credible set sizes for ",nRegs," regions at level ",b))

        toPlot = log(x[2,])-log(x[1,])
        seqs = seq(0,max(abs(toPlot))+0.5,1)+0.5
        hist(toPlot,xlab=paste0("Log of the size ratio (",sourceNames[type2], " / ",sourceNames[type1],")"),breaks=c(-rev(seqs),seqs),main = paste0("Ratio of credible set size for ",nRegs," regions at level ",b))


         # barplot of the sizes
        specialBinsBarplot(x[1,],main=paste0("Credible set sizes in ",sourceNames[type1]," for ",nRegs," regions"),ylab="Number of regions",xlab=paste0("Number of markers in ",b*100,"% credible set"),col=colours[1])
        
        specialBinsBarplot(x[2,],main=paste0("Credible set sizes in ",sourceNames[type2]," for ",nRegs," regions"),ylab="Number of regions",xlab=paste0("Number of markers in ",b*100,"% credible set"),col=colours[2])
        
        # together
        specialBinsBarplot2(x,main=paste0("Size of ",b*100,"% credible sets \nfor ",nRegs," regions genome-wide significant in GIANT"),ylab="Number of regions",xlab=paste0("Number of markers in ",b*100,"% credible set"),col=colours,legend=sourceNames[c(type1,type2)])
        
    }

        # overlaps
    for(b in plotbins){
        
        x = overlapABBA[[which(binsN==b)]]
        plot(x=x[1,],y=x[2,],
             xlab=paste0(sourceNames[type2]," posterior for markers in the ",sourceNames[type1]," credible set"),
             ylab=paste0(sourceNames[type1]," posterior for markers in the ",sourceNames[type2]," credible set"),
             col=add.alpha("blue",0.5),
             main=paste0("Credible sets of ",nRegs," regions at level ",b),...)
        abline(0,1,col="red")

        m = paste0("Posterior in ",sourceNames[type2]," \nfor markers in the ",sourceNames[type1]," ",b*100,"% credible set")
        hist(x[1,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type2]),breaks=10,main=m)

        m = paste0("Posterior in ",sourceNames[type1]," \nfor markers in the ",sourceNames[type2]," ",b*100,"% credible set")
        hist(x[2,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type1]),breaks=10,main=m)


        # Histgram of the values
        x = overlapABnotBA[[which(binsN==b)]]

        m = paste0("Markers in ",sourceNames[type2]," ",b*100,"% credible set \nbut not in the ",sourceNames[type1]," ",b*100,"% credible set")
        hist(x[1,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type2]) ,breaks=10,main=m)

        m = paste0("Markers in ",sourceNames[type1]," ",b*100,"% credible set \nbut not in the ",sourceNames[type2]," ",b*100,"% credible set")
        hist(x[2,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type1]) ,breaks=10,main=m)

       # Histgram of the overlap as fraction of the union
        x = overlapsizes[[which(binsN==b)]]/unionsizes[[which(binsN==b)]]
        m = paste0("Overlap of ",b*100,"% credible sets")
        hist(x,ylab="Number of regions",xlab=paste0("Fraction overlapping (markers in both studies)/(markers in either)") ,breaks=10,main=m)
        

    }
    
}



specialBinsBarplot <- function(s,myBins=c(1:10,11,21,51),...){
                                        #s = x[1,]
                                        # myBins = c(1:10,11,21,51)
        myBins = c(myBins,max(s)+1)
        nBins = length(myBins)-1
        tab = sapply(1:nBins,function(r) sum(( s >= myBins[r] ) & ( s < myBins[r+1] )))
        
        labels = sapply(1:nBins,function(r) {
            if(myBins[r+1]-myBins[r]==1) o=myBins[r] else o=paste0(myBins[r],"-",myBins[r+1]-1)
            if(r==nBins) o = paste0(myBins[r],"+")
            return(o)
            })
        if( ! sum(tab) == length(s) ) print("ERROR: with your bin sizes!")
        
        barplot(tab,names.arg=labels,...)


}




specialBinsBarplot2 <- function(s,myBins=c(1:10,11,21),...){
                                        #s = x[1,]
                                        # myBins = c(1:10,11,21,51)
        myBins = c(myBins,max(s)+1)
        nBins = length(myBins)-1
        tab1 = sapply(1:nBins,function(r) sum(( s[1,] >= myBins[r] ) & ( s[1,] < myBins[r+1] )))
        tab2 = sapply(1:nBins,function(r) sum(( s[2,] >= myBins[r] ) & ( s[2,] < myBins[r+1] )))

        tab = cbind(tab1,tab2)
        labels = sapply(1:nBins,function(r) {
            if(myBins[r+1]-myBins[r]==1) o=myBins[r] else o=paste0(myBins[r],"-",myBins[r+1]-1)
            if(r==nBins) o = paste0(myBins[r],"+")
            return(o)
            })
        if( ! sum(tab) == length(s) ) print("ERROR: with your bin sizes!")
        
        barplot(t(tab),names.arg=labels,beside=TRUE,...)

}


plotCredibleSizes <- function(DATA,type1,type2,b=0.95,Expand=1.5,ma=17){
    
#    par(cex.axis=1.5,cex.main=2,cex.lab=2,mar=c(6,6,2,2))
    sizes = DATA$sizes
    totals = DATA$totals
    overlapsizes = DATA$overlapsizes
    unionsizes = DATA$unionsizes
    overlapABBA = DATA$overlapABBA
    overlapABnotBA = DATA$overlapABnotBA
    binsN = DATA$binsN
    nRegs = DATA$nRegs
    
    ind = which(binsN==b)

    x = sizes[[ind]]
    x1=x[1,];x2=x[2,]
    x1persnp = x1/totals[1,]
    x2persnp = x2/totals[2,]

    #ma = quantile(c(x1,x2),0.75)
    #ma=18
    zoomed = (x1<=ma)&(x2<=ma)

    sizeFreqs1 = table(x1)
    sizeFreqs2 = table(x2)

    plot(x=x1,y=x2,xlab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]),ylab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type2]),pch=16,
         col=add.alpha("blue",0.3),cex=0.8)
    abline(0,1,col="red")

    
    plot(x=x1[zoomed],y=x2[zoomed],ylim=c(0,ma),xlim=c(0,ma),xlab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]),ylab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type2]),pch=16,
         col=add.alpha("blue",0.1),cex=0.8)
    abline(0,1,col="red")

    oo = abind(sapply(1:length(my_breaks),function(w) rep(w,my_breaks[w])),along=1)
    plot(rep(1,length(oo)),oo,xlab=NA,axes=FALSE,ylab=NA,col=add.alpha("blue",0.1),main="Count",pch=16,cex=8,xpd=NA)
    text(x=rep(1,length(my_breaks)),y=1:length(my_breaks),labels=my_breaks,cex=2,xpd=NA,adj=0.5)

### sized by number in group
    toplot = table(x1,x2)
    p = expand.grid(1:ma,1:ma)
    sizeScale = apply(p,1,function(x){
        toplot[rownames(toplot)==as.character(x[1]),colnames(toplot)==as.character(x[2])]
    })
    sizeScale = Expand*log2(sizeScale+1)

    maxBreaks = max(toplot)+1
    my_breaks = 2^c(0:ceiling(log2(maxBreaks)))

    
    plot(p[,1],p[,2],ylim=c(0,ma),xlim=c(0,ma),xlab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]),ylab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type2]),col=NA,
         main=paste0(sum(zoomed)," regions with credible sets with \nfewer than ",ma+1," markers in both studies"),axes=FALSE)
    box()

    axis(1,at=1:ma)
    axis(2,at=1:ma)
    
    abline(v=1:20,h=1:20,col="gray",lwd=0.5)
    abline(v=seq(5,20,5),h=seq(5,20,5),col="gray",lwd=1.5)

    abline(0,1,col="red")
    points(p[,1],p[,2],pch=16,cex=sizeScale,col="blue")

    text(y=rep(0,length(sizeFreqs1)),x = as.numeric(names(sizeFreqs1)),labels=sizeFreqs1,font=3)
    text(x=rep(0,length(sizeFreqs2)),y = as.numeric(names(sizeFreqs2)),labels=sizeFreqs2,font=3)

### Without numbers
    plot(p[,1],p[,2],ylim=c(0,ma),xlim=c(0,ma),xlab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]),ylab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type2]),col=NA,
         main=paste0(sum(zoomed)," regions with credible sets with \nfewer than ",ma+1," markers in both studies"),axes=FALSE,xaxs="i",yaxs="i")
    box()
    
    axis(1,at=1:ma)
    axis(2,at=1:ma)

    abline(v=1:(ma-1),h=1:(ma-1),col="gray",lwd=0.5)
    #abline(v=seq(5,ma-1,5),h=seq(5,20,5),col="gray",lwd=0.8)

    abline(0,1,col="red")
    points(p[,1],p[,2],pch=16,cex=sizeScale,col="blue")


    # Legend.
    oo = Expand*log2(my_breaks+1)
    plot(p[,1],p[,2],ylim=c(0,ma),xlim=c(0,ma),col=NA,axes=FALSE,xlab=NA,ylab=NA,main="Count")

    points(rep(1,length(oo)),seq(1,2*length(oo),by=2),col="blue",pch=16,cex=oo,xpd=NA)
    text(x=rep(1,length(oo))-1,y=seq(1,2*length(oo),by=2),labels=my_breaks,cex=1.5,xpd=NA)


    # on log scale
#    plot(x=x1,y=x2,xlab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]),ylab=paste0("Size of ",b*100,"% credible set in ",sourceNames[type2]),pch=16,
#         col=add.alpha("blue",0.1),main=,cex=3)
#    abline(0,1,col="red")

    log10x1=log10(x1persnp)
    log10x2=log10(x2persnp)

    ranges = range(c(log10x1,log10x2))
        
    plot(x=log10x1,y=log10x2,xlim=ranges,ylim=ranges,xlab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type1]),ylab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type2]),pch=16,xaxt="n",yaxt="n",
         col=add.alpha("blue",0.8),cex=1)
    axis(1,at=axTicks(1),labels=10^axTicks(1))
    axis(2,at=axTicks(2),labels=10^axTicks(2))

    abline(0,1,col="red")

    plot(x=x1persnp,y=x2persnp,xlim=c(0,0.1),ylim=c(0,0.1),xlab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type1]),ylab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type2]),pch=16,
         col=add.alpha("blue",0.3),cex=0.8)

    abline(0,1,col="red")

    plot(x=x1persnp[zoomed],y=x2persnp[zoomed],xlab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type1]),ylab=paste0("Fraction of markers in ",b*100,"% credible set for ",sourceNames[type2]),pch=16,
         col=add.alpha("blue",0.3),cex=0.8)

    abline(0,1,col="red")

    
    #### overlaps

    x = overlapABBA[[ind]]
    

    m = paste0("Posterior in ",sourceNames[type2]," \nfor markers in the ",sourceNames[type1]," ",b*100,"% credible set")
    hist(x[1,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type2]),breaks=10,main=m,col="darkgray")

    plot(x[1,],totals[2,],xlab=paste0("Credible set size in ",sourceNames[type2]),ylab=paste0("Posterior in ",sourceNames[type2]))

    m = paste0("Posterior in ",sourceNames[type1]," \nfor markers in the ",sourceNames[type2]," ",b*100,"% credible set")
    hist(x[2,],ylab="Number of regions",xlab=paste0("Posterior in ",sourceNames[type1]),breaks=10,main=m,col=giantColour)


    plot(x[2,],totals[1,],xlab=paste0("Credible set size in ",sourceNames[type1]),ylab=paste0("Posterior in ",sourceNames[type1]))

    print(toplot[1,1])
    
}
