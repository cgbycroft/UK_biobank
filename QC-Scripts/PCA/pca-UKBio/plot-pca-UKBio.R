## Script to plot output of submit-pca-UKBio.sh. This will be projections based on data from each batch.
library(zoo)
args=commandArgs(trailingOnly=T)

#args = c('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R','/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/readPSperformance.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/auxFunctions.R','-in','/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init2.snpload.map','-out','b1__b11-b001__b095-pca-UKbio-proj5-init2')


print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}


#################  SnpLoads should be the same file used in the submit-pca-UKBio.sh script.
BatchesFile = paste0(baseSampleQCDir,"/QC-Scripts/batchList.txt")
OutputFile = args[which(args%in%c("-out"))+1]
OutDir =  paste0(baseSampleQCDir,"/QC-Scripts/PCA/pca-UKBio")
SnpLoads = args[which(args%in%c("-in"))+1]
 
#################

system(paste0('mkdir ',OutDir,'/plots'))

BatchList = read.table(BatchesFile,header = FALSE, stringsAsFactors = FALSE)
nBatches = nrow(BatchList)

Info=read.multiple.batch.info(batchInfoFields[c(45,43,36,55:59)])
Info = dplyr::tbl_df(Info)

Rfile=paste0(baseSampleQCDir,"/data/PCA/",OutputFile,".RData")

if( Rfile %in% list.files(path=paste0(baseSampleQCDir,"/data/PCA"),full.names=T) ){
    
    load(Rfile,verbose=T)
    
} else {
    
    PCs = NULL

    for (b in seq(nBatches)) {
        
        ## PCsFile should be a file _without_ a header line,
        ## with nPCs + 1 columns. The first column should list the sample IDs,
        ## the other columns should correspond to PC1 to PCk.

        ## Csv file is the augmented sample table file which also contains
        ## the CTSU information, most importantly, the self-reported ethnicity.

        Batch = BatchList[b,1]
        writeLines(Batch)

        PCsFile = paste0(baseSampleQCDir,"/data/PCA/",Batch,"-PCA.pcs")
                                        #CsvFile = paste0(baseDataDir,"/",get.sample.table.name(Batch))
        
                                        #BatchInfo = get.batch.info(Batch,CsvFile)
        BatchPCs = read.table(PCsFile,header = FALSE,stringsAsFactors = FALSE)
                                        #Info = rbind(Info,BatchInfo)
        PCs = rbind(PCs,BatchPCs)
    }

    ## It is convenient to label the columns
    nPCs = ncol(PCs) - 1
    names(PCs) = c("PIID",paste0("PC",seq(nPCs)))

    PCs = dplyr::tbl_df(PCs)
    Info2 = dplyr::select(Info, PIID, Pops, Chars, Colors)
    PCs = dplyr::left_join(PCs, Info2, by = c("PIID" = "PIID"))

    save(PCs,file=Rfile)
    # print out text file with just 20 PCs and self-reported ethnicity
    write.table(PCs[,c("PIID",paste0("PC",1:20),"Pops")],file=paste0(baseSampleQCDir,"/data/PCA/",OutputFile,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    
}

loads = dplyr::tbl_df(read.table(SnpLoads,header=FALSE,stringsAsFactors=FALSE))

nPCs = sum(grepl("PC",colnames(PCs)))

print(nPCs)
head(PCs)


## Plot the first XX PCs by ethnic background

Colors = PCs[["Colors"]]
Chars = PCs[["Chars"]]
Order = order.by.number.occurrences(PCs[["Pops"]])

######

## Plot using a grid

png(paste0(OutDir,'/plots/',OutputFile,'-grid-%02d.png'),width=1500,height=1500,res=150)
# if number of PCs is less than 10, and also uneven...
if(nPCs <= 10) {
    
    plot.grid.PCs(as.data.frame(PCs[Order,grepl("PC",colnames(PCs))]),n=5,shapes=Chars[Order],colours=Colors[Order])
 
} else {
    sequence = seq(1,nPCs,10)

    for(l in sequence){
        print(l)
        plot.grid.PCs(as.data.frame(PCs[Order,paste0("PC",c(l:(l+min(9,(nPCs-l)) )))]),shapes=Chars[Order],colours=Colors[Order])
    }
}
dev.off()

#remove odd-numbered pages!
sapply(seq(1,nPCs,by=2),function(i) system(paste0("rm ",OutDir,'/plots/',OutputFile,'-grid-0',i,'.png')))

#quit()


png(paste0(OutDir,'/plots/',OutputFile,'-%02d.png'),width=1500,height=1500,res=150)

for (pc in 1:nPCs){
    
    if(pc%%2==0) next

    print(paste0('Plotting PCs',pc,' and ',pc + 1))

    x = PCs[[paste0("PC",pc)]]
    y = PCs[[paste0("PC",pc+1)]]

    plot(x[Order],y[Order],xlab=paste0('PC ',pc),ylab=paste0('PC ',pc+1),
         col=Colors[Order],pch=Chars[Order],axes=FALSE)    
    axis(1,cex.lab=1.3)
    axis(2,cex.lab=1.3)

}

dev.off()

# plot legend

png(paste0(OutDir,'/plots/',OutputFile,'-Legend.png'),width=1500,height=1500,res=150)
#png(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/plots/b1__b11-b001__b095-pca-UKbio-proj40-Legend.png'),width=1500,height=1500,res=150)
plot.new()
eths = ethnicities[ethnicityOrder]
legend("topleft",legend=eths,col=ethnicity2col[eths],pch=ethnicity2char[eths],pt.lwd=2,horiz=F,bty="n")
dev.off()

# Table of counts
png(paste0(OutDir,'/plots/',OutputFile,'-LegendCounts.png'),width=1500,height=1500,res=150)
png(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/plots/b1__b11-b001__b095-pca-UKbio-proj40-LegendCounts.png'),width=1500,height=1500,res=150)
plot.new()
tab = sort(table(PCs$Pops),decreasing=T)
eths = names(tab)
perc = 100*tab/sum(tab)
perc = c(paste0("\t(",round(perc,2),"%)\t"),"\t\t")
legend("topleft",legend=paste0(c(tab,sum(tab)),perc,c(eths,"Total")),col=c(ethnicity2col[eths],"transparent"),pch=c(ethnicity2char[eths],1),pt.lwd=2,horiz=F,bty="n")
dev.off()



## Plot the first XX PCs by place of birth
Info$Country.of.birth=get.place.of.birth(Info)
Info$Eth2 = ethnicity2pop(Info$Ethnic.background)
# get some colours --> get most common ancestry
# plot each continent - i.e ethnic background - separately (and make everyone else gray)

ethTable = table(Info$Eth2)
InfoSorted = left_join(select(PCs,PIID),Info,by=c("PIID"="PIID"))

if(paste0("./pca-UKBio-countryColours.RData")%in%list.files(full.names=T)) {
    load(file=paste0("./pca-UKBio-countryColours.RData"),verbose=T)
} else {

    continents = list("Europe"=ethnicities[c(10,15,20,6)],
        "Africa"=ethnicities[c(2,4,9,11)],
        "Asia"=ethnicities[c(3,7,8,12,14,18)],
        "Mixed"=ethnicities[c(5,16,21,22,23)],
        "OtherUnknown"=ethnicities[c(1,13,17,19)])

    contColours=list()    
    for(cont in names(continents)){
        Eths = continents[[cont]]
                                        # main colour based on largest ethnicity
                                        #mainColour = ethnicity2col[names(sort(ethTable[Eths],decreasing=T)[1])]
        print(cont)
        countries = table(InfoSorted$Country.of.birth[InfoSorted$Eth2%in%Eths])
        if(cont=="Europe") countryToCol = sort(countries[countries > 100],decreasing=T)
        if(cont!="Europe") countryToCol = sort(countries[countries > 20],decreasing=T)
        cols = generateColours(names(countryToCol))
        cols[names(countryToCol)=="Prefer_not_to_answer"] = "gray"    
        contColours[[cont]] = cols
    }
                                        # save colours for countries
    save(contColours,continents,file=paste0("./pca-UKBio-countryColours.RData"))    
}

for(cont in names(continents)){

    cols=contColours[[cont]]
    Colors = cols[[1]][InfoSorted$Country.of.birth]
    Chars = cols[[2]][InfoSorted$Country.of.birth]
    Colors[!InfoSorted$Country.of.birth%in%names(cols[[1]])] = "black"
    Colors[!InfoSorted$Eth2%in%continents[[cont]]] = "gray"
    Colors[is.na(Chars)] = "gray"
    Chars[is.na(Chars)] = 1
    
    png(paste0(OutDir,'/plots/',OutputFile,'-',cont,'-byCountry-%02d.png'),width=1500,height=1500,res=150)

    for (pc in 1:nPCs){
        
        if(pc%%2==0) next

        Order = order.by.number.occurrences(Colors)
        print(paste0('Plotting PCs',pc,' and ',pc + 1))

        x = PCs[[paste0("PC",pc)]]
        y = PCs[[paste0("PC",pc+1)]]

        plot(x[Order],y[Order],xlab=paste0('PC ',pc),ylab=paste0('PC ',pc+1),
             col=Colors[Order],pch=Chars[Order],axes=FALSE)    
        axis(1,cex.lab=1.3)
        axis(2,cex.lab=1.3)

    }

    dev.off()

    # same thing, but on a grid
    png(paste0(OutDir,'/plots/',OutputFile,'-',cont,'-byCountry-grid-%02d.png'),width=1500,height=1500,res=150)
                                        # if number of PCs is less than 10, and also uneven...
    if(nPCs <= 10) {
        
        plot.grid.PCs(as.data.frame(PCs[Order,grepl("PC",colnames(PCs))]),n=5,shapes=Chars[Order],colours=Colors[Order])
        
    } else {
        sequence = seq(1,nPCs,10)

        for(l in sequence){
            print(l)
            plot.grid.PCs(as.data.frame(PCs[Order,paste0("PC",c(l:(l+min(9,(nPCs-l)) )))]),shapes=Chars[Order],colours=Colors[Order])
        }
    }
    dev.off()

                                        #remove odd-numbered pages!
    sapply(seq(1,nPCs,by=2),function(i) system(paste0("rm ",OutDir,'/plots/',OutputFile,'-',cont,'-byCountry-grid-0',i,'.png')))

                                        #quit()

    
    png(paste0(OutDir,'/plots/',OutputFile,'-',cont,'-byCountry-Legend.png'),width=1500,height=1500,res=150)
    plot.new()
    leg = c(names(cols[[1]]),"Other continents")
    if(length(leg)>35) {
        legend("topleft",legend=leg[1:30],col=c(cols[[1]],"gray")[1:30],pch=c(cols[[2]],1)[1:30],pt.lwd=2,horiz=F)
        legend("topright",legend=leg[31:length(leg)],col=c(cols[[1]],"gray")[31:length(leg)],pch=c(cols[[2]],1)[31:length(leg)],pt.lwd=2,horiz=F)   
    } else {
        legend("topleft",legend=leg,col=c(cols[[1]],"gray"),pch=c(cols[[2]],1),pt.lwd=2,horiz=F)
    }
    dev.off()
}


## boxplot by batch

b = dplyr::left_join(PCs, Info, by = c("PIID" = "PIID"))
b = b[["Batch"]]

png(paste0(OutDir,'/plots/',OutputFile,'-Batch-boxplot-%02d.png'),width=1500,height=1500,res=150)

for (pc in 1:20){
    
    print(paste0('Plotting boxplot for PC',pc))

    x = PCs[[paste0("PC",pc)]]
        
    boxplot(x~b,ylab=paste0('PC ',pc),axes=FALSE,col=c("black","blue"),outcol=c("black","blue"),border=NA)    
    axis(1,cex.lab=1.3)
    axis(2,cex.lab=1.3)
}

dev.off()


#quit()


## Plot the first 20 PCs SNP loads
print('Plotting snp-loads for pc...')

loads$index = loads$V4
loads$col = loads$V1
index=0

for(chr in 1:22){
    these = loads$V1==chr
    loads$index[these] = loads$index[these] + index
    index = max(loads$index[these])
    loads$col[these] = c("blue","darkblue")[chr%%2 + 1]
}

x = loads[["index"]]
cols = loads[["col"]]

png(paste0(OutDir,'/plots/',OutputFile,'-loads-%02d.png'),width=1500,height=1000,res=150)

for (pc in 1:nPCs){
    print(pc)
    y = loads[[paste0("V",pc+7)]]
    plot(x,y,xlab='position',ylab=paste0('PC ',pc),col=cols)
}

dev.off()


# Find regions of high SNP-loads

png(paste0(OutDir,'/plots/',OutputFile,'-loads-ksmooth-%02d.png'),width=1500,height=1000,res=150)
par(mfrow=c(2,1),mar=c(3,4,1,1))
for (pc in 1:nPCs){
#for (pc in c(1,17,21,29) ){
    print(pc)    
    y = loads[[paste0("V",pc+7)]]
    ysd = sd(y)
    qs = quantile(y,c(0.001,0.999))

    # where are the high loads bunched together?
#    t=which(abs(y)>sd(y))
    
 #   r = rle(y)
  #  w = cumsum(r$lengths)
   # table()
    #t = table(loads$V1[which(y > ysd)])/table(loads$V1)
    #print(t)
    #print(mean(t))
    plot(NULL,pch=NA,ylim=c(0,0.05),xlim=c(min(loads[["index"]]),max(loads[["index"]])),xlab=NA,ylab="Kernel smoothed")
    
    for(i in 1:22){
        ind = loads$V1==i
        x = loads$index[ind]
#        newx = x - min(x) + 1
#        newy = rep(NA,max(newx))
        Y = y[ind]
#        newy[newx] = Y
      
#        r2 = sapply(Y,function(j){
#            todo = c(x[j] - window,x[j] + window)
#            sd(Y[(x < todo[2])&(x > todo[1])])
#            })
        p = ksmooth(x,abs(Y),kernel="normal",bandwidth=1000000)
        lines(p[[1]],p[[2]],col=c("blue","darkblue")[i%%2 + 1])
    }
    
    inover = (y < qs[1])|(y > qs[2])
    cols=loads$col
    cols[inover] = "red"
   # inoverSpaces = diff(loads[["index"]][inover])    
    
    plot(loads[["index"]],y,xlab='position',ylab=paste0('PC ',pc),col=cols)
    
}

dev.off()


## finished

print(paste0('Finished: Plots saved to: ',OutDir,OutputFile,'*.png'))
