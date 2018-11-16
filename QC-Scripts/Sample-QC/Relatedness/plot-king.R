## Script to plot output from King, following filtering of the output from filter-king-output.R

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R",
#"-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.kin0","-out","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/Relatedness/plots/","-ibs0","0.0012")

print(args)
h = args[-c(which(args%in%c("-in","-out","-ibs0")),1+which(args%in%c("-in","-out","-ibs0")))]
for(helperScript in h){
    source(helperScript)
}

KingOutFile = args[which(args=="-in")+1]
OutputFile = gsub(".kin0","",basename(KingOutFile))
OutDir = args[which(args=="-out")+1]
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])

## FUNCTIONS
otherInfo = read.multiple.batch.info(c("Ethnic.background","Place.of.birth.in.UK...north.co.ordinate","Pops","Chars","Colors"))

getPopColors <- function(data){
    popn1 = as.character(otherInfo$Pops[match(data$ID1,otherInfo$PIID)])
    popn2 = as.character(otherInfo$Pops[match(data$ID2,otherInfo$PIID)])
    pops = as.character(popn1)

    orders = sapply(1:length(popn1),function(x) order(c(popn1[x],popn2[x])))
    swap = orders[1,]==2
    popnTemp=popn1
    popn1[swap]=popn2[swap]
    popn2[swap]=popnTemp[swap]

    pops[popn1!=popn2]=paste0(popn1[popn1!=popn2],"--",popn2[popn1!=popn2])

    cols = ethnicity2col[pops]
    cols[is.na(cols)]="black"
    shapes = ethnicity2char[pops]
    shapes[is.na(shapes)]=1
    return(list(cols,shapes,pops))
}


## KING kinship coefficients for each degree
## get kinship data
print( paste0("Reading kinship data... ",KingOutFile) )
kin = read.table(KingOutFile,header=T,stringsAsFactors=F)


## make sure to remove references (these shouldn't be in this set, but in case they are...)
referenceList=paste0( baseSampleQCDir,"/QC-Scripts/referenceList.txt" )
references = read.table(referenceList,header=F,stringsAsFactors=F)[,1]

nrefs = sum(( (kin$ID1%in%references) | (kin$ID2%in%references) ))

if(nrefs > 0 ){
    print(paste0("WARNING: ",nrefs," pairs involving reference sample are found in kinship file ",KingOutFile,". Do you expect this?" ))
    print("Removing them from output in any case...")
    kin = kin[!( (kin$ID1%in%references) | (kin$ID2%in%references) ),]
}

print(paste0(nrow(kin), " pairs of 3rd-degree (or closer) relatives found." ))

## get indexes of kinship bins
classes = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)

kin1 = classes == "dupe/twins"
kin2 = classes %in% c("sibs","parent/child")
kin2a = classes == "parent/child"
kin2b = classes == "sibs"
kin3 = classes == "2nd degree"
kin4 = classes == "3rd degree"

print(paste0(sum(kin1)," twins (or duplicates)"))
print(paste0(sum(kin2)," siblings or parent-child pairs"))
print(paste0(sum(kin2a)," parent-child pairs"))
print(paste0(sum(kin2b)," siblings"))
print(paste0(sum(kin3)," 2nd-degree relatives"))
print(paste0(sum(kin4)," 3rd-degree relatives"))

# write out list of duplicates/twins

write.table(kin[kin1,],file=paste0(OutDir,"/",OutputFile,"-duplicates-twins.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)


## bar plot of samples and numbers of related pairs.

nKins = table(c(kin$ID1,kin$ID2))
bins = c(0:10,20,30,40,100,500,Inf)
hists = sapply(2:length(bins),function(n){
    print(n)
    N = c(bins[n-1],bins[n])
    inds = names(nKins)[( nKins > N[1] )&( nKins <= N[2] )]
    total = length(inds)
    t1 = sum(inds%in%c(kin$ID1[kin1],kin$ID2[kin1]))
    t2 = sum(inds%in%c(kin$ID1[kin2],kin$ID2[kin2]))
    t3 = sum(inds%in%c(kin$ID1[kin3],kin$ID2[kin3]))
    t4 = sum(inds%in%c(kin$ID1[kin4],kin$ID2[kin4]))
    div = c(t1,t2,t3,t4)
    sum(div) == total
    return(div)
})

nKinsNot4 = table(c(kin$ID1[!kin4],kin$ID2[!kin4]))
histsNot4 = sapply(2:length(bins),function(n){
    print(n)
    N = c(bins[n-1],bins[n])
    inds = names(nKins)[( nKins > N[1] )&( nKins <= N[2] )]
    total = length(inds)
    t1 = sum(inds%in%c(kin$ID1[kin1],kin$ID2[kin1]))
    t2 = sum(inds%in%c(kin$ID1[kin2],kin$ID2[kin2]))
    t3 = sum(inds%in%c(kin$ID1[kin3],kin$ID2[kin3]))
    div = c(t1,t2,t3)
    sum(div) == total
    return(div)
})

png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
barplot(hists,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green","blue"))
legend("top",legend=rev(c("3rd","2nd","1st","twin")),col=c("gray","red","green","blue"),pch=15)
dev.off()

png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives-Log.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
barplot(hists,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green","blue"),log="y")
legend("top",legend=rev(c("3rd","2nd","1st","twin")),col=c("gray","red","green","blue"),pch=15)
dev.off()


png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives-fractions.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
hists2 = t(t(hists)/colSums(hists))
barplot(hists2,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green","blue"))
legend("top",legend=rev(c("3rd","2nd","1st","twin")),col=c("gray","red","green","blue"),pch=15)
dev.off()

# excluding 3rd degree relatives
png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives-no3rdDegree.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
barplot(histsNot4,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green"))
legend("top",legend=rev(c("2nd","1st","twin")),col=c("gray","red","green"),pch=15)
dev.off()

png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives-no3rdDegree-Log.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
barplot(histsNot4,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green"),log="y")
legend("top",legend=rev(c("2nd","1st","twin")),col=c("gray","red","green"),pch=15)
dev.off()

png(paste0(OutDir,"/",OutputFile,"-barplotNRelatives-fractions-no3rdDegree.png"),width=1000,height=1000,res=150)
par(cex.axis=0.6)
histsNot4b = t(t(histsNot4)/colSums(histsNot4))
barplot(histsNot4b,names.arg=bins[-1],xlab="number of relatives",ylab="sample count",col=c("gray","red","green"))
legend("top",legend=rev(c("2nd","1st","twin")),col=c("gray","red","green"),pch=15)

dev.off()



## plot kinship x IBS0
y = kin$Kinship
x = kin$IBS0

colors = colDef[classes]

shapes = rep(1,length(x))
shapes[kin2a] = 2

# histogram
png(paste0(OutDir,"/",OutputFile,"-kinship-hist.png"),width=1000,height=1000,res=150)
hist(y,breaks=1000,xlab="Kinship")
abline(v=c(1/2,1/4,1/8,1/16),col="red")
dev.off()

# barplot
png(paste0(OutDir,"/",OutputFile,"-kinship-classes-barplot.png"),width=1000,height=1000,res=150)
barplot(table(classes,useNA="ifany"),las=2)
dev.off()


png(paste0(OutDir,"/",OutputFile,".png"),width=1000,height=1000,res=150)

par(bty="n")
plot(NULL,ylab="Estimated kinship coefficient",xlab = "Observed proportion IBS0",xlim = c(min(x),max(x) + max(x)/5),ylim=c(min(y),max(y)),yaxt = "n")

abline(h = degrees,lty=3,lwd=0.5)
axis(2,at = degrees,labels=round(degrees,3),lty=0,las=1,cex.axis=0.8)
yText = rev(degrees) + diff(rev(degrees)[-length(degrees)])/2
text(x = rep(max(x) + max(x)/5,4),y = yText,labels=c("3rd","2nd","1st","twin"),pos=4,cex=0.8,xpd=NA)

points(x,y,col=colors,pch=shapes)

dev.off()


# plot kinship x IBS0 by ethnicity
png(paste0(OutDir,"/",OutputFile,"-byEthnicity.png"),width=1000,height=1000,res=150)

par(bty="n")
plot(NULL,ylab="Estimated kinship coefficient",xlab = "Observed proportion IBS0",xlim = c(min(x),max(x) + max(x)/5),ylim=c(min(y),max(y)),yaxt = "n")

abline(h = degrees,lty=3,lwd=0.5)
axis(2,at = degrees,labels=round(degrees,3),lty=0,las=1,cex.axis=0.8)
yText = rev(degrees) + diff(rev(degrees)[-length(degrees)])/2
text(x = rep(max(x) + max(x)/5,4),y = yText,labels=c("3rd","2nd","1st","twin"),pos=4,cex=0.8,xpd=NA)

colors=getPopColors(kin)[[1]]
shapes=getPopColors(kin)[[2]]

Order = order.by.number.occurrences(colors)
points(x[Order],y[Order],col=colors[Order],pch=shapes[Order])

dev.off()

png(paste0(OutDir,"/",OutputFile,"-Parent-child-byEthnicity.png"),width=1000,height=1000,res=150)
y1 = y; y1[!kin2a]=NA
x1 = x; x1[!kin2a]=NA 
plot(x1[Order],y1[Order],col=colors[Order],pch=shapes[Order])
dev.off()
png(paste0(OutDir,"/",OutputFile,"-Sibs-byEthnicity.png"),width=1000,height=1000,res=150)
y1 = y; y1[!kin2b]=NA
x1 = x; x1[!kin2b]=NA 
plot(x1[Order],y1[Order],col=colors[Order],pch=shapes[Order])
dev.off()
png(paste0(OutDir,"/",OutputFile,"-2ndDegree-byEthnicity.png"),width=1000,height=1000,res=150)
y1 = y; y1[!kin3]=NA
x1 = x; x1[!kin2]=NA 
plot(x1[Order],y1[Order],col=colors[Order],pch=shapes[Order])
dev.off()


## plot number of 'relatives' per sample by heterozygosity
hetFile = "../../../data/Combined/b1__b11-b001__b095-autosome-sampleqc.het.gz"
hetero = read.table(gzfile(hetFile),stringsAsFactors=FALSE,header=T)

missFile = "../../../data/Combined/b1__b11-b001__b095-autosome-sampleqc.imiss.gz"
miss = read.table(gzfile(missFile),stringsAsFactors=FALSE,header=T)

if( sum(miss$IID!=hetero$IID) > 0 ) print("WARNING: samples in het and imiss files are not in the same order...")

hetero$het = (hetero$N.NM. - hetero$O.HOM.)/hetero$N.NM.
hetero$miss = miss$F_MISS
hetero$logit.miss = log(hetero$miss/(1 - hetero$miss))

# Is there PC-corrected heterozygosity values available?
hetFilePCcorrect = "../../../data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"
if( basename(hetFilePCcorrect)%in%list.files(dirname(hetFilePCcorrect)) ) {
    load(hetFilePCcorrect,verbose=T)
    hetero=Table
}



# exclude reference samples (if necessary)
nrefs = sum(hetero$IID%in%references)
if(nrefs > 0 ){
    print(paste0("WARNING: ",nrefs," reference samples are found in heterozygosity file ",hetFile,". Do you expect this?" ))
    print("Removing them from output in any case...")
    hetero = hetero[!hetero$IID%in%references,]
}

all = c(kin$ID1[kin1],kin$ID2[kin1])
totalKin1 = table(all)

all = c(kin$ID1[kin2],kin$ID2[kin2])
totalKin2 = table(all)

all = c(kin$ID1[kin3],kin$ID2[kin3])
totalKin3 = table(all)

all = c(kin$ID1[kin4],kin$ID2[kin4])
totalKin4 = table(all)

all = c(kin$ID1[!kin4],kin$ID2[!kin4])
totalKinNot4 = table(all)


totalKin = table(c(kin$ID1,kin$ID2))



# get ethnicities 

#hetero$Batch=find.my.batch(hetero$IID,batchFile=paste0(baseSampleQCDir,"/QC-Scripts/sampleBatchList.txt"))

indexHet = match(names(totalKin),hetero$IID)
indexOther = match(names(totalKin),otherInfo$PIID)

if(sum(hetero$IID[indexHet]!=names(totalKin)) > 0) print("WARNING: something wrong with plotting of heterozygosity and kinship.")

Order = order.by.number.occurrences(otherInfo$Pops[indexOther])

y = as.vector(totalKin)
cap = 250
y[y>=cap] = cap

cols = otherInfo$Colors[indexOther]
shapes = otherInfo$Chars[indexOther]

# plot heterozygosity - all kinships
x = hetero$het[indexHet]
png(paste0(OutDir,"/",OutputFile,"-byHet.png"),width=1000,height=1000,res=150)
    plot(x[Order],y[Order],ylab="number of relatives" ,xlab="heterozygosity",col=cols[Order],pch=shapes[Order])
dev.off()

# use difference cap
y = as.vector(totalKin)
cap = 2000
y[y>=cap] = cap
png(paste0(OutDir,"/",OutputFile,"-byHet2.png"),width=1000,height=1000,res=150)
    plot(x[Order],y[Order],ylab="number of relatives" ,xlab="heterozygosity",col=cols[Order],pch=shapes[Order])
dev.off()

# plot missing rate
x = hetero$miss[indexHet]
#png(paste0(OutDir,"/",OutputFile,"-byMissing.png"),width=1000,height=1000,res=150)
plot(x[Order],y[Order],ylab="number of relatives" ,xlab="missing rate",col=cols[Order],pch=shapes[Order])
dev.off()

png(paste0(OutDir,"/",OutputFile,"-byMissing.png"),width=1000,height=1000,res=150)
x = hetero$logit.miss[indexHet]
plot(x[Order],y[Order],ylab="number of relatives" ,xlab="logit missing rate",col=cols[Order],pch=shapes[Order])
dev.off()



# use het.corrected
if("het.corrected"%in%colnames(hetero)){

    y = as.vector(totalKin)
    cap = 2000
    y[y>=cap] = cap

    x = hetero$het.corrected[indexHet]   
    png(paste0(OutDir,"/",OutputFile,"-byHet2Corrected.png"),width=1000,height=1000,res=150)
    plot(x[Order],y[Order],ylab="number of relatives" ,xlab="PC-corrected heterozygosity",col=cols[Order],pch=shapes[Order])
    dev.off()

    png(paste0(OutDir,"/",OutputFile,"-byHet2Corrected-no3rdDegree.png"),width=1000,height=1000,res=150)
    y = as.vector(totalKinNot4)
    x = hetero$het.corrected[match(names(totalKinNot4),hetero$IID)]
    Order = order.by.number.occurrences(hetero$Pops[match(names(totalKinNot4),hetero$IID)])
    cols = otherInfo$Colors[match(names(totalKinNot4),otherInfo$PIID)]
    shapes = otherInfo$Chars[match(names(totalKinNot4),otherInfo$PIID)]
    plot(x[Order],y[Order],ylab="number of relatives" ,xlab="PC-corrected heterozygosity",col=cols[Order],pch=shapes[Order])
    dev.off()
    
}


# plot heterozygosity - only 1st & 2nd degree
png(paste0(OutDir,"/",OutputFile,"-byHet-no3rdDegree.png"),width=1000,height=1000,res=150)
    y = as.vector(totalKinNot4)
    x = hetero$het[match(names(totalKinNot4),hetero$IID)]
    cols = otherInfo$Colors[match(names(totalKinNot4),otherInfo$PIID)]
    shapes = otherInfo$Chars[match(names(totalKinNot4),otherInfo$PIID)]
    plot(x,y,ylab="number of relatives" ,xlab="heterozygosity",col=cols,pch=shapes)
dev.off()



# plot legend
png(paste0(OutDir,"/",OutputFile,"-EthnicityLegend.png"),width=1000,height=1000,res=150)

plot.new()

legend("topleft",legend=ethnicities,col=ethnicity2col,pch=ethnicity2char,pt.lwd=2,horiz=F)

dev.off()



print('Finished plotting King. ')
print(paste0('Saved here: ',OutDir,"/",OutputFile,".png"))
