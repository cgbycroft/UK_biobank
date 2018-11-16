## Script to plot output of submit-pca-UKBio.sh. This will be projections based on data from each batch.
library(zoo)
library(grid)
args=commandArgs(trailingOnly=T)

args = c('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R','/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R','-in','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British.snpload.map','-out','b1__b11-b001__b095-pca-UKbio-White_British')


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

Info=read.multiple.batch.info(batchInfoFields)
Info = dplyr::tbl_df(Info)

Rfile=paste0(baseSampleQCDir,"/data/PCA/",OutputFile,".RData")
load(Rfile,verbose=T)
nPCs = sum(grepl("PC",colnames(PCs)))


# White British subset
wb = paste0(baseSampleQCDir,"/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt")
wbList = read.table(wb,header=FALSE)[,1]

# make RData file with only WB set
PCs = PCs[PCs$PIID%in%wbList,]
save(PCs,file=paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-White_British_subset.RData") )


# get mapping shape files etc.
load(mapFileUK0,verbose=TRUE)
uk0@bbox = rbind(c(-10000,665695),c(-10000,1325640))
# original coordinates
#          min       max
#x -296677.815  655695.6
#y    5408.466 1295640.3

baseMap = spplot(uk0,zcol="ID_0",col.regions="transparent",colorkey=FALSE,main="hello")

png("plots/test.png",height=1500,width=1000,res=150)
print(baseMap)
dev.off()

# plot on map
x = Info$Place.of.birth.in.UK...east.co.ordinate[match(PCs$PIID,Info$PIID)]
y = Info$Place.of.birth.in.UK...north.co.ordinate[match(PCs$PIID,Info$PIID)]    
notNA = (!is.na(x))&(!is.na(y))
colors = c("blue","purple","red")
mapHeight = 0.7

for(i in 1:10){
    png(paste0(OutDir,'/plots/',OutputFile,'-POB-UK-pc',i,'.png'),height=1500,width=1000,res=150)
    pc = paste0("PC",i)
    cols = colour.scale(PCs[[pc]],colourSet = colors)
    Points = make.points(x[notNA],y[notNA],col=cols[notNA],pch=16)
    map = update(baseMap,sp.layout=list(Points),main=)
    plot.new()
    grid.newpage()
    pushViewport(viewport(y=1,height=mapHeight,width=1,just="top"))
    print(map)
    upViewport(1)
    pushViewport(viewport(y=0,height=(1-mapHeight),width=1,just="bottom"))    
    makeScaleWithHist(PCs[[pc]][notNA],colourSet = colors)
    dev.off()
}


