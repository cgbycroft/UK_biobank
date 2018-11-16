source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')

library(dplyr)
library(qqman)
library(stringr)

sexChroms = c(23,24,25,26)
names(sexChroms) = c("X","Y","XY","MT")

#args = commandArgs(TRUE)
#dataFile = args[1]

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

for(s in h){
    source(s)
}


dataFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-quant-Age-v1.out"

outname = gsub(".out","",basename(dataFile))

######## Read in some SNP-QC stats
hwe = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc.hwe",stringsAsFactors=FALSE,header=TRUE)
hwe = tbl_df(hwe)
hwe = filter(hwe,TEST=="ALL")

miss = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc.imiss",stringsAsFactors=FALSE,header=TRUE)
miss = tbl_df(miss)

maf = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc.frq",stringsAsFactors=FALSE,header=TRUE)
maf = tbl_df(maf)

# image artefact SNPs
art = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/SelectSNPs/probesets_in_subimage_37555.txt",stringsAsFactors=FALSE,header=TRUE)[,2]


######## Read in GWAS results

if(chroms=="all") chroms = 1:22 else chroms = parse.range.string(chroms)

if(!grepl("%%",dataFile)){
    DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
}


######## Set usual filtering

DF = dplyr::tbl_df(DFraw)

if("P_BOLT_LMM_INF" %in% colnames(DF)) {
    print("using BOLT_LMM_INF")
    DF = dplyr::rename(DF, P = P_BOLT_LMM_INF) } else {
    DF = dplyr::rename(DF, P = P_LINREG)
}

# get hwe into D
hwe2 = select(hwe,SNP,P); colnames(hwe2)[2] = "HWE"
DF = left_join(DF,hwe2,by=c("SNP"="SNP"))


######## What are the properties of the highest SNPs?
high = DF$P < 5e-8

png(paste0( "plots/",outname,"-snpProperties-pval8-%02d.png" ),height=500,width=500,res=100 ) 
for (var in c("F_MISS","MAF","HWE")){

    x = DF[[var]][high]
    nSNPS = length(x)
    hist(x,breaks=nSNPS/20,xlab=var)

}
dev.off()


######## Apply some filtering ==> what's left??


minmaf = 0.001
mininfo = 0.3
maxmiss = 0.05  # maximum 5% missing data
minhwe = 10^-20
    
DF.Filt = dplyr::filter(DF, MAF > minmaf & F_MISS < maxmiss & HWE > minhwe)

high.Filt = DF.Filt$P < 10^-8
sum(high.Filt)


######### Run a plot
Pvalset = paste(outname,".maf",minmaf,".miss",maxmiss,".hwe",-log10(minhwe),".pruned",sep="")

Ymax = ceiling(max(-log10(DF.Filt$P[DF.Filt$P!=0]),na.rm=T)) + 10
Ymax = min(Ymax,50)

png(paste("plots/",Pvalset,"-manhattan%02d.png",sep=""),width=100,height=12,units="in",res=150)
par(las=1,font.main=1,cex.axis=2,cex.lab=2,mar=c(7 ,7, 5 ,2))
myManhattan(DF.Filt,ymax=Ymax)
dev.off()


summary(DF.Filt)
