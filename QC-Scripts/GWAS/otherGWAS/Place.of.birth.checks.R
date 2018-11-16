h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
QCexclude=c()

args = commandArgs(trailingOnly=TRUE)

args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v14/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v14.out","genome","plots","-ymax","50","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v13/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v13-chr16.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr16.hrc+uk10k.I4.v1.1.sample","-dontComputeLD","-phased","/well/ukbiobank/expt/V2_QCed.export/data/imputation_pipeline_input/v4/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.phasing_QC.chr","-ldRData","Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v14.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-hitsAndLDcalc.RData")



dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chroms = parse.range.string(chroms)
#chroms=16

bgenData=args[which(args=="-bgenFile")+1]

# phased file?
phaseFile = "/well/ukbiobank/expt/V2_QCed.export/data/imputation_pipeline_input/v4/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.phasing_QC.chr"

if("-phased"%in%args){
    phaseFile = args[which(args=="-phased")+1]
} 





# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

if("-lreg" %in% args ) {
    useLmmInf = FALSE
    extraTitle = paste0(extraTitle,"-lreg")
} else {
    useLmmInf=TRUE
}

###  Set QC thresholds!
if("-qc" %in% args){
                                        # QC THRESHOLDS  (used in plots made before 6-07-2016)
    minmaf = 0.001
    mininfo = 0.3
    maxmiss = 0.05  # maximum 5% missing data    
} else {
    minmaf=0
    mininfo=0
    maxmiss=1
}


if("-minmaf" %in% args) minmaf = as.numeric(args[which(args=="-minmaf")+1])
if("-mininfo" %in% args) mininfo = as.numeric(args[which(args=="-mininfo")+1])
if("-maxmiss" %in% args) maxmiss = as.numeric(args[which(args=="-maxmiss")+1])



###################################
# read GWAS results!

resultsGENO = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF1 = resultsGENO$DF

resultsIMP = read.gwas(bgenData,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF2 = resultsIMP$DF


### read in recombination rates
print('Reading recombination rate files from /well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr*')

recombrates = list()
for(i in c(1:22,"X_nonPAR","X_PAR1","X_PAR2")){    
    print(i)
    r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
    recombrates[[i]] = r
}



###################################
# Read in set of snps used in phasing!

phasedBims = list()
for(i in c(1:22)){
    # NOTE: sex chromosomes not there yet!
    print(i)
    r = read.table(paste0(phaseFile,i,".bim"),header=FALSE,stringsAsFactors=FALSE)
    phasedBims[[i]] = r
}




###################################
# Plot manhattans/qqplots with rare/non-rare SNPs.

#bins=seq(0,0.5,by=0.001)
bins=c(0,0.001,0.01,0.05,0.5)
binColours = c("NA","blue","purple","red","orange")

for(chrom in 16:22){
                                        #chrom=chroms
    print(chrom)
for(i in 2:length(bins)){

    bin1 = bins[i-1]
    bin2 = bins[i]
    print(bin1);
    print(bin2);

##### 1. genotypes    
    resultsGENOsubset = resultsGENO
    resultsGENOsubset$DF = resultsGENO$DF[(resultsGENO$DF$MAF >= bin1)&(resultsGENO$DF$MAF < bin2),]
                                        # snps with MAF >= bin1 AND < bin2
    totalMarkers = sum((resultsGENO$DF$MAF >= bin1)&(resultsGENO$DF$MAF < bin2)&(resultsGENO$DF$CHR==chrom))

    print(paste0("Snps in range: ",bin1,"<= X < ",bin2,":"))
    print(totalMarkers)
    
    plot.BOLT.pvalues(Data=resultsGENOsubset,chrom,plotOutDir="plots",plotQQ=TRUE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=paste0(".MAFBin",bin1,"-",bin2))
    
    mtext(3,text=paste0(totalMarkers," markers in range: ",bin1,"<= maf < ",bin2 ),cex=3)

    dev.off()



##### 2. Imputed
    resultsIMP = read.gwas(gsub("16",chrom,bgenData),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF2 = resultsIMP$DF

    
    resultsIMPsubset = resultsIMP
    resultsIMPsubset$DF = resultsIMP$DF[(resultsIMP$DF$MAF >= bin1)&(resultsIMP$DF$MAF < bin2),]
                                        # snps with MAF >= bin1 AND < bin2
    totalMarkers = nrow(resultsIMPsubset$DF)

    print(paste0("Snps in range: ",bin1,"<= X < ",bin2,":"))
    print(totalMarkers)
    print("out of")
    print(nrow(resultsIMP$DF))
    
    plot.BOLT.pvalues(Data=resultsIMPsubset,chrom,plotOutDir="plots",plotQQ=TRUE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=paste0(".MAFBin",bin1,"-",bin2))
    mtext(3,text=paste0(totalMarkers," markers in range: ",bin1,"<= maf < ",bin2 ),cex=3)

    dev.off()

}
}


# qqplots with bins =======> GENOTYPE DATA
chroms=16:22
    
all = myQQplot(resultsGENO$DF$P[resultsGENO$DF$CHR %in% chroms])

Pvalset =  resultsGENO$Pvalset

png(paste(plotOutDir,"/",Pvalset,"-MAFBins-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
plot(NULL,xlab = expression(Expected ~ ~-log[10](italic(p))), 
     ylab = expression(Observed ~ ~-log[10](italic(p))),
     xlim=c(0,max(all$e)),
     ylim=c(0,max(all$o)))

abline(0, 1, col = "red")

for(i in 2:length(bins)){
    
    bin1 = bins[i-1]
    bin2 = bins[i]
    print(bin1);
    print(bin2);

    resultsGENOsubset = resultsGENO
    resultsGENOsubset$DF = resultsGENO$DF[(resultsGENO$DF$MAF >= bin1)&(resultsGENO$DF$MAF < bin2)&(resultsGENO$DF$CHR %in% chroms),]
    DF = resultsGENOsubset$DF
    
    qqData = myQQplot(DF$P)
    
    points(qqData$e,qqData$o,col=binColours[i],pch=1)
    
}

legend("topleft",legend=bins[-1],col=binColours[-1],pch=1,bty="n")
dev.off()



for(chrom in 16:22){
                                        #chrom=chroms
    print(chrom)
    Pvalset =  gsub("genome",chrom,resultsGENO$Pvalset)

    png(paste(plotOutDir,"/",Pvalset,"-MAFBins-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
    
    plot(NULL,xlab = expression(Expected ~ ~-log[10](italic(p))), 
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         xlim=c(0,max(all$e)),
         ylim=c(0,max(all$o)))

    abline(0, 1, col = "red")

    for(i in 2:length(bins)){
        
        bin1 = bins[i-1]
        bin2 = bins[i]
        print(bin1);
        print(bin2);

        resultsGENOsubset = resultsGENO
        resultsGENOsubset$DF = resultsGENO$DF[(resultsGENO$DF$MAF >= bin1)&(resultsGENO$DF$MAF < bin2)&(resultsGENO$DF$CHR == chrom),]
        DF = resultsGENOsubset$DF
        
        chisq1 = qchisq(DF$P,1,lower.tail=FALSE)
        lambda1 = median(chisq1)/qchisq(0.5,1)    # careful with LMM...

        qqData = myQQplot(DF$P)
        
        points(qqData$e,qqData$o,col=binColours[i],pch=1)
        
    }

    legend("topleft",legend=bins[-1],col=binColours[-1],pch=1,bty="n")
    dev.off()
}








# qqplots with bins =======> IMPUTED DATA
chroms=16:22
    
all = myQQplot(resultsIMP$DF$P[resultsIMP$DF$CHR %in% chroms])

Pvalset =  resultsIMP$Pvalset

png(paste(plotOutDir,"/",Pvalset,"-MAFBins-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
plot(NULL,xlab = expression(Expected ~ ~-log[10](italic(p))), 
     ylab = expression(Observed ~ ~-log[10](italic(p))),
     xlim=c(0,max(all$e)),
     ylim=c(0,max(all$o)))

abline(0, 1, col = "red")

for(i in 2:length(bins)){
    
    bin1 = bins[i-1]
    bin2 = bins[i]
    print(bin1);
    print(bin2);

    resultsIMPsubset = resultsIMP
    resultsIMPsubset$DF = resultsIMP$DF[(resultsIMP$DF$MAF >= bin1)&(resultsIMP$DF$MAF < bin2)&(resultsIMP$DF$CHR %in% chroms),]
    DF = resultsIMPsubset$DF
    
    qqData = myQQplot(DF$P)
    
    points(qqData$e,qqData$o,col=binColours[i],pch=1)
    
}

legend("topleft",legend=bins[-1],col=binColours[-1],pch=1,bty="n")
dev.off()



for(chrom in 16:22){
                                        #chrom=chroms
    print(chrom)
    Pvalset =  gsub("genome",chrom,resultsIMP$Pvalset)

    png(paste(plotOutDir,"/",Pvalset,"-MAFBins-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
    
    plot(NULL,xlab = expression(Expected ~ ~-log[10](italic(p))), 
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         xlim=c(0,max(all$e)),
         ylim=c(0,max(all$o)))

    abline(0, 1, col = "red")

    for(i in 2:length(bins)){
        
        bin1 = bins[i-1]
        bin2 = bins[i]
        print(bin1);
        print(bin2);

        resultsIMPsubset = resultsIMP
        resultsIMPsubset$DF = resultsIMP$DF[(resultsIMP$DF$MAF >= bin1)&(resultsIMP$DF$MAF < bin2)&(resultsIMP$DF$CHR == chrom),]
        DF = resultsIMPsubset$DF

        
        chisq1 = qchisq(DF$P,1,lower.tail=FALSE)
        lambda1 = median(chisq1)/qchisq(0.5,1)    # careful with LMM...

        qqData = myQQplot(DF$P)
        
        points(qqData$e,qqData$o,col=binColours[i],pch=1)
        
    }

    legend("topleft",legend=bins[-1],col=binColours[-1],pch=1,bty="n")
    dev.off()
}
