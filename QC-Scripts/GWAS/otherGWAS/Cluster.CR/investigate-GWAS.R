library(dplyr)
library(qqman)
library(stringr)

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

h=args
for(s in h){
    source(s)
}



# get qc info
QCSNPList = read.SNPQC.files(justSNPs=TRUE)
qcexclude = unique(unlist(QCSNPList))


# ps-performance: just read one batch
psPerfFile = "/well/ukbiobank/source/V2_QCed/Probeset_Metrics/V2_All.UKBiobank.b001-b095.probeset-performance-metrics.txt"
header = t(read.table(psPerfFile,header=FALSE,nrow=1))[,1]
colclasses = rep("NULL",length(header))
colclasses[header%in%c("probeset_id","AvgFLDoverBatches")] = c("character","numeric")
psPerf = read.table(psPerfFile,header=TRUE,colClasses=colclasses)

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
snpIDs = BB.ps2snp$AffySNPID[match(psPerf$probeset_id,BB.ps2snp$ProbeSetID)]
snpIDs2 = BL.ps2snp$AffySNPID[match(psPerf$probeset_id,BL.ps2snp$ProbeSetID)]
psPerf$AffySNPID = snpIDs
psPerf$AffySNPID[is.na(psPerf$AffySNPID)] = snpIDs2[is.na(psPerf$AffySNPID)]


###############
pheno = 'Cluster.CR'
pheno = 'Cluster.CR.qnorm'
pheno = 'dQC'
#vers = 'v2'
vers = 'v3'
###############

for(pheno in c('Cluster.CR','Cluster.CR.qnorm','dQC','dQC.qnorm')){

    print(pheno)
    
    GWASOutFile = paste0(pheno,"-BOLT-LMM-",vers)
    setwd(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/',pheno))

                                        # get GWAS results
    dataFile = paste0("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/",pheno,"/BOLTLMM.",vers,"/",GWASOutFile,".out")
    DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
    DF = tbl_df(DFraw)
    if("P_BOLT_LMM_INF" %in% colnames(DF)) {
        print("using BOLT_LMM_INF")
        DF = dplyr::rename(DF, P = P_BOLT_LMM_INF) } else {
            DF = dplyr::rename(DF, P = P_LINREG)
        }
    

    colors = rep("black",nrow(DF))
    colors[DF$SNP%in%QCSNPList$batchHweSNPs] = "green"
    colors[DF$SNP%in%c(QCSNPList$plateSNPs,QCSNPList$batchSNPs)] = "purple"
    colors[DF$SNP%in%c(QCSNPList$imageSNPs)] = "orange"
    colors[DF$SNP%in%c(QCSNPList$arraySNPs)] = "red"
    colors[DF$SNP%in%c(QCSNPList$concordanceSNPs)] = "blue"

    DF$colors = colors


    # set some qc metrics and filter results
    minmaf = 0.001
    mininfo = 0.3
    maxmiss = 0.05  # maximum 5% missing data    
    chrom="genome"

    DF = dplyr::filter(DF, MAF > minmaf & F_MISS < maxmiss)

    assign(paste0("DF.",pheno),DF)

}


for(pheno in c('Cluster.CR','Cluster.CR.qnorm','dQC','dQC.qnorm')){

    print(pheno)

    GWASOutFile = paste0(pheno,"-BOLT-LMM-",vers)
    setwd(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/',pheno))

    DF = get(paste0("DF.",pheno))


                                        # set plot names
    Pvalset = paste0(GWASOutFile,".chr",chrom,".maf",minmaf,".miss",maxmiss,".pruned")


                                        # plot allele freq by pvalue for highest SNPs
    subset = (DF$P < 5e-8) & (!DF$SNP %in% qcexclude)
    x = DF$MAF[subset]
    xmiss = DF$F_MISS[subset]
    y = DF$P[subset]
    y[y < 10^-50] = 10^-50
    cols = DF$colors[subset]
    Order = order.by.number.occurrences(cols)

    png(paste0("plots/",Pvalset,"-MAFbyPval.png"),width=1000,height=1000,res=150)
    plot(x[Order],-log10(y[Order]),xlab="MAF",ylab="-log10(p)",col=cols[Order])
    dev.off()
    png(paste0("plots/",Pvalset,"-MAFhist.png"),width=1000,height=1000,res=150)
    hist(x,xlab="MAF",col=add.alpha("red",0.5),xlim=c(0,max(DF$MAF)),freq=FALSE,breaks=50,main="minor allele frequency")
    hist(DF$MAF[!subset],xlab=NA,col=add.alpha("blue",0.5),add=TRUE,freq=FALSE,breaks=50)
    dev.off()
    png(paste0("plots/",Pvalset,"-CRhist.png"),width=1000,height=1000,res=150)
    hist(DF$F_MISS[!subset],xlab="SNP missing rate",col=add.alpha("blue",0.5),xlim=c(0,max(DF$F_MISS)),freq=FALSE,breaks=50,main="missing rates")
    hist(xmiss,xlab=NA,col=add.alpha("red",0.5),add=TRUE,freq=FALSE,breaks=50)
    dev.off()

                                        # does it correlate with FLD on SNPs?

    y = DF$P[subset]
    x = psPerf$AvgFLDoverBatches[match(DF$SNP[subset],psPerf$AffySNPID)]
    x2 = psPerf$AvgFLDoverBatches[match(DF$SNP[!subset],psPerf$AffySNPID)]
    y2 = DF$P[!subset]
    cols2 = DF$colors[!subset]
    Order2 = order.by.number.occurrences(cols2)

    png(paste0("plots/",Pvalset,"-AvgFLDbyPval.png"),width=1000,height=1000,res=150)
    plot(x[Order],-log10(y[Order]),col=cols[Order],xlab="AvgFLDoverBatches - UKBiobank",ylab="-log10(p)")
    dev.off()
    png(paste0("plots/",Pvalset,"-AvgFLDhist.png"),width=1000,height=1000,res=150)
    hist(x,xlab="AvgFLDoverBatches",col=add.alpha("red",0.5),xlim=c(0,max(psPerf$AvgFLDoverBatches,na.rm=TRUE)),freq=FALSE,breaks=50,main="Average FLD over batches (UKBiobank)")
    hist(x2,xlab=NA,col=add.alpha("blue",0.5),add=TRUE,freq=FALSE,breaks=50)
    dev.off()

                                        # Save table of SNPs above 10^-8
    write.table(DF[subset,c("SNP","CHR","BP","P")],file=paste0(GWASOutFile,"-significantSnps.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE )


}


cluster.plots <- function(pheno,nsnps){
    print(pheno)

    GWASOutFile = paste0(pheno,"-BOLT-LMM-",vers)
    wd = paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/',pheno)

    DF = get(paste0("DF.",pheno))

                                        # set plot names
    Pvalset = paste0(GWASOutFile,".chr",chrom,".maf",minmaf,".miss",maxmiss,".pruned")

                                        # plot allele freq by pvalue for highest SNPs
    subset = (DF$P < 5e-8) & (!DF$SNP %in% qcexclude)

                                            # select some SNPs to have a look at
    SNPs = DF$SNP[order(DF$P)]
    SNPs = SNPs[!SNPs%in%qcexclude]
    SNPs = SNPs[1:nsnps]

    write.table(SNPs,file=paste0(wd,"/",GWASOutFile,"-top",nsnps,"Snps.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE )

                                        # plot cluster plots for the top n snps
   # forSystem = paste0("Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",wd,"/",GWASOutFile,"-top",nsnps,"Snps.txt ",wd,"/",GWASOutFile,"-clusterPlots-top",nsnps,"-orig -b Batch_b001,Batch_b064,Batch_b092,Batch_b095,UKBiLEVEAX_b7 -orig")

   # system(forSystem)

    return(NULL)
}

for(pheno in c('Cluster.CR','Cluster.CR.qnorm','dQC','dQC.qnorm')){
    cluster.plots(pheno,20)
}



