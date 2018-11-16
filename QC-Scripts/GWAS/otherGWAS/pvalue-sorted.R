h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


args = commandArgs(trailingOnly=TRUE)

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Rnorm.0.1/BOLTLMM.v17/Rnorm.0.1-BOLT-LMM-v17-chr%%.out")

print(args)

inData = args[1]
outFile = gsub(".out",".pvalues.RData",gsub("%%","genome",basename(inData)))


pvals = list()

for( i in 1:22 ){

    print(i)
    RnormFileImputation=paste0(gsub("%%",i,inData))

    RnormResults = read.gwas(RnormFileImputation,chrom=i)

   # sort by pvalue!
    psortedLREG = sort(RnormResults$DFraw$P_LINREG)
    psortedLMM = sort(RnormResults$DFraw$P_BOLT_LMM_INF)

    pvals[[as.character(i)]] = list("psortedLREG"=psortedLREG,"psortedLMM"=psortedLMM)
}

save(pvals,file=outFile)

print(paste0("DONE! ",outFile ))
