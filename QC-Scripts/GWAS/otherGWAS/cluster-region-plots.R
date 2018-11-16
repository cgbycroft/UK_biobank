source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


args = commandArgs(TRUE)

#args = c("<pheno>","-vers","-outdir","-snps","-rsids")
print(args)
vers = args[which(args=="-vers")+1]
pheno = args[1]
setwd(paste0(baseSampleQCDir,"/QC-Scripts/GWAS/otherGWAS/",pheno))
outdir = paste0(getwd(),"/",args[which(args=="-outdir")+1])
snps = args[which(args=="-snps")+1] # a comma-separated list of SNPs, or a file of snps.

if("-rsids"%in%args){
                                        #read in annotation files (only necessary if snp region plots are in rsid form)
    print("reading in snp annotation...")
    BB.ps2snp = ukbiobank.ps2snp("autosome")
    BL.ps2snp = ukbileve.ps2snp("autosome")
    BB.ps2snpSex = ukbiobank.ps2snp("sexchrom")
    BL.ps2snpSex = ukbileve.ps2snp("sexchrom")

}


if( grepl(".txt",snps) ) {
    snpList = read.table(snps,header=FALSE,stringsAsFactors=FALSE)[,1]
} else {
    snpList = sapply(snps,function(x) str_split(x,",")[[1]])
}

otherClusterPlotCommands = c()

if( "-otherCommands" %in% args ){
    otherClusterPlotCommands = args[( which(args == "-otherCommands")+1):length(args)]
    otherClusterPlotCommands = paste(otherClusterPlotCommands, collapse=" ")
}


# which snps don't have cluster plots?
hasClusterPlots = sapply(snpList,function(s) length(list.files(path=outdir,pattern=paste0("clusterplot.*",s),recursive=TRUE) > 0))

snps2 = paste(snpList[!hasClusterPlots],collapse=",")
print(snps2)

if( sum(!hasClusterPlots) > 0){
    print(paste0('making cluster plots for ',sum(!hasClusterPlots),' snps.'))
    forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps2,' ',outdir,' ',otherClusterPlotCommands)
    system(forSystem)
}


# combine with region plots
print('combining with region plots')
print(snpList)

if("-prefix" %in% args) Prefix = args[which(args=="-prefix")+1] else  Prefix=paste0(pheno,".*",vers,".*")

sapply(snpList,combine.cluster.regions.plots,prefix=Prefix,clusterPlotDir=outdir)

