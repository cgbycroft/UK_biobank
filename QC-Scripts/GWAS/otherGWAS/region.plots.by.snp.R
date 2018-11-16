##### This script requires input files from region.plots

source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


genotypeFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc"


args = commandArgs(TRUE)

#args = c("Array.as.binary-BOLT-LMM-v11a-23.out.chr23.maf0.miss1.pruned-recombRegions-raw-hitsAndLDcalc.RData","allsnps","23")

print(args)

hitsFile = args[1]
snpFile = args[2]
chrom = args[3]

load(hitsFile,verbose=TRUE)
HITS = HITSchrom[[chrom]]
DF = HITS$DFthese
ld = LDchrom[[chrom]]
    
if( grepl(".txt",snpFile)) snps = read.table(snpFile,header=FALSE,stringsAsFactors=FALSE)[,1]

if( snpFile == "allsnps" ) snps = DF$SNP

if( "-snps" %in% args ){
    snps = unique( str_split(args[which(args=="-snps")+1],",")[[1]] )
}

if( "-range" %in% args ){
    Range = args[which(args=="-range")+1]
    sf = parse.range.string(Range)
    snps = DF$SNP[( DF$BP >= sf[1] ) & ( DF$BP <= sf[2] )]
}


snps = DF$SNP[DF$SNP%in%snps]

if( length(snps) == 0 ) quit() else print( paste0( length(snps)," SNPs to plot regions for.") ) 


# read recombination maps
sex = FALSE
if( ("-sex" %in% args) | (chrom %in% c("23","24","25","26")) ) sex =TRUE
recombrates = read.recomb(sex,chrom)

if( sex ) type = "sexchrom" else type = "autosome"
BB.ps2snp = ukbiobank.ps2snp(type)
BL.ps2snp = ukbileve.ps2snp(type)


Pvalset = gsub("-hitsAndLDcalc.RData","",hitsFile)
    
# run region.plots for specified snps
recomb = recombrates[[chrom]]
sapply(1:length(snps),FUN=function(s) {
    
    pos = DF$BP[DF$SNP==snps[s]]
    r = which( (HITS$regs[,1] <= pos) & (HITS$regs[,2] >= pos) )[1] # just pick the first hit region that this snp appears in.    
    regionToPlot = c(HITS$regs[r,1],HITS$regs[r,2])
    plot.gwas.region(DF,snps[s],catFile=NULL,ld,plotDir=paste0("plots/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=50,width=10^6,moreToPlot=TRUE,Pvalset=Pvalset)

    outWidth = diff(regionToPlot)
    abline(v=regionToPlot,col="darkgray",lty=3,lwd=3)
    text(regionToPlot[1],y=par("usr")[4],labels=paste0(round(outWidth/1000000,1),"MB region"),xpd=NA,cex=2,pos=3)
    dev.off()
})
