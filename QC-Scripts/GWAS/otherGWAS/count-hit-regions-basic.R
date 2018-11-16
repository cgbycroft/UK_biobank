#####################
# Script to count how many Mb windows contain at least one significant hit.
#####################

#####################
# Preliminaries

args = commandArgs(TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

# args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16","-chr","2","-title","1Mb-regs")

# args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt","-chr","2","-title",".1Mb-regs")

print(args)

gwasFile = args[1] # This is the GWAS results
regionFile = args[2]
outDir = dirname(gwasFile)

if("-outdir"%in%args){
    outDir = args[which(args=="-outdir")+1]
}
print(outDir)

if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chrom = as.numeric(args[which(args=="-chr")+1])
#chroms = parse.range.string(args[which(args=="-chr")+1])


#####
 # SET QC SUBSETTING

if("-qc" %in% args){
# DEFAULT QC THRESHOLDS  (used in plots made before 6-07-2016)
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

#plotOutDir = paste0(baseSampleQCDir,"/QC-Scripts/GWAS/otherGWAS/plots")


## Read in region file
regions = read.table(regionFile,header=TRUE)

regions = regions[regions[,1]==chrom,]

nRegions = nrow(regions)
print(nRegions)


if(nRegions==0){
    print("No regions actually found for this chromosome. Not reporting anything.")
    quit()
}

## Read in GWAS results
results = read.gwas(gwasFile,chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=50,QCexclusions=c(),extraTitle=extraTitle,useLmmInf=TRUE)

DF = results$DF

# subset for this chromosome, and snps over significance threshold
DF = DF[(DF$CHR==chrom)&(DF$P<signifLevel),]
   
# count how many in each region
counts = apply(regions,1,function(r){
    print(r)
    these = (DF$BP >= r[2])&(DF$BP <= r[3])
    tot = sum(these)
    max_INFO=NA
    max_MAF=NA
    if(tot>0){
        if("INFO" %in% colnames(DF))  max_INFO = max(DF$INFO[these])
        if("MAF" %in% colnames(DF))  max_MAF = max(DF$MAF[these])
    }
    return(c("counts"=tot,"max_INFO"=max_INFO,"max_MAF"=max_MAF))
})
counts = t(counts)

out = cbind(regions,counts)

#outFileName = paste0("hitCounts-",signifLevel,".txt")
outFileName = paste0(results$Pvalset,"-",signifLevel,".txt")

write.table(out,file=paste0(outDir,"/",outFileName),quote=FALSE,col.names=TRUE,row.names=FALSE)

print("DONE! See:")
print(paste0(outDir,"/",outFileName))
