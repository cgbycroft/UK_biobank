#####################
# Script to computing posteriors for a given GWAS at a given set of regions.
#####################

#####################
# Preliminaries

args = commandArgs(trailingOnly=TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

# args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out","-regions","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts","-chr","19","-prior","0.2")

# args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz","-regions","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts","-chr","2")

print(args)

dataFile = args[which(args=="-gwas")+1]
regionFile = args[which(args=="-regions")+1]
outDir = dirname(dataFile)

if("-outdir"%in%args){
    outDir = args[which(args=="-outdir")+1]
}
print(outDir)

if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chrom = as.numeric(args[which(args=="-chr")+1])
print(paste0("Processing chromosome ",chrom))

if(grepl("%%",dataFile)) dataFile = gsub("%%",chrom,dataFile)

#####
 # What is the prior for the bayse factors?

prior = as.numeric(args[which(args=="-prior")+1])

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


## Read in region file
regions = read.table(regionFile,header=TRUE)

regions = regions[regions[,1]==chrom,]

nRegions = nrow(regions)

if(nRegions==0){
    print("No regions actually found for this chromosome. Not reporting anything.")
    quit()
}

print("This many regions: ")
print(nRegions)

# NOTE: Some regions might overlap. This is because of the region definition.
overlappingRegions = apply(regions,1,function(x) which(! ( (regions[,"end"] < x["start"]) | (regions[,"start"] > x["end"]) )) )

overlaps = overlappingRegions[which(sapply(overlappingRegions,function(x) length(x)>0))]

print(paste0(length(overlaps)," regions actually overlap one-another!"))



## Read in bayes factors
bfFile = paste0(dataFile,".",prior,".bf")

## Read in GWAS results (just to get the order of the SNPs)
results = read.gwas(dataFile,chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=50,QCexclusions=c(),extraTitle=extraTitle,useLmmInf=TRUE,bayesFactorFile=bfFile)

if(is.null(results)){
    print("No GWAS data available for this chromosome.")
    quit()
}

DF = results$DF
DF = DF[DF$CHR==chrom,]

# Exclude any markers with MAF=0 (they won't have a Bayes Factor).
noBeta = is.na(DF[,paste0("bf.",prior)])
DF = DF[!noBeta,]

print(paste0(sum(noBeta)," markers on this chromosome filtered out because no BF (MAF=0)."))
print(paste0(nrow(DF)," (filtered) markers in this chromosome."))

## Subset by the regions

DFsubset = lapply(1:nRegions,function(i){
    r = regions[i,]
    #print(r)
    these = (DF$BP >= r[,2])&(DF$BP <= r[,3])
    if(sum(these)==0) {
        print(paste0("No data found for the region: ",r[,2]," to ",r[,3]))
        return(NULL)
    }
    out = cbind(i,DF[these,])
    colnames(out)[1]="region.index"
    return(out)
})

DFsubset = abind(DFsubset,force.array=FALSE,along=1)

print(head(DFsubset))

print(paste0(nrow(DFsubset)," SNPs found in these regions to compute posteriors."))
print(paste0(length(unique(DFsubset$SNP2))," of these are unique (due to overlapping regions)."))

noSNPs = which( !1:nRegions %in% unique(DFsubset[,1]))
noSignifSNPs = which( !1:nRegions %in% unique(DFsubset[DFsubset$P < signifLevel,1]))

print(paste0(length(noSNPs)," regions have no snps in them."))
print(paste0(length(noSignifSNPs)," regions have no significant snps in them, at ",signifLevel))


#######
# Compute posteriors for each region

# First compute the totals
#DFsubsetRaw=DFsubset
#DFsubset=DFsubset[sample(1:nrow(DFsubset),nrow(DFsubset),replace=FALSE),]

posteriors = lapply(1:nRegions,function(i){
    in.region = DFsubset$region.index == i
    BF = DFsubset[in.region,paste0("bf.",prior)]
    post = compute.BF.posteriors(BF)
    out = cbind(which(in.region),post)
    return(out)
})

Posteriors = abind(posteriors,along=1,force.array=FALSE)

print(head(Posteriors))

DFsubset$Posterior = NA
DFsubset$Posterior[Posteriors[,1]] = Posteriors[,2]

print( head(DFsubset) )

####### save the relevant files
outFile = paste0(outDir,"/Posteriors-",results$Pvalset,".",prior,".RData")
save(DFsubset,regions,chrom,regionFile,dataFile,overlaps,file=outFile)


print("DONE! see:")
print(outFile)

