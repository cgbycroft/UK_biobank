##########
# Script to define a set of regions for comparing results across studies, using results from a given study.
##########
# Output is a regions file: chrom start end

#####################
# Preliminaries

args = commandArgs(TRUE)

h=c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


# args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out.signif","-chr","1-22")

# args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz.signif","-chr","19-22")


print(args)

dataFile = args[which(args=="-gwas")+1]
chromsString = args[which(args=="-chr")+1]

CM_RANGE=0.125
KB_BUFFER=25

if("-width" %in% args ) CM_RANGE = as.numeric(args[which(args=="-width")+1])
if("-buffer" %in% args ) KB_BUFFER = as.numeric(args[which(args=="-buffer")+1])

outDir = dirname(dataFile)
if("-outdir"%in%args){
    outDir = args[which(args=="-outdir")+1]
}
print(outDir)

if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chroms = parse.range.string(chromsString)


#####
# SET QC SUBSETTING

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


## Read in GWAS results
results = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=50,QCexclusions=c(),extraTitle=extraTitle,useLmmInf=TRUE)

DF = results$DF

outFile = paste0(outDir,"/regions-",CM_RANGE,"cM-",KB_BUFFER,"KB-",results$Pvalset,".txt")


### read in recombination rates
recombrates = read.recomb(chrom=chromsString)

## Find the hits and the region around it

HITSchrom = vector("list", 22)
OUT = vector("list",22)

nRegions = 0
for(i in chroms){
    print(paste0("chr:",i))
    
    recomb = recombrates[[as.character(i)]]
    HITS = find.hits.recomb(DF,i,minGap=CM_RANGE,buffer=KB_BUFFER,recombRates=recomb)
    HITSchrom[[as.character(i)]] = HITS
    OUT[[as.character(i)]] = cbind(i,HITS$regs,HITS$topVariants,HITS$topVariants2)
}


pos = abind(OUT,along=1,force.array=FALSE)
colnames(pos)[1] = "chr"
posF = format(pos,scientific=FALSE,trim=TRUE)

    
nRegions = nrow(pos)

write.table(posF,file=outFile,quote=FALSE,row.names=FALSE,col.names=TRUE)


print("DONE! See:")
print(outFile)
print(paste0("Should be ",nRegions," regions in total."))
