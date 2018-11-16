## Script to summarise cnv data for sex chromosomes.
## Output is for each sample, the mean, median, and sd of l2r for xxx width chunks of the X chromosome.
#NOTE: R library rhdf5 must be installed for this to work. This is a part of bioConductor. To install:
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-window","1000000","-out","b1__b11-b001__b095-sexchrom-sampleqc-cnv-Xchrom-summaries-females")


print(args)
h = args[-c(which(args%in%c("-in","-out","-window")),1+which(args%in%c("-in","-out","-window")))]
for(helperScript in h){
    source(helperScript)
}

windowSize = as.numeric(args[which(args=="-window")+1])
skip = windowSize/5

outFile = paste0(args[which(args=="-out")+1],"-win",windowSize/1000,"KB")

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/ExtractCNV"))
batches= all.batches()

# read in Snp-qc files
snpQC = read.SNPQC.files(type="sexchrom",justSNPs=TRUE)
QCexclude = unique(c(snpQC$arraySNPs,snpQC$concordanceSNPs))

    
# get snp annotations (position etc)
UKBB = ukbiobank.ps2snp("sexchrom")
UKBL = ukbileve.ps2snp("sexchrom")
UKBL = tbl_df(UKBL)
UKBB = tbl_df(UKBB)

ps2snp = merge(UKBL, UKBB, all = TRUE)
ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]


# read in .fam file to get list of samples
fam = read.table(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc.fam"),header=FALSE)[,2]

# only include females
females =otherInfo$PIID[ otherInfo$Inferred.Gender=="F"]
samples = intersect(fam,females)

# Just use SNPs in the sampleQC bim files (created using create-plink-subset.sh)
snpsToInclude = read.table(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc.bim"),header=FALSE)[,2]


# Any other QC snps?
sapply(snpQC,function(sn) length(intersect(sn,snpsToInclude)))

# read in other info
otherInfo = read.multiple.batch.info(batchInfoFields)

# read in cnv data and compute means and sd for all individuals. Takes about 1.5mins per batch ~= 2.5hrs!
ps2snp$Chromosome2 = names(sexChroms)[match(ps2snp$Chromosome,sexChroms)]
ps2snp$Chromosome2[ps2snp$IsInPAR1==1]="XPAR1"
ps2snp$Chromosome2[ps2snp$IsInPAR2==1]="XPAR2"
chromosomes = unique(ps2snp$Chromosome2)

#windowSize=1000000
#skip = windowSize/5

# get chunks
X = ps2snp$Chromosome%in%c(23,25)
chunkStarts = seq(min(ps2snp$Position[X]),max(ps2snp$Position[X]),by=skip)
paste0( length(chunkStarts)," chunks to compute")
allSNPs  = ps2snp[(ps2snp$AffySNPID%in%snpsToInclude) & X,]

sizes = sapply(chunkStarts,function(s) {
    sum( (allSNPs$Position >= s) & (allSNPs$Position < (s + skip) ) )
})

chunks = chunkStarts[sizes>1]
chunkSizes = sizes[sizes>1]
    
Summaries = sapply(batches,function(batch){
    
    dataSubset = read.cnv.hdf5(inds=samples,snps=allSNPs$ProbeSetID,type="log2ratio",otherInfo=otherInfo,batch=batch)

                                        # compute means over each chunk
    print("computing means...")
    Summary = sapply(chunks,function(s){
        
        theseSNPs = allSNPs$ProbeSetID[(allSNPs$Position >= s) & (allSNPs$Position < (s + skip) )]
        means = colMeans(dataSubset[[1]][theseSNPs,],na.rm=TRUE)            
        return(means)
    })
    return(Summary)
    
},simplify=FALSE)


SummariesBatches = abind(Summaries,along=1)  # this step takes a few minutes

# get mean position of snps used in each chunk (in the same order as the columns of SummariesBatches)
centres = sapply(chunks,function(s){
    theseSNPsPos = allSNPs$Position[(allSNPs$Position >= s) & (allSNPs$Position < (s + skip) )]
    centre = min(theseSNPsPos) + sum( (theseSNPsPos - min(theseSNPsPos) ))/length(theseSNPsPos)
    return(centre)
})


print("saving...")

save(SummariesBatches,centres,chunks,chunkSizes,file = paste0(baseSampleQCDir,"/data/CNV/",outFile,".RData"))
