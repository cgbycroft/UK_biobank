# Script to process and plot window-based-averages of l2r in females

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-window","1000000","-in","b1__b11-b001__b095-sexchrom-sampleqc-cnv-Xchrom-summaries-females")


print(args)
h = args[-c(which(args%in%c("-in","-out","-window")),1+which(args%in%c("-in","-out","-window")))]
for(helperScript in h){
    source(helperScript)
}


windowSize = as.numeric(args[which(args=="-window")+1])
skip = windowSize/5

outFile = paste0(args[which(args=="-in")+1],"-win",windowSize/1000,"KB")

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/sexChroms"))
batches= all.batches()

inputData = paste0(baseSampleQCDir,"/data/CNV/",outFile,".RData")


# use same criteria for snps and sample inclusions as 'gender-checks.R'

# read in Snp-qc files
snpQC = read.SNPQC.files(type="sexchrom",justSNPs=FALSE)

# get snp annotations (position etc)
UKBB = ukbiobank.ps2snp("sexchrom")
UKBL = ukbileve.ps2snp("sexchrom")
UKBL = tbl_df(UKBL)
UKBB = tbl_df(UKBB)

ps2snp = merge(UKBL, UKBB, all = TRUE)
ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]
ps2snp$Chromosome2 = names(sexChroms)[match(ps2snp$Chromosome,sexChroms)]
ps2snp$Chromosome2[ps2snp$IsInPAR1==1]="XPAR1"
ps2snp$Chromosome2[ps2snp$IsInPAR2==1]="XPAR2"

# duplicates list
dupFile = "/well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt"
dup = read.table(dupFile,header=FALSE)[,1]

# get other Info
otherInfo = read.multiple.batch.info(batchInfoFields)
otherInfo = tbl_df(otherInfo)
otherInfoRaw = otherInfo

# get het/miss outliers
hetmissFile = paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt")
hetmissOutliers = read.table(hetmissFile,header=FALSE,stringsAsFactors=FALSE)[,1]

# exclude these from the analysis
otherInfo = filter(otherInfo,!PIID%in%hetmissOutliers)

# male/females
males = otherInfo$Inferred.Gender=="M"
females = otherInfo$Inferred.Gender=="F"

# mismatches
mismatches = otherInfo$Inferred.Gender!=otherInfo$Submitted.Gender
mismatchesM = mismatches & males
mismatchesF = mismatches & females

# duplicates (and exclude)
dupes = otherInfo$PIID%in%dup
otherInfo = otherInfo[!dupes,]


# get data from cnv summaries
load(paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries.RData"),verbose=TRUE)

# order by otherInfo samples
bafsum = outData[["baf"]][otherInfo$PIID,]
l2rsum = outData[["log2ratio"]][otherInfo$PIID,]

sum( rownames(bafsum)!=otherInfo$PIID ) 
sum( rownames(l2rsum)!=otherInfo$PIID )

# get list of SNPs used in cnv summaries
snpsToInclude = read.table(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc.bim"),header=FALSE)[,2]
nsnps = sapply(names(sexChroms2),function(j) length(ps2snp$ProbeSetID[(ps2snp$AffySNPID%in%snpsToInclude) & (ps2snp$Chromosome2==j) ]))


# set up colors for plotting (MUST APPLY THESE IN THE RIGHT ORDER. I.E Duplicates last!)
colDefinition = c("green","purple","blue","red","black"); names(colDefinition)=c("Submitted and inferred female","Submitted and inferred male","Submitted male and inferred female","Submitted female and inferred male","Mismatch")
colors = rep(NA,nrow(otherInfo))
colors[females] = colDefinition["Submitted and inferred female"]
colors[males] = colDefinition["Submitted and inferred male"]
#colors[mismatches] = colDefinition["Mismatches"]
colors[mismatchesM] = colDefinition["Submitted female and inferred male"]
colors[mismatchesF] = colDefinition["Submitted male and inferred female"]
colors[dupes] = "transparent" # just plot duplicates transparent

Order = order.by.number.occurrences(colors)

# ethnicity colors
colorsE = otherInfo$Colors
shapes = otherInfo$Chars
#colors[mismatches] = "black"
OrderE = order.by.number.occurrences(colorsE)


######################### End of preliminaries

# load in mean data (computed in ../ExtractCNV/compute-l2r-means-across-X.R )
load( inputData, verbose = TRUE)


######################### Run linear model on chunked data

means = SummariesBatches[rownames(SummariesBatches)%in%otherInfo$PIID,]
otherInfoIndex =  match(rownames(means), otherInfo$PIID) 
x = otherInfo$Age.when.attended.assessment.centre[otherInfoIndex]

# this takes about 15 seconds per 100 chunks ==> about 10mins for 4000 chunks
lms = apply(means,2,function(y){
    t = lm( y ~ x )
    return(t)    
})

# this is fast
betas = sapply(lms,function(l){
    coef(l)[2]
})

# this is a bit slower
pvals = sapply(lms,function(l){
    ls = summary(l)
    pf(ls$fstatistic[1], ls$fstatistic[2], ls$fstatistic[3],
       lower.tail = FALSE)  # not sure if lower tail is what I want here...
})

# save results (just betas and pvals)
save(betas,pvals,file=paste0(outFile,"-lm-test-results.RData"))


######################### plot results
chunkCols = colour.scale(log10(chunkSizes))

x = centres
y = -log10(pvals)
png( paste0("plots/",outFile,"-lm-test-pvalues.png") ,width=5000,height=1000,res=150)
plot(x,y,xlab="Position centre for mean",col=chunkCols,pch=6)
abline(v = c(par1,par2),lty=3)
dev.off()

y = betas
png( paste0("plots/",outFile,"-lm-test-betas.png") ,width=5000,height=1000,res=150)
plot(x,y,xlab="Position centre for mean")
abline(v = c(par1,par2),lty=3)
abline(h = 0,lty=3)
dev.off()

