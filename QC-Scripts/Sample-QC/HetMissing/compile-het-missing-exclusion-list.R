## Script to plot (raw) heterozygosity and missingness
library(dplyr)

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss")

print(args)

h = args[-c(which(args%in%c("-in","-out","-npcs")),1+which(args%in%c("-in","-out","-npcs")))]
for(helperScript in h){
    source(helperScript)
}


input = args[which(args=="-in")+1]
outname = basename(input)

# read in aberrant results from BIWI analysis
load(paste0(outname,"-BIWI-lambda120.RData"),verbose=T)

# read in het-missing stats
load(paste0(input,".RData"),verbose=TRUE)


# Set missing threshold for all samples
######################
missThreshold = logit(0.05)
#missThreshold = logit(0.06)  # this was the interim release threshold
######################

# apply threshold across-the-board
exclude.miss.samples = Table$IID[Table$logit.miss > missThreshold]


# get list of BIWI outliers. "group" is aberrant's list of outliers.
mean.het.correct = mean(aberrant.hetmiss$x[,"het.corrected"],na.rm = TRUE)

aberrant.outlier = (aberrant.hetmiss$group == 1 & aberrant.hetmiss$x[,"het.corrected"] > mean.het.correct) | 
		(aberrant.hetmiss$x[,"logit.miss"] > missThreshold)

# this is the order of samples in the aberrant analysis
ToUse.aberrant = dplyr::filter(Table, Pops %in% c("British","Irish","Any other white background","Indian") )
aberrant.outlier.samples = ToUse.aberrant$IID[aberrant.outlier]


# get all unique samples
hetmiss.outliers.all = unique(c(exclude.miss.samples,aberrant.outlier.samples))

Table$hetmiss.outlier = FALSE
Table$hetmiss.outlier[Table$IID%in%hetmiss.outliers.all] = TRUE

table(Table$Pops,Table$hetmiss.outlier)
table(Table$hetmiss.outlier)


# save a list of these samples
write.table(hetmiss.outliers.all,file=paste0(outname,"-all-hetmiss-outliers.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)


## plot these (red and black)
Order = c(which(!Table$hetmiss.outlier),which(Table$hetmiss.outlier))
png(paste0("plots/",outname,"-het.missing-",round(missThreshold,3),"-all-outliers-logitScale.png"),height=1000,width=1000,res=150)
plot(Table$logit.miss[Order],Table$het.correct[Order],col = Table$hetmiss.outlier[Order] + 1,xaxt="n",
     xlab="missing rate (logit scale)",ylab="PC-corrected heterozygosity",
     main=paste0(length(hetmiss.outliers.all)," outliers") )
axis(1,at = -9:-2,labels=round(inv.logit(c(-9:-2)),4))
abline(v = missThreshold,col="blue",lty=3)
dev.off()

png(paste0("plots/",outname,"-het.missing-",round(missThreshold,3),"-all-outliers-rateScale.png"),height=1000,width=1000,res=150)
plot(Table$miss[Order],Table$het.correct[Order],col = Table$hetmiss.outlier[Order] + 1,
     xlab="missing rate",ylab="PC-corrected heterozygosity",
     main=paste0(length(hetmiss.outliers.all)," outliers") )
abline(v = inv.logit(missThreshold),col="blue",lty=3)
dev.off()

