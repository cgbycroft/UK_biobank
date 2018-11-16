## Script to plot X/Y intenisty data in order to check sex determination by Affy

source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-out","b1__b11-b001__b095-sexchrom-sampleqc-sexCheck")

print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/sexChroms"))

outFile = args[which(args=="-out")+1]

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

# get PAR individuals missing rate
PARmiss = read.table(gzfile(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc-25.imiss.gz")),header=TRUE,stringsAsFactors=FALSE)
PARmiss$logit.miss = log( PARmiss$F_MISS/(1 - PARmiss$F_MISS) )

# get Y individuals missing rate
Ymiss = read.table(gzfile(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc-24.imiss.gz")),header=TRUE,stringsAsFactors=FALSE)
Ymiss$logit.miss = log( Ymiss$F_MISS/(1 - Ymiss$F_MISS) )

# get X het rate
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc-23-hetcorrected-6PCs-imiss.RData"),verbose=TRUE)
XTable = Table

apply.het.correction.6PCs = function(hetmisspcs,het.correct) {
  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4","PC5","PC6") %in% names(hetmisspcs))) {
    return(NULL)
  }
  het = hetmisspcs$het
  lmfit = het.correct$lmfit
  betahat = het.correct$betahat
  Intercept = het.correct$Intercept
  ## het.corrected = het - (model.matrix(lmfit) %*% betahat - Intercept)
  het.corrected = het - (predict(lmfit, newdata = hetmisspcs) - Intercept)
  return (het.corrected)
}
Xhet.corrected = do.call( paste0("apply.het.correction.6PCs"), list(XTable,het.correct) )
Xhet.corrected[XTable$het==0]=0
XTable$het.corrected=Xhet.corrected

# get adjustment factors for heterozygosity (from autosomes model)
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"),verbose=TRUE)



# add adjusted het to otherInfo
otherInfo = left_join(otherInfo,select(XTable,het.corrected,logit.miss,het,miss,IID),by=c("PIID"="IID") )

# some derived variables
otherInfo$XY = otherInfo$Y.intensity/otherInfo$X.intensity
otherInfo$hetAutosome = Table$het[match(otherInfo$PIID,Table$IID)]
otherInfo$het.correctedAutosome = Table$het.corrected[match(otherInfo$PIID,Table$IID)]
otherInfo$missAutosome = Table$miss[match(otherInfo$PIID,Table$IID)]
otherInfo$logit.missAutosome = Table$logit.miss[match(otherInfo$PIID,Table$IID)]
otherInfo$XAuto.het = otherInfo$het/otherInfo$hetAutosome
otherInfo$XAuto.het.corrected = otherInfo$het.corrected/otherInfo$het.correctedAutosome
otherInfo$XAuto.miss = otherInfo$miss/otherInfo$missAutosome
otherInfo$XAuto.logit.miss = otherInfo$logit.miss/otherInfo$logit.missAutosome
otherInfo$XAuto.diff.logit.miss = otherInfo$logit.miss - otherInfo$logit.missAutosome

otherInfo$XPARmiss = PARmiss$F_MISS[match(otherInfo$PIID,PARmiss$IID)]
otherInfo$XPARlogit.miss = PARmiss$logit.miss[match(otherInfo$PIID,PARmiss$IID)]
otherInfo$XPARAuto.diff.logit.miss = otherInfo$XPARlogit.miss - otherInfo$logit.missAutosome

otherInfo$Ymiss = Ymiss$F_MISS[match(otherInfo$PIID,Ymiss$IID)]
otherInfo$Ylogit.miss = Ymiss$logit.miss[match(otherInfo$PIID,Ymiss$IID)]
otherInfo$YAuto.diff.logit.miss = otherInfo$Ylogit.miss - otherInfo$logit.missAutosome


# normalise x and y intensities for non-polymorphic probes
medianX = median(otherInfo$X.intensity[otherInfo$Inferred.Gender=="F"])
medianY = median(otherInfo$X.intensity[otherInfo$Inferred.Gender=="M"])
otherInfo$l2rNPX = log2(otherInfo$X.intensity) - log2(medianX)
otherInfo$l2rNPY = log2(otherInfo$Y.intensity) - log2(medianY)

# apply het correction to the ratio
Table$hetRatio = XTable$het[match(Table$IID,XTable$IID)]/Table$het
Table$hetRatio[Table$hetRatio==0]=NA
frmla = formula("hetRatio ~ (PC1 + PC2 + PC3 + PC4 + PC5 + PC6)^2 + I(PC1*PC1) + I(PC2*PC2) + I(PC3*PC3) + I(PC4*PC4) + I(PC5*PC5) + I(PC6*PC6)")
lmfit = lm(formula = frmla, data = Table)
betahat = coefficients(lmfit)
Intercept = betahat[1]
Table$hetRatio.corrected = Table$hetRatio - (predict(lmfit, newdata = Table) - Intercept)

otherInfo$hetRatio.corrected = Table$hetRatio.corrected[match(otherInfo$PIID,Table$IID)]
otherInfo$hetRatio.corrected[is.na(otherInfo$hetRatio.corrected)] = 0

# male/females
males = otherInfo$Inferred.Gender=="M"
females = otherInfo$Inferred.Gender=="F"

# mismatches
mismatches = otherInfo$Inferred.Gender!=otherInfo$Submitted.Gender
mismatchesM = mismatches & males
mismatchesF = mismatches & females

# duplicates
dupes = otherInfo$PIID%in%dup

print( paste0(sum(mismatches)," gender mismatches." ))
print( paste0(sum(mismatches&males)," male inferred and female submitted." ))
print( paste0(sum(mismatches&females)," female inferred and male submitted."))
print( paste0(sum(mismatches&!dupes)," gender mismatches not otherwise duplicates." ))


# get data from cnv summaries
load(paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries.RData"),verbose=TRUE)

# order by otherInfo samples
bafsum = outData[["baf"]][otherInfo$PIID,]
l2rsum = outData[["log2ratio"]][otherInfo$PIID,]

sum( rownames(bafsum)!=otherInfo$PIID ) 
sum( rownames(l2rsum)!=otherInfo$PIID )

# get list of SNPs used in cnv summaries (these files were created using ../SelectSNPs/create-plink-subsets.sh)
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

#### PLOTS on HET and Y:X intensity combinations
png( paste0("plots/",outFile,"-mismatches-Legend.png") ,width=1000,height=1000,res=150,bg="transparent")
plot.new()
legend("topleft",col=colDefinition,legend=names(colDefinition),bty="n",pch=4)
dev.off()

# XY ratio by raw het + corrected het
x = otherInfo$XY
y = otherInfo$het
png( paste0("plots/",outFile,"-YXintensityByHet-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=4,xlab="Y:X intensity",ylab="Heterozygosity in non-PAR X")
dev.off()

y = otherInfo$het.corrected
y[is.na(y)] = otherInfo$het[is.na(y)]
png( paste0("plots/",outFile,"-YXintensityByHetCorrected-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=4,xlab="Y:X intensity",ylab="PC-corrected heterozygosity in non-PAR X")
dev.off()

x = otherInfo$XY
y = otherInfo$het
png( paste0("plots/",outFile,"-YXintensityByHet-ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=shapes[OrderE],xlab="Y:X intensity",ylab="Heterozygosity in non-PAR X")
dev.off()

y = otherInfo$het.corrected
y[is.na(y)] = otherInfo$het[is.na(y)]
png( paste0("plots/",outFile,"-YXintensityByHetCorrected-ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=shapes[OrderE],xlab="Y:X intensity",ylab="PC-corrected heterozygosity in non-PAR X")
dev.off()


# plot ratio of X heterozygosity to autosomal heterozygosity

x = otherInfo$hetAutosome
y = otherInfo$het
Order = order.by.number.occurrences(colors)
png( paste0("plots/",outFile,"-XhetByAutoHet-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="Autosome het",ylab="non-PAR X heterozygosity")
dev.off()

y = otherInfo$XAuto.het
x = otherInfo$X.intensity
png( paste0("plots/",outFile,"-XHetAutoHetRatioByXintensity-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="X intensity",ylab="non-PAR X het / autosome het")
dev.off()


# plot X:Autosome heterozygosity by missing rate on X
y = otherInfo$XAuto.het
x = otherInfo$logit.miss
png( paste0("plots/",outFile,"-XHetAutoHetRatioByXmissing-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="non-PAR X missing rate",ylab="non-PAR X het / autosome het")
dev.off()


# plot X:Autosome heterozygosity by XY intensity ratio
y = otherInfo$XAuto.het
x = otherInfo$XY
png( paste0("plots/",outFile,"-XHetAutoHetRatioByYXintensity-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="Y:X intensity",ylab="non-PAR X het / autosome het")
dev.off()

x = otherInfo$XAuto.diff.logit.miss
y = otherInfo$XY
png( paste0("plots/",outFile,"-YXintensityByXAutoMissingRatio-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="Y:X intensity")
dev.off()



# plot XhetRate(Affy):Autosome heterozygosity by missing
y = otherInfo$Het.Rate/otherInfo$hetAutosome
x = otherInfo$logit.miss
png( paste0("plots/",outFile,"-XHet.RateAutoHetRatioByXmissing-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="non-PAR X logit missing",ylab="non-PAR X Het.Rate / autosome het")
dev.off()



# plot X:Autosome heterozygosity by X:Autosome missing
y = otherInfo$XAuto.het.corrected
x = otherInfo$XAuto.diff.logit.miss
png( paste0("plots/",outFile,"-XHetCorrectedAutoHetRatioByXAutoMissingRatio-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="non-PAR X / Autosome - PC-corrected heterozygosity")
abline(v=0,h=1,col="red",lty=3)
dev.off()

png( paste0("plots/",outFile,"-XHetCorrectedAutoHetRatioByXAutoMissingRatio-Ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=shapes[OrderE],xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="non-PAR X / Autosome - PC-corrected heterozygosity")
abline(v=0,h=1,col="red",lty=3)
dev.off()

# uncorrected het ratio x diff missing
y = otherInfo$XAuto.het
x = otherInfo$XAuto.diff.logit.miss
png( paste0("plots/",outFile,"-XHetAutoHetRatioByXAutoMissingRatio-Ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=shapes[OrderE],xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="non-PAR X / Autosome - heterozygosity")
abline(v=0,h=1,col="red",lty=3)
dev.off()

# het ratio - post-ratio correction
y = otherInfo$hetRatio.corrected
png( paste0("plots/",outFile,"-XAutoHetRatioCorrectedByXAutoMissingRatio-Ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=shapes[OrderE],xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="PC-corrected: (non-PAR X / Autosome heterozygosity)")
abline(v=0,h=1,col="red",lty=3)
dev.off()


# X intensity by Y intensity
x = otherInfo$Y.intensity
y = otherInfo$X.intensity
png( paste0("plots/",outFile,"-XintensityByYintensity-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="Y intensity (NP probes)",ylab="X intensity (NP probes)")
dev.off()

# X x Y l2r non-polymorphic probes (is this more resolution than cnv data?)
# some odd samples here
oddIntensitySamples = c(otherInfo$PIID[( otherInfo$l2rNPY < -4) & ( otherInfo$l2rNPX < -3) & (otherInfo$Inferred.Gender=="F")],otherInfo$PIID[( otherInfo$l2rNPY < -3) & ( otherInfo$l2rNPX < -3.5) & (otherInfo$Inferred.Gender=="M")])

x = otherInfo$l2rNPY
y = otherInfo$l2rNPX
png( paste0("plots/",outFile,"-l2rXintensityByl2rYintensity-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="log2ratio Y intensity (NP probes)",ylab="log2ratio X intensity (NP probes)")

points(x[otherInfo$PIID%in%oddIntensitySamples],y[otherInfo$PIID%in%oddIntensitySamples],pch=shapes[otherInfo$PIID%in%oddIntensitySamples],col="red")

dev.off()



# X l2r non-polymorphic probes by cnv data
y = rowMeans(cbind(otherInfo$X.intensity,otherInfo$Y.intensity))
x = l2rsum[,"means.X"]
png( paste0("plots/",outFile,"-meanXYintensityByL2RmeansX-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="non-PAR X mean log2ratio",ylab="mean X andd Y intensity (NP probes)")
dev.off()


# average intensity across both X and Y
y = otherInfo$X.intensity
x = rowMeans(cbind(otherInfo$l2rNPX,otherInfo$l2rNPY))
png( paste0("plots/",outFile,"-l2rXintensityByMeanXandYl2rintensity-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="mean YandX l2r intensity (NP probes)",ylab="log2ratio X intensity (NP probes)")
dev.off()


######## PLOT CNV DATA


# mean nonPAR-X x mean Y log2ratio (this should be similar to the X.intensity Y.intensity values)
x = l2rsum[,"means.X"]
y = l2rsum[,"means.Y"]
png( paste0("plots/",outFile,"-L2Rmean-nonPARXByY-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="mean log2ratio non-PAR X",ylab="mean log2ratio Y")
dev.off()

x = otherInfo$XY
y = l2rsum[,"means.X"]
png( paste0("plots/",outFile,"-L2Rmean-YXintensityBynonPARXbaf-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab="Y:X intensity (non-polymorphic SNPs)",ylab="mean log2ratio non-PAR X")
dev.off()

# plot combinations of chromosomes - l2r mean - quite informative
for( i in c(names(sexChroms2))){
    print(i)
    for( j in names(sexChroms2)[(which(names(sexChroms2)==i)+1):length(sexChroms2)]){
        png( paste0("plots/",outFile,"-L2Rmean-",i,"By",j,"-mismatches.png") ,width=1000,height=1000,res=150)
        x = l2rsum[,paste0("means.",i)]
        y = l2rsum[,paste0("means.",j)]
        plot(x[Order],y[Order],col=colors[Order],pch=4,xlab=paste0(i," mean log2ratio"),ylab=paste0(j," mean log2ratio"))
        dev.off()
    }
}

# plot baf missing rates  (not that informative...)
for( i in c(names(sexChroms2))){
    print(i)
    for( j in names(sexChroms2)[(which(names(sexChroms2)==i)+1):length(sexChroms2)]){
        png( paste0("plots/",outFile,"-bafMiss-",i,"By",j,"-mismatches.png") ,width=1000,height=1000,res=150)
        x = bafsum[,paste0("miss.",i)]/nsnps[i]
        y = bafsum[,paste0("miss.",j)]/nsnps[j]
        plot(x[Order],y[Order],col=colors[Order],pch=4,xlab=paste0(i," fraction missing"),ylab=paste0(j," fraction missing"))
        dev.off()
    }
}

# x call missing rates / autosome missing rates
x = otherInfo$XAuto.diff.logit.miss
for( i in c(names(sexChroms2))){
    png( paste0("plots/",outFile,"-L2Rmean-",i,"ByXAutoMissingRatio-mismatches.png") ,width=1000,height=1000,res=150)
    y = l2rsum[,paste0("means.",i)]
    plot(x[Order],y[Order],col=colors[Order],pch=4,xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab=paste0(i," mean log2ratio"))
    dev.off()
}


# standard deviations of baf
for( i in c(names(sexChroms2))){
    print(i)
    for( j in names(sexChroms2)[(which(names(sexChroms2)==i)+1):length(sexChroms2)]){
        png( paste0("plots/",outFile,"-BAFsds-",i,"By",j,"-mismatches.png") ,width=1000,height=1000,res=150)
        x = bafsum[,paste0("sds.",i)]
        y = bafsum[,paste0("sds.",j)]
        plot(x[Order],y[Order],col=colors[Order],pch=4,xlab=paste0(i," standard deviation BAF"),ylab=paste0(j," standard deviation BAF"))
        dev.off()
    }
}

# mean l2r x pico-green concentration (this might make a difference!!)
png( paste0("plots/",outFile,"-L2Rmean-YbyPicoGreen-mismatches.png") ,width=1000,height=1000,res=150)
x = l2rsum[,"means.Y"]
y = otherInfo$Internal.Pico..ng.uL.
plot(x[males],y[males],col=colors[males],pch=4,ylab="PicoGreen concentration",xlab=paste0("Y mean log2ratio"))
dev.off()
png( paste0("plots/",outFile,"-L2Rmean-adjustedYbyPicoGreen-mismatches.png") ,width=1000,height=1000,res=150)
l = lm(x[males] ~ y[males])
bh = coefficients(l)
Intercept = bh[1]
x2 = x[males] - (predict(l, newxreg = x[males]) - Intercept)
plot(x2,y[males],col=colors[males],pch=4,ylab="PicoGreen concentration",xlab=paste0("Y mean after adjusting for picogreen"))
dev.off()


# histograms of mean l2r and baf
png( paste0("plots/",outFile,"-L2Rmean-hist.png") ,width=1000,height=2000,res=150)
par(mfrow=c(5,1))
for(chr in names(sexChroms2)){
    x = l2rsum[,paste0("means.",chr)]
    hist(x,breaks=1000,xlab=paste0(chr," mean log2ratio"),main="")
}
dev.off()

png( paste0("plots/",outFile,"-BAFmean-hist.png") ,width=1000,height=2000,res=150)
par(mfrow=c(5,1))
for(chr in names(sexChroms2)){
    x = bafsum[,paste0("means.",chr)]
    hist(x,breaks=1000,xlab=paste0(chr," mean baf"),main="")
}
dev.off()


# histograms of sd l2r and baf
png( paste0("plots/",outFile,"-L2Rsd-hist.png") ,width=1000,height=2000,res=150)
par(mfrow=c(5,1))
for(chr in names(sexChroms2)){
    x = l2rsum[,paste0("sds.",chr)]
    hist(x,breaks=1000,xlab=paste0(chr," standard deviation log2ratio"),main="")
}
dev.off()

png( paste0("plots/",outFile,"-BAFsd-hist.png") ,width=1000,height=2000,res=150)
par(mfrow=c(5,1))
for(chr in names(sexChroms2)){
    x = bafsum[,paste0("sds.",chr)]
    hist(x,breaks=1000,xlab=paste0(chr," standard deviation baf"),main="")
}
dev.off()


# mean l2r x mean baf
for(chr in names(sexChroms2)){
    png( paste0("plots/",outFile,"-BAFmeanByL2Rmean-",chr,"mismatches.png") ,width=1000,height=2000,res=150)
    x = bafsum[,paste0("means.",chr)]
    y = l2rsum[,paste0("means.",chr)]
    plot(x[Order],y[Order],col=colors[Order],pch=shapes[Order],xlab=paste0(chr," mean baf"),ylab=paste0(chr," mean log2ratio"))
    dev.off()
}

# PCA on all mean values?? (didn't really work)
dataSum= cbind(l2rsum[,grep("means",colnames(bafsum))],bafsum[,grep("miss",colnames(bafsum))])
dataSum=dataSum[,!grepl("MT",colnames(dataSum))]

pca = princomp(dataSum)

png( paste0("plots/",outFile,"-pca-%02d-mismatches.png") ,width=1000,height=1000,res=150)
for(i in 1:ncol(dataSum)){
    x = pca$scores[,i]
    y = pca$scores[,(i+1)]
    plot(x[Order],y[Order],col=colors[Order],pch=1,xlab=paste0("PC ",i),ylab=paste0("PC ",(i+1)))
}
dev.off()


# plot age vs j chromosome intensity
x = otherInfo$Age.when.attended.assessment.centre
for( j in c(names(sexChroms2))){
    print(j)
    png( paste0("plots/",outFile,"-L2Rmean-",j,"ByAge-mismatches.png") ,width=1000,height=1000,res=150)
    y = l2rsum[,paste0("means.",j)]
    plot(x[Order],y[Order],col=colors[Order],pch=4,xlab="Age at assessment centre",ylab=paste0(j," mean log2ratio"))
    dev.off()
}

# just plot inferred females with age on X chroms
colors2 = colors; colors2[!females] = "transparent"
for( j in names(sexChroms2)[c(1,3,4)]){    
    png( paste0("plots/",outFile,"-L2Rmean-",j,"ByAge-mismatches.png") ,width=1000,height=1000,res=150)
    y = l2rsum[,paste0("means.",j)]
    plot(x[Order],y[Order],col=colors2[Order],pch=4,xlab="Age at assessment centre",ylab=paste0(j," mean log2ratio"))
    dev.off()
}


# LM for age/PAR1
y = l2rsum[,paste0("means.XPAR1")]
y = l2rsum[,paste0("means.X")]
y = l2rsum[,paste0("means.Y")]
x = otherInfo$Age.when.attended.assessment.centre
tF = lm(y[females]~x[females])
tM = lm(y[males]~x[males])


# 

# What if there are different SNPs with different loss-prevalences (e.g on the Y)?
# Look at loss as a function of age - are there regions that have higher correlations?




############ SELECT SET OF SAMPLES TO EXCLUDE FROM PHASING

X = l2rsum[,"means.X"]
Y = l2rsum[,"means.Y"]
Ynorm = Y - X

## criteria1:
criteria=1
# based only on mean X and Y cnv values.
odd1 = females & ( X < -0.17 )  # genotyped as XX but low l2r
odd2 = females & ( X > 0.145 ) # genotyped as XX but high l2r
odd3 = ( Y >= -1 ) & ( X > -0.2 ) # potentially XXY
odd4 = ( Y >= 0.23 ) # potentially XYY or XXYY
odd = odd1 | odd2 | odd3 | odd4

## criteria2:
criteria=2
# add missing rate criteria
# also normalise males for Y chrom l2r values
missThreshold = 0.05  #  (this threshold matched the autosome one. see compile-het-missing-exclusion-list.R)
missx = otherInfo$miss > missThreshold
# get list of samples with odd properties for PAR1
misspar = otherInfo$XPARmiss > missThreshold
# 401 samples with high missing rates on PAR
# 85 of these are also excluded on X from criteria 1
lowYthreshold = mean(Ynorm[males]) + 3*sd(Ynorm[males])*c(1,-1)
lowY = ( Ynorm <= lowYthreshold[2]) & males

odd4 = ( Ynorm >= lowYthreshold[1] ) # potentially XYY or XXYY
odd1 = females & ( X < -0.17 )  # genotyped as XX but low l2r
odd2 = females & ( X > 0.145 ) # genotyped as XX but high l2r
odd3 = ( Y >= -1 ) & ( X > -0.2 ) # potentially XXY

odd = odd1 | odd2 | odd3 | odd4 | missx  # this is for the X
oddPar = odd1 | odd2 | odd3 | odd4 | lowY | missx | misspar # this is for the PAR

## criteria3:
criteria=3
# as for 2, but don't apply the lowY threshold
odd = odd1 | odd2 | odd3 | odd4 | missx  # this is for the X
oddPar = odd1 | odd2 | odd3 | odd4 | missx | misspar # this is for the PAR



######## plot these

shapesEx = rep(1,length(X))
shapesEx[!odd] = 4
shapesEx = rep(1,length(X))
shapesEx[!odd] = 4
colorsEx = colors
colorsEx[odd] = "gray"
colorsEx[(!odd)&oddPar] = "pink"

x = X; y = Y
png( paste0("plots/",outFile,"-L2Rmean-nonPARXByY-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="mean log2Ratio non-PAR X",ylab="mean log2Ratio Y")
abline(v = c(-0.17,0.15,-0.2),lty=3)
abline(h = c(-1,0.23),lty=3)
dev.off()


# Y chromosome loss - normalised by X
x = X; y = Ynorm
png( paste0("plots/",outFile,"-L2Rmean-nonPARXByYnorm-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="mean log2Ratio non-PAR X",ylab="mean log2Ratio Y/X")
abline(v = c(-0.17,0.15,-0.2),lty=3)
abline(h = c(-1,0.23),lty=3)
dev.off()


# with missing rates
x = otherInfo$XAuto.diff.logit.miss
y = X
png( paste0("plots/",outFile,"-L2Rmean-nonPARXByXAutoMissingRatio-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="mean log2Ratio X")
dev.off()

x = otherInfo$logit.miss
y = X
png( paste0("plots/",outFile,"-L2Rmean-nonPARXByXMissing-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="non-PAR X logit missing",ylab="mean log2Ratio non-PAR X")
abline(v = logit(0.05),lty=3)
dev.off()

y = Y
png( paste0("plots/",outFile,"-L2Rmean-YByXAutoMissingRatio-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="( non-PAR X logit missing ) - ( Autosome logit missing )",ylab="mean log2Ratio Y")
dev.off()
x = otherInfo$XPARAuto.diff.logit.miss
png( paste0("plots/",outFile,"-L2Rmean-YByXPARAutoMissingRatio-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="( PAR X logit missing ) - ( Autosome logit missing )",ylab="mean log2Ratio Y")
dev.off()
x = otherInfo$XPARlogit.miss
png( paste0("plots/",outFile,"-L2Rmean-YByXPARMissing-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="PAR X logit missing",ylab="mean log2Ratio Y")
abline(v = logit(0.05),lty=3)
dev.off()

# Y chromosome missing rates
x = otherInfo$Ylogit.miss
png( paste0("plots/",outFile,"-L2Rmean-YByYMissing-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="Y logit missing",ylab="mean log2Ratio Y")
dev.off()
png( paste0("plots/",outFile,"-L2Rmean-YByYMissing-Ethnicity.png") ,width=1000,height=1000,res=150)
plot(x[OrderE],y[OrderE],col=colorsE[OrderE],pch=4,xlab="Y logit missing",ylab="mean log2Ratio Y")
dev.off()


# XPAR1 l2r with Y chrom rates
x = Ynorm
y = l2rsum[,"means.XPAR1"]
png( paste0("plots/",outFile,"-L2Rmean-YnormByXPAR1-exclusionsCriteria",criteria,"-mismatches.png") ,width=1000,height=1000,res=150)
plot(x[Order],y[Order],col=colorsEx[Order],pch=4,xlab="mean log2Ratio Y/X",ylab="mean log2Ratio XPAR1")
dev.off()


# get list of samples with odd karyotype (for non-PAR X)
x0 = rownames(l2rsum)[odd1]
xxx = rownames(l2rsum)[odd2]
xxy = rownames(l2rsum)[odd3]
yyy = rownames(l2rsum)[odd4]
lowy = names(Ynorm)[lowY]
hiMissX = otherInfo$PIID[missx] # only four extra samples are excluded here
hiMissPAR = otherInfo$PIID[misspar] # only four extra samples are excluded here
allPhaseExclusionsX = rownames(l2rsum)[odd]
allPhaseExclusionsPAR = rownames(l2rsum)[oddPar] 
    
if(criteria == 1) {

    save(x0,xxx,xxy,yyy,allPhaseExclusionsX,file = paste0(outFile,"-phaseExclusionsCriteria",criteria,".Rdata"))
    write.table(allPhaseExclusionsX, file = paste0(outFile,"-phaseExclusionsCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
                                        # 667 samples excluded with Criteria1
}

if(criteria == 2) {
    
    save(x0,xxx,xxy,yyy,lowy,hiMissX,hiMissPAR,allPhaseExclusionsX,allPhaseExclusionsPAR,file = paste0(outFile,"-phaseExclusionsCriteria",criteria,".Rdata"))
    
                                        # write list of lowY samples
    write.table(lowy,file = paste0(outFile,"-lowYchromosomeMalesCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(hiMissPAR,file = paste0(outFile,"-hiMissPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

    write.table(allPhaseExclusionsX, file = paste0(outFile,"-phaseExclusionsXnonPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(allPhaseExclusionsPAR, file = paste0(outFile,"-phaseExclusionsPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

    
}

if(criteria == 3) {
    # difference is no lowy exclusion here
    save(x0,xxx,xxy,yyy,hiMissX,hiMissPAR,allPhaseExclusionsX,allPhaseExclusionsPAR,file = paste0(outFile,"-phaseExclusionsCriteria",criteria,".Rdata"))
    
    write.table(hiMissPAR,file = paste0(outFile,"-hiMissPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

    write.table(allPhaseExclusionsX, file = paste0(outFile,"-phaseExclusionsXnonPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(allPhaseExclusionsPAR, file = paste0(outFile,"-phaseExclusionsPARCriteria",criteria,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

    
}


# NOTE: for 2 and 3, XnonPAR exclusion lists are the same. 673 sample (667 for criteria 1)
# However, the PAR exclusions differ. For 2: 4221 samples in total; 988 for criteria 3.











############ SELECT SET OF SAMPLES TO CALL HMM ON (NOT THE SAME SET FOR PHASING EXCLUSIONS!!)
# for each chromsome pick individuals outside of normal range on mean log2ratio

## X by Y log2ratio
# ellipse around females
x = l2rsum[,"means.X"]; y = l2rsum[,"means.Y"]
w=0.2 # the radius of the ellipse
isInEllipse = is.in.ellipse(l2rsum[,"means.X"], l2rsum[,"means.Y"] ,x0=mean(x[y < -1]),vx= (w/2)^2,y0=mean(y[y < -1]),vy = w^2,cov=0)
keep1 = (!isInEllipse) & (y < -1)
keep2 = (y > -0.1) & (x > -0.25)
keep3 = y > 0.18

## XPAR1 x XPAR2
x = l2rsum[,"means.XPAR1"]; y = l2rsum[,"means.XPAR2"]
keep4 = (x > -0.3) & (y < -0.3)
keep5 = (y > 0.3) | (x < -0.3 )
keep51 = x > 0.18

## XPAR2 & X
x = l2rsum[,"means.X"]; y = l2rsum[,"means.XPAR2"]
keep6 = ( y < -0.3 ) & ( x > -0.2 )

keepAll = (keep1 + keep2 + keep3 + keep4 + keep5 + keep51 + keep6) > 0 

## plot outlier samples
colorsKeep = rep("gray",length(keepAll)); colorsKeep[keepAll] = colors[keepAll]
orderKeep = order.by.number.occurrences(colorsKeep)

for( i in c(names(sexChroms2))){
    print(i)
    for( j in names(sexChroms2)[(which(names(sexChroms2)==i)+1):length(sexChroms2)]){
        extra = NULL
        if(( "X" %in% c(i,j) ) & ("Y" %in% c(i,j)) ) extra = (keep1 + keep2 + keep3) > 0
        if(( "XPAR1" %in% c(i,j) ) & ("XPAR2" %in% c(i,j)) ) extra = (keep4 + keep5 + keep51) > 0
        if(( "X" %in% c(i,j) ) & ("XPAR2" %in% c(i,j)) ) extra = keep6
        x = l2rsum[,paste0("means.",i)]
        y = l2rsum[,paste0("means.",j)]

        png( paste0("plots/",outFile,"-L2Rmean-",i,"By",j,"-oddOnly-mismatches.png") ,width=1000,height=1000,res=150)
        plot(x[keepAll],y[keepAll],col=colors[keepAll],pch=4,xlab=paste0(i," mean log2ratio"),ylab=paste0(j," mean log2ratio"),xlim=c(min(x),max(x)),ylim=c(min(y),max(y)))
        if(!is.null(extra)) points(x[extra],y[extra],pch=4,col="red")
        dev.off()
        
        png( paste0("plots/",outFile,"-L2Rmean-",i,"By",j,"-oddColored-mismatches.png") ,width=1000,height=1000,res=150)        
        plot(x[orderKeep],y[orderKeep],col=colorsKeep[orderKeep],pch=4,xlab=paste0(i," mean log2ratio"),ylab=paste0(j," mean log2ratio"))        
        dev.off()

    }
}


## ~1,300 samples in keepAll.
write.table(otherInfo$PIID[keepAll],file=paste0(outFile,"-oddSamples.txt"),quote=FALSE, col.names=FALSE,row.names=FALSE)

## Add some more 'normal' samples as controls in the HMM
set.seed(1234567)
ms = sample(otherInfo$PIID[(!keepAll) & (males)],500,replace=FALSE)
fs = sample(otherInfo$PIID[(!keepAll) & (females)],500,replace=FALSE)
write.table(c(ms,fs),file=paste0(outFile,"-controlSamples.txt"),quote=FALSE, col.names=FALSE,row.names=FALSE)


## create list of classifications based on mean values

x1init = rownames(l2rsum)[ (l2rsum[,"means.X"] < -0.2)]  # one X
x2init = rownames(l2rsum)[ (l2rsum[,"means.X"] > -0.2) & (meanL2r[,"means.X"] < 0.12)]  # two X
x3init = rownames(l2rsum)[l2rsum[,"means.X"] >= 0.12] # three X
y2init = rownames(l2rsum)[ (l2rsum[,"means.Y"] >= -1) & (meanL2r[,"means.Y"] < 0.2)]  # one Y
y3init = rownames(l2rsum)[l2rsum[,"means.Y"] >= 0.2] # two Y

par1init = rownames(l2rsum)[(l2rsum[,"means.XPAR1"] < -0.4)] # one PAR1
par3init = rownames(l2rsum)[(l2rsum[,"means.XPAR1"] >= 0.1)] # three PAR1s
par1init2 = rownames(l2rsum)[(l2rsum[,"means.XPAR2"] < -0.3)] # one PAR2
par3init2 = rownames(l2rsum)[(l2rsum[,"means.XPAR2"] >= 0.2)] # three PAR2s

save(x1init,x2init,x3init,y2init,y3init,par1init,par3init,par1init2,par3init2,file=paste0(outFile,"-initialClassificationsOddSamples.Rdata"))
