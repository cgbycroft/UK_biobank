#############
# This script takes a list of phenotype fields, an optional file with sample IDs, and returns a nice text file for running SNP-test or BOLT, or any other test for that matter.
# It uses the extraction script created by W. Rayner (wrayner@well.ox.ac.uk), A. Mahajan (anubha@well.ox.ac.uk) 2015
#############

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

h=args
for(s in h){
    source(s)
}


setwd('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/YXintensity')

#############
vers="v1"  # Y:X ratio insteat of Y intensity!! Age is important here, first just leave it out. no obvious binary trait.
#############
vers="v2"  # as v1, but include mean l2ratio and ratio of Y:X mean l2r as phenotype.
#############
vers="v3"  # as v2 but include females!  just X chrom intensity analysis. Also include mean non-PAR and PAR X l2r
#############


# l2r only actually used for v2 and v3
# get l2r info (extracted in compute-l2r-summaries.R or something)
load(paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries.RData"),verbose=TRUE)
meanL2rY = outData[["log2ratio"]][,"means.Y"]
meanL2rX = outData[["log2ratio"]][,"means.X"]


## if vers == v2 skip between %%%%%%%%%
#%%%%%%%%%%
genoInfo = read.multiple.batch.info(fields = batchInfoFields)

# get Array
genoInfo$Array = "UKBB"
genoInfo$Array[genoInfo$Batch%in%ukbileve.batches()] = "UKBL"

# get Y:X ratio
genoInfo$Y.X.intensity.ratio = genoInfo$Y.intensity/genoInfo$X.intensity

# read in data for het/missing
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"),verbose=T)

# read in data for PCs (%%%%%% CHANGE THIS TO FULL SET LATER %%%%%%%%)
load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-round2.RData"),verbose=T)
# USED THIS FOR V1: load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData"),verbose=T)

# fold in the PCs
genoInfo = tbl_df(genoInfo)
genoInfo = left_join(genoInfo,PCs,by=c("PIID"="PIID"))

genoInfo$X.mean.l2r = meanL2rX[genoInfo$PIID]
genoInfo$Y.mean.l2r = meanL2rY[genoInfo$PIID]
genoInfo$YXratio.mean.l2r = genoInfo$Y.mean.l2r - genoInfo$X.mean.l2r
genoInfo$XYratio.mean.l2r = genoInfo$X.mean.l2r - genoInfo$Y.mean.l2r
genoInfo$XYdiff.mean.l2r = log2(2^(genoInfo$X.mean.l2r) - 2^(genoInfo$Y.mean.l2r))

sum(is.na(genoInfo$Y.mean.l2r))


###### get list of samples to use in GWAS
# Keep just males, (start with White British set, too) and exclude any related samples

# get list of QC'd samples
## NOTE: this is preliminary as of May 14th. I.e haven't run aberrant.
hetMiss = Table
highHet = hetMiss$IID[hetMiss$het.corrected > 0.2]
highMiss = hetMiss$IID[hetMiss$logit.miss > -3]
sexMismatch = genoInfo$PIID[genoInfo$Submitted.Gender!=genoInfo$Inferred.Gender]
other = c()

if(vers%in%c("v1")) other = genoInfo$PIID[genoInfo$Inferred.Gender=="F"]
if(vers%in%c("v3")) other = genoInfo$PIID[ (genoInfo$Inferred.Gender=="M") | ( genoInfo$Y.mean.l2r > -1) ] # i.e exclude males or females with XXX

excludeSamples = unique(c(highHet,highMiss,sexMismatch,other))
write.table(excludeSamples,file="excludeSamples.tmp",quote=FALSE,row.names=FALSE,col.names=FALSE)

# get unrelated samples within the QC'd set (this script also finds individuals with many relatives)
args <- c(h,"-in",paste0(baseSampleQCDir,"/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0"), "-ibs0","0.0012", "-exclude", "excludeSamples.tmp","-outdir",".","-outfile",paste0("b1__b11-b001__b095-pair_batches.filtered-",vers) )

source( paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/postking/get-unrelated.R")) 


# subset genoInfo data
goodSamples = read.table(paste0("b1__b11-b001__b095-pair_batches.filtered-",vers,"-unrelated.txt"),header=FALSE,stringsAsFactors=FALSE)[,1]

genoInfoGood = filter(genoInfo,genoInfo$PIID%in%goodSamples)


# Now write phenotype file for plink (or BOLT-LMM). E.g continuous vs categorical; missing value code etc.
BOLTCols = c("PIID","PIID",grep("PC",colnames(genoInfo),value=TRUE),batchInfoFields[c(1:8,10,11,17,19,20,23,24,28,29,32:52,54,56)],"Y.X.intensity.ratio","Array","X.mean.l2r","Y.mean.l2r","YXratio.mean.l2r","XYratio.mean.l2r","XYdiff.mean.l2r")

toPrint = genoInfoGood[,BOLTCols]
colnames(toPrint)[1:2] = c("FID","IID")

# Write any blank spaces as NA!!!
toPrint[toPrint==""] = NA

p1 = sapply(colnames(toPrint),function(x) sum(toPrint[[x]] < 0,na.rm=TRUE))
p2 = sapply(colnames(toPrint),function(x) sum(is.na(toPrint[[x]])))
p3 = sapply(colnames(toPrint),function(x) sum(grepl(" ",toPrint[[x]]),na.rm=TRUE))
p4 = sapply(colnames(toPrint),function(x) sum(toPrint[[x]]=="",na.rm=TRUE) )

# add in mosaics binary phenotypes
if( vers %in% c("v3") ) {

    # these values match those set in gender-checks.R
    toPrint$XXX = 0
    toPrint$XXX[toPrint$X.mean.l2r > 0.145] = 1
    toPrint$X0 = 0
    toPrint$X0[toPrint$X.mean.l2r < -0.17] = 1
    toPrint$X0[toPrint$XXX==1] = NA # exclude samples with XXX

}


#%%%%%%%%%%

if( vers %in% c("v2")){
    genoInfoGood = read.table(paste0(baseSampleQCDir,"/data/GWAS/YXintensity/YchromPhenotypesForBOLT-v1.txt"),header=TRUE)
    
    genoInfoGood$X.mean.l2r = meanL2rX[genoInfoGood$IID]
    genoInfoGood$Y.mean.l2r = meanL2rY[genoInfoGood$IID]

    genoInfoGood$YXratio.mean.l2r = genoInfoGood$Y.mean.l2r - genoInfoGood$X.mean.l2r

    # replace with better PCs
    load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-round2.RData"),verbose=T)
    genoInfoGood2 = genoInfoGood[,!grepl("PC",colnames(genoInfoGood))]
    genoInfoGood2 = left_join(genoInfoGood2,PCs,by=c("IID"="PIID"))

    genoInfoGood2 = genoInfoGood2[,!colnames(genoInfoGood2)%in%c( "Pops" ) ]
    
    # for Xchrom GWAS exclude males with high Ychrom intensity, and with XXY karyotype
    exc = (genoInfoGood$YXratio.mean.l2r > 0.7) | (genoInfoGood$Y.mean.l2r > 0.22) | ( (genoInfoGood$Y.mean.l2r > -0.3 ) & (genoInfoGood$YXratio.mean.l2r < 0) ) | ( (genoInfoGood$Y.mean.l2r > -0.15 ) & (genoInfoGood$YXratio.mean.l2r < 0.15) )  # 183 males excluded
    colors = rep("black",dim(genoInfoGood)[1]); colors[exc] = "red"

    y = genoInfoGood$YXratio.mean.l2r
    x = genoInfoGood$Y.mean.l2r
    png( paste0("plots/Y-X-l2r-by-YmeanL2r-",vers,"-beforeExclusions.png"),height=1000,width=1000, res=150 )
    plot(x,y,xlab="Y mean l2r",ylab="Y/X mean l2r",main=paste0(n," QC'd males"),col=colors)
    dev.off()
    
    toPrint = genoInfoGood2[!exc,]  # 205333 individuals
}

    
write.table(toPrint,file=paste0(baseSampleQCDir,"/data/GWAS/YXintensity/YchromPhenotypesForBOLT-",vers,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

# can we read it? i.e are the black spaces sorted out?
test = read.table(paste0(baseSampleQCDir,"/data/GWAS/YXintensity/YchromPhenotypesForBOLT-",vers,".txt"),header=TRUE)


###### SOME PLOTS
toPrint = test

x = toPrint$Age.when.attended.assessment.centre
y = toPrint$Y.X.intensity.ratio
n = sum(!is.na(x)&(!is.na(y)))

# phenotype by age
png( paste0("plots/Y-X-intensity-by-Age-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Age.when.attended.assessment.centre",ylab="Y/X average intensities",main=paste0(n," QC'd samples"))
dev.off()

y = toPrint$YXratio.mean.l2r
n = sum(!is.na(x)&(!is.na(y)))
png( paste0("plots/Y-X-l2r-by-Age-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Age.when.attended.assessment.centre",ylab="Y/X mean l2r",main=paste0(n," QC'd samples"))
dev.off()

y = toPrint$XYratio.mean.l2r # used for female version
n = sum(!is.na(x)&(!is.na(y)))
png( paste0("plots/X-Y-l2r-by-Age-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Age.when.attended.assessment.centre",ylab="X/Y mean l2r",main=paste0(n," QC'd samples"))
dev.off()

y = toPrint$XYdiff.mean.l2r # used for female version
n = sum(!is.na(x)&(!is.na(y)))
png( paste0("plots/X-Y-diff-l2r-by-Age-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Age.when.attended.assessment.centre",ylab="(X - Y) mean l2r",main=paste0(n," QC'd samples"))
dev.off()


x = toPrint$Y.X.intensity.ratio
y = toPrint$YXratio.mean.l2r
png( paste0("plots/Y-X-l2r-by-YX.intensity-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Y:X intensity (Affy)",ylab="Y/X mean l2r",main=paste0(n," QC'd samples"))
dev.off()

x = toPrint$Y.mean.l2r
png( paste0("plots/Y-X-l2r-by-YmeanL2r-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Y mean l2r",ylab="Y/X mean l2r",main=paste0(n," QC'd samples"))
dev.off()

x = toPrint$XYratio.mean.l2r
y = toPrint$X.mean.l2r
png( paste0("plots/X-Y-l2r-by-XmeanL2r-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="X/Y mean l2r",ylab="X mean l2r",main=paste0(n," QC'd samples"))
dev.off()

x = toPrint$XYdiff.mean.l2r
y = toPrint$X.mean.l2r
png( paste0("plots/X-Y-diff-l2r-by-XmeanL2r-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="(X - Y) mean l2r",ylab="X mean l2r",main=paste0(n," QC'd samples"))
dev.off()


x = toPrint$Y.X.intensity.ratio
y = toPrint$XYratio.mean.l2r
png( paste0("plots/X-Y-l2r-by-YX.intensity-",vers,".png"),height=1000,width=1000, res=150 )
plot(x,y,xlab="Y:X intensity (Affy)",ylab="X/Y mean l2r",main=paste0(n," QC'd samples"))
dev.off()




# histogram of phenotype
png( paste0("plots/Y-X-intensity-by-Age-",vers,".png"),height=1000,width=1000, res=150 )
hist(y,xlab="Y/X average intensities",main=paste0(n," QC'd males"))
dev.off()



####### BELOW NOT DONE YET #######


# Create files for snp.test. They must have the same set of samples as in the plink files (just put them in the same order for each chromosome, for safety!)
toPrint = genoInfo[,BOLTCols] # must include everyone. Add extra lines with NA for extra samples (they'll be excluded anyway)
fam = read.table("/well/donnelly/ukbiobank_project_8874/interim_genotype_data/binaryPed/ukb8874.chr22.calls.fam",header=FALSE)
extraSamples =fam[,1][!fam[,1]%in%toPrint[,1]]
extra = toPrint[1:3,]; extra[,1] = extraSamples
toPrint = rbind(toPrint,extra)

toPrint$missing = 0
colnames(toPrint)[1:2] = c("ID_1","ID_2")
toPrint = toPrint[,c(1,2,ncol(toPrint),3:(ncol(toPrint)-1))]
#extraHeader = c(0,0,0,"B",rep("C",15),"P","P","B","C","D") version 2 only
extraHeader = c(0,0,0,"B",rep("C",15),"P","P","B","C","D","P")
cbind(extraHeader,colnames(toPrint))

rownames(toPrint) = toPrint$ID_1
# make separate phenotype file for each chromosome, in order of .fam file
Order22 = as.character(fam[,1])

for (chr in 1:22){
    print(chr)
    f = read.table(paste0("/well/donnelly/ukbiobank_project_8874/interim_genotype_data/binaryPed/ukb8874.chr",chr,".calls.fam"),header=FALSE)
    Order = as.character(f[,1])
    toPrintThisChrom = toPrint[Order,]
    print(sum(toPrintThisChrom$ID_1!=Order))
    print(sum(toPrintThisChrom$ID_1!=Order22))
    
    File = paste0("data/YchromPhenotypesForSNPTest-",vers,"-chrom",chr,".sample")
                                        #File = paste0("data/YchromPhenotypesForSNPTest-",vers,".sample")
    write.table(t(colnames(toPrint)),file=File,col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(t(extraHeader),file=File,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(toPrintThisChrom,file=File,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

exclusions = fam[,1][!fam[,1]%in%genoInfoGood[,1]]

write.table(exclusions,file=paste0("data/YchromPhenotypesForSNPTest-",vers,".exclusions"),col.names=FALSE,row.names=FALSE,quote=FALSE)


# ====> now run BOLT, or snptest
