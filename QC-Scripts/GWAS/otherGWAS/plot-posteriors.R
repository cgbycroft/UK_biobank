#####################
# Script to plot posterior distributions and calculate some summaries.
#####################

#####################
# Preliminaries

args = commandArgs(trailingOnly=TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)
library(plyr)

# args = c("-mergeddata","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/Posteriors-merged-IMP_GENO_GIANT.version1.0.2.RData","-plotdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts","-prior","0.2")

sourceNames = c("IMP"="UK Biobank (imputed data)","IMPSUB"="UK Biobank (imputed data subset)","GENO"="UK Biobank (genotyped data)","GIANT"="GIANT (2014)")

print(args)

title = args[which(args=="-title")+1]
dataFile = str_split(args[which(args=="-mergeddata")+1],",")[[1]]

plotdir=args[which(args=="-plotdir")+1]
outDir = args[which(args=="-outdir")+1]
print(outDir)

#chroms = parse.range.string(args[which(args=="-chr")+1])

###############
# What is the prior for the bayse factors (if not the default)?

if("-prior"%in%args) prior = as.numeric(args[which(args=="-prior")+1]) else prior=NA

###############
# QC THRESHOLDS?

if("-qc" %in% args){

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

title = paste0(title,".maf",minmaf,".info",mininfo)



############################
##### Load in the data.

# Takes about 30secs to load in
load(dataFile,verbose=TRUE)

nChroms = length(MergedData)
MergedDataChromOrder = sapply(MergedData,function(x) x$CHR[1])


############################
##### Filter for QC

print("Applying QC...")

QCNone = getQCsubset(MergedData,minmaf=0,mininfo=0,maxmiss=1) # Actually no qc
MergedDataQCNone = QCNone[[1]]

#QC2 = getQCsubset(MergedData,minmaf=0,mininfo=0.3,maxmiss=1) # Some things actually get excluded if MAF == 0
#MergedDataQC2 = QC2[[1]]
#print( colSums(QC4[[2]]) )
#IMP   GENO  GIANT 
#663469      0    249

QC3 = getQCsubset(MergedData,minmaf=0.001,mininfo=0.3,maxmiss=1)
MergedDataQC3 = QC3[[1]]
#IMP    GENO   GIANT 
#3390055    8879     249

#QC4 = getQCsubset(MergedData,minmaf=0.001,mininfo=0.3,maxmiss=0.05)
#MergedDataQC4 = QC4[[1]]
#print( colSums(QC4[[2]]) )
#    IMP    GENO   GIANT 
#3390055   12099     249

QC5 = getQCsubset(MergedData,minmaf=0.001,mininfo=0.9,maxmiss=0.05)
MergedDataQC5 = QC5[[1]]
print( colSums(QC5[[2]]) )
#    IMP    GENO   GIANT 
#3721738   12099     249

qcNonetitle = paste0("noqc")
#qc2title = paste0(".maf",0,".info",0.3)
qc3title = paste0(".maf",0.001,".info",0.3)
#qc4title = paste0(".maf",0.001,".info",0.3,".maxmiss",0.05)
qc5title = paste0(".maf",0.001,".info",0.9,".maxmiss",0.05)

print(qcNonetitle)
print( colSums(QCNone[[2]]) )
#print(qc2title)
#print( colSums(QC2[[2]]) )
print(qc3title)
print( colSums(QC3[[2]]) )
#print(qc4title)
#print( colSums(QC4[[2]]) )
print(qc5title)
print( colSums(QC5[[2]]) )



############################
##### Pick a set of regions to compare

# Regions based on QC4

topRegionsGIANT = getTopHits(type="GIANT",MergedData,MergedDataQC3)
topRegionsIMP = getTopHits(type="IMP",MergedData,MergedDataQC3)
#topRegionsGENO = getTopHits(type="GENO",MergedData,MergedDataQC3)
#topRegionsGENO2 = getTopHits(type="GENO",MergedData,MergedDataQC3) # <=== QC 3 makes not difference


exclIMP = topRegionsIMP$region_ID[which(topRegionsIMP$IMP..P >= signifLevel)] # This excludes 56 regions
#exclGENO = topRegionsGENO$region_ID[which(topRegionsGENO$GENO..P >= signifLevel)] # This excludes an extra 37 regions
#exclGENO2 = topRegionsGENO2$region_ID[which(topRegionsGENO2$GENO..P >= signifLevel)] # This excludes an extra 37 regions
 
theseRegions = getNonOverlappingHits(topRegionsGIANT,AllRegions)
exclOverlapping = topRegionsGIANT$region_ID[!topRegionsGIANT$region_ID%in%theseRegions$region_ID]

theseRegionsGIANT_and_IMP = getNonOverlappingHits(topRegionsGIANT[!topRegionsGIANT$region_ID%in%exclIMP,],AllRegions)

#theseRegionsGIANT_and_GENO = getNonOverlappingHits(topRegionsGIANT[!topRegionsGIANT$region_ID%in%exclGENO,],AllRegions)

#theseRegionsGIANT_and_GENO_and_IMP = getNonOverlappingHits(topRegionsGIANT[!topRegionsGIANT$region_ID%in%c(exclGENO,exclIMP),],AllRegions)  # <== The same as GIANT_AND_GENO

#test = getTopHits(type="GIANT",MergedData,MergedDataQC4)

# HOW OFTEN IS THE TOP SNP in the Region not actually the one we originally picked?
chrs = theseRegionsGIANT_and_IMP$CHR[theseRegionsGIANT_and_IMP$SNP2 %in% SNPs]
index = theseRegionsGIANT_and_IMP$region.index[theseRegionsGIANT_and_IMP$SNP2 %in% SNPs]

a = sapply(1:length(chrs),function(i){
    AllRegions[[chrs[i]]][index[i],5]
    
})
CH = cbind(a,t)
k1 = sapply(CH[,1],function(s) {
    l = str_split(s,":")[[1]][2]
    sapply(l,str_split,"_")[[1]][1]
})
k2 = sapply(CH[,2],function(s) {
    l = str_split(s,":")[[1]][2]
    sapply(l,str_split,"_")[[1]][1]
})

print(  sum(k1!=k2) )  # <-======== once! "17:54694286_A_G"



##########################################
### Get scale factor for the Betas.

# Read in phenotype file
ph = read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt"),header=TRUE,stringsAsFactors=FALSE)

regressHeight = lm(Standing.height~Array+Age.when.attended.assessment.centre+Inferred.Gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=ph)

scaleFactor = sd(regressHeight$residuals)

################################
#scaleFactor = sd(regressHeight$residuals)#6.322236
#prior = 0.2
################################

columns = c("CHR","BP","BETA","SE","SNP","SNP2.swap","SNP2","P","MAF","swapped","INFO","F_MISS")

#theseRegions = theseRegions[1:5,]

#TOPLOT_NOQC = extractInfoForRegions(chrs=theseRegions$CHR,region.indices=theseRegions$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQCNone,scaleFactor=scaleFactor,prior=prior)

#TOPLOT_QC2 = extractInfoForRegions(chrs=theseRegions$CHR,region.indices=theseRegions$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC2,scaleFactor=scaleFactor,prior=prior)

#TOPLOT_QC3 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior)

#TOPLOT_QC4 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC4,scaleFactor=scaleFactor,prior=prior)

# only the intersections
#TOPLOT_NOQC_allStudiesOnly = extractInfoForRegions(chrs=theseRegions$CHR,region.indices=theseRegions$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQCNone,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE)

#TOPLOT_NOQC_GiantImpOverlap = extractInfoForRegions(chrs=theseRegions$CHR,region.indices=theseRegions$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQCNone,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","IMP"))

#TOPLOT_NOQC_GiantGenoOverlap = extractInfoForRegions(chrs=theseRegions$CHR,region.indices=theseRegions$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQCNone,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","GENO"))


#TOPLOT_QC3_GiantImpOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","IMP"))

#TOPLOT_QC3_GiantGenoOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_GENO$CHR,region.indices=theseRegionsGIANT_and_GENO$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","GENO"))


#TOPLOT_QC4_GiantImpOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC4,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","IMP"))

#TOPLOT_QC4_GiantGenoOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_GENO$CHR,region.indices=theseRegionsGIANT_and_GENO$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC4,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","GENO"))


# TEST2 = extractInfoForRegions(chrs=theseRegions$CHR[2],region.indices=theseRegions$region.index[2],columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC,scaleFactor=1,prior=0.2)

# TEST = extractInfoForRegions(chrs=theseRegions$CHR[1:2],region.indices=theseRegions$region.index[1:2],columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC,scaleFactor=scaleFactor,prior=0.2)

# TEST3 = extractInfoForRegions(chrs=theseRegions$CHR[2],region.indices=theseRegions$region.index[2],columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC,scaleFactor=scaleFactor,prior=0.2,allOverlapping=TRUE)

#save(TOPLOT_NOQC,TOPLOT_QC2,TOPLOT_QC3,TOPLOT_QC4,theseRegions,scaleFactor,prior,file=paste0(outDir,"/",outName))


TOPLOT_QC3_GiantImpOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","IMP"))


TOPLOT_QC3_GiantImpOverlap_prior0.02 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=0.02,allOverlapping=TRUE,types=c("GIANT","IMP"))

TOPLOT_QC3_GiantImpOverlap_prior20 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=20,allOverlapping=TRUE,types=c("GIANT","IMP"))



#TOPLOT_QC3_GiantGenoOverlap = extractInfoForRegions(chrs=theseRegionsGIANT_and_GENO$CHR,region.indices=theseRegionsGIANT_and_GENO$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior,allOverlapping=TRUE,types=c("GIANT","GENO"))

TOPLOT_QC3 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior)

TOPLOT_QC3_prior0.02 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=0.02)

TOPLOT_QC3_prior20 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=20)


# ===> Very little impact from the prior!


TOPLOT_QC3_pvals = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC3,scaleFactor=scaleFactor,prior=prior,mustBeSignificant=TRUE)


# info>0.9
TOPLOT_QC5 = extractInfoForRegions(chrs=theseRegionsGIANT_and_IMP$CHR,region.indices=theseRegionsGIANT_and_IMP$region.index,columns=columns,MergedData=MergedData,MergedDataQC=MergedDataQC5,scaleFactor=scaleFactor,prior=prior)

# ===> Not much change.


#########################
# Plot a region and the effect sizes.


topXhits=nrow(theseRegionsGIANT_and_IMP)

#chr=regInfo$CHR;
#reg=regInfo$region.index



#PlotPosteriorsAndEffectSizes(X = TOPLOT_NOQC,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",topXhits,"-",qcNonetitle),plotdir=plotdir,nRegions=20,effectSizesOnly=TRUE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_QC2,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",topXhits,"-",qc2title),plotdir=plotdir,nRegions=20,effectSizesOnly=TRUE)

PlotPosteriorsAndEffectSizes(X = TOPLOT_QC3,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title),plotdir=plotdir,nRegions=c(1:20),effectSizesOnly=FALSE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_QC4,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",topXhits,"-",qc4title),plotdir=plotdir,nRegions=20,effectSizesOnly=TRUE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_NOQC_allStudiesOnly,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",topXhits,"-noqc-overlappingOnly"),plotdir=plotdir,nRegions=20,effectSizesOnly=TRUE)

# NO QC overlaps

#PlotPosteriorsAndEffectSizes(X = TOPLOT_NOQC_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",nrow(theseRegionsGIANT_and_IMP),"-",qcNonetitle,"-GIANT-IMP-overlappingOnly"),plotdir=plotdir,nRegions=20,types=c("GIANT","IMP"),cols=c(giantColour,impColour),shapes=c(24,21),effectSizesOnly=TRUE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_NOQC_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,plotName = paste0("GIANT-top",nrow(theseRegionsGIANT_and_GENO),"-",qcNonetitle,"-GIANT-GENO-overlappingOnly"),plotdir=plotdir,nRegions=20,types=c("GIANT","GENO"),cols=c(giantColour,genoColour),shapes=c(24,23),effectSizesOnly=TRUE)

#  QC3 overlaps

PlotPosteriorsAndEffectSizes(X = TOPLOT_QC3_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-IMP-top",nrow(theseRegionsGIANT_and_IMP),"-",qc3title,"-GIANT-IMP-overlappingOnly"),plotdir=plotdir,nRegions=c(50:100),types=c("GIANT","IMP"),cols=c(giantColour,impColour),shapes=c(24,21),effectSizesOnly=FALSE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_QC3_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,plotName = paste0("GIANT-top",nrow(theseRegionsGIANT_and_GENO),"-",qc3title,"-GIANT-GENO-overlappingOnly"),plotdir=plotdir,nRegions=20,types=c("GIANT","GENO"),cols=c(giantColour,genoColour),shapes=c(24,23),effectSizesOnly=TRUE)

#  QC4 overlaps

#PlotPosteriorsAndEffectSizes(X = TOPLOT_QC4_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-top",nrow(theseRegionsGIANT_and_IMP),"-",qc4title,"-GIANT-IMP-overlappingOnly"),plotdir=plotdir,nRegions=20,types=c("GIANT","IMP"),cols=c(giantColour,impColour),shapes=c(24,21),effectSizesOnly=TRUE)

#PlotPosteriorsAndEffectSizes(X = TOPLOT_QC4_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,plotName = paste0("GIANT-top",nrow(theseRegionsGIANT_and_GENO),"-",qc4title,"-GIANT-GENO-overlappingOnly"),plotdir=plotdir,nRegions=20,types=c("GIANT","GENO"),cols=c(giantColour,genoColour),shapes=c(24,23),effectSizesOnly=TRUE)




#########################
# Compute and plot credible sets.


###
#X=TOPLOT_NOQC_GiantImpOverlap
#theseTypes=c("GIANT","IMP")
###


#test = summarisePosteriorComparison(TOPLOT_NOQC,theseRegions,theseTypes=c("GIANT","IMP"))


#CredibleSets_NOQC_GiantImpOverlap = summarisePosteriorComparison(TOPLOT_NOQC_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_NOQC_GiantGenoOverlap = summarisePosteriorComparison(TOPLOT_NOQC_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))

#CredibleSets_NOQC_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_NOQC,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_NOQC_GIANTvsGENO = summarisePosteriorComparison(TOPLOT_NOQC,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))

# Apply QC3
CredibleSets_QC3_GiantImpOverlap = summarisePosteriorComparison(TOPLOT_QC3_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

CredibleSets_QC3_GiantImpOverlap_prior0.02 = summarisePosteriorComparison(TOPLOT_QC3_GiantImpOverlap_prior0.02,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

CredibleSets_QC3_GiantImpOverlap_prior20 = summarisePosteriorComparison(TOPLOT_QC3_GiantImpOverlap_prior20,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_QC3_GiantGenoOverlap = summarisePosteriorComparison(TOPLOT_QC3_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))

CredibleSets_QC3_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC3,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

exclb = is.na(TOPLOT_QC3$IMP..MyNewPosteriors) & (!is.na(TOPLOT_QC3$GIANT..MyNewPosteriors))
CredibleSets_QC3b_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC3[!exclb,],theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

CredibleSets_QC3_prior0.02_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC3_prior0.02,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

CredibleSets_QC3_prior20_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC3_prior20,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_QC3_pvals_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC3_pvals,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))


#CredibleSets_QC3_GIANTvsGENO = summarisePosteriorComparison(TOPLOT_QC3,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))


# Apply QC4
#CredibleSets_QC4_GiantImpOverlap = summarisePosteriorComparison(TOPLOT_QC4_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_QC4_GiantGenoOverlap = summarisePosteriorComparison(TOPLOT_QC4_GiantGenoOverlap,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))

#CredibleSets_QC4_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC4,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

#CredibleSets_QC4_GIANTvsGENO = summarisePosteriorComparison(TOPLOT_QC4,theseRegions=theseRegionsGIANT_and_GENO,theseTypes=c("GIANT","GENO"))

# Apply QC5
CredibleSets_QC5_GIANTvsIMP = summarisePosteriorComparison(TOPLOT_QC5,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))


#### Get the data for all markers and whether they're in the 95% credible set

inCredibleSets_QC3_GIANTvsIMP = getIn95CredibleSet(TOPLOT_QC3,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))

inCredibleSets_QC3_GiantImpOverlap = getIn95CredibleSet(TOPLOT_QC3_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,theseTypes=c("GIANT","IMP"))


save(TOPLOT_QC3,TOPLOT_QC3_GiantImpOverlap,
     CredibleSets_QC3_GiantImpOverlap,CredibleSets_QC3_GiantImpOverlap_prior0.02,CredibleSets_QC3_GiantImpOverlap_prior20,
     CredibleSets_QC3_GIANTvsIMP,CredibleSets_QC3_prior0.02_GIANTvsIMP,CredibleSets_QC3_prior20_GIANTvsIMP,
     CredibleSets_QC5_GIANTvsIMP,
     theseRegionsGIANT_and_IMP,
     AllRegions,
     QC3,
     prior,
     scaleFactor,dataFile,
     inCredibleSets_QC3_GIANTvsIMP,inCredibleSets_QC3_GiantImpOverlap,
     file = paste0(outDir,"/CredibleSetsAnalysis.QC3.RData"))

# Plotted to here.







############
# For paper.


PlotPosteriorsAndEffectSizes(X = TOPLOT_QC3,theseRegions=theseRegionsGIANT_and_IMP,plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title),plotdir="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/forPaper/",nRegions=c(1:20),effectSizesOnly=TRUE,cex.lab=1.5,cex.axis=1.3,cex=1.3,mar=c(5, 5, 4, 2) + 0.1)



#t = inCredibleSets[[1]][inCredibleSets[[1]]$region_ID=="1_1",]


#### TEST THE PLOTS

#plotName = paste0("GIANT-top",topXhits,"-",qcNonetitle,"..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_NOQC_GIANTvsGENO,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qcNonetitle,"..GIANTvsIMP")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_NOQC_GIANTvsIMP,theseRegions,type1="GIANT",type2="IMP",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qcNonetitle,"-GIANT-GENO-overlappingOnly..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_NOQC_GiantGenoOverlap,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qcNonetitle,"-GIANT-IMP-overlappingOnly..GIANTvsIMP")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_NOQC_GiantImpOverlap,theseRegions,type1="GIANT",type2="IMP",pch=16)

#dev.off()



### After applying QC3

#plotName = paste0("GIANT-top",topXhits,"-",qc3title,"..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC3_GIANTvsGENO,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"..GIANTvsIMP")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#DATA = extractCredibleData(CredibleSets_QC3_GIANTvsIMP)
PlotCredibleSets(CredibleSets_QC3_GIANTvsIMP,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "Size-1 numbers"
#[1] 78
#[1] 85
#[1] "median size"
#[1] 6
#[1] 8
#[1] "median size per snp"
#[1] 0.04700855
#[1] 0.01038961
#[1] "total number markers"
#[1] 106263
#[1] 768502
#[1] "sizes"

dev.off()

# Only allowing SNPS in GIANT if they're in IMP
plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"b..GIANTvsIMP")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3b_GIANTvsIMP,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

dev.off()

# DIFFERENT PRIORS
plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"..GIANTvsIMP.prior0.02")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3_prior0.02_GIANTvsIMP,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "Size-1 numbers"
#[1] 75
#[1] 86
#[1] "median size"
#[1] 6
#[1] 8
#[1] "median size per snp"
#[1] 0.05102041
#[1] 0.01075269
#[1] "total number markers"
#[1] 106263
#[1] 768502
#[1] "sizes"
dev.off()

plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"..GIANTvsIMP.prior20")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3_prior20_GIANTvsIMP,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "Size-1 numbers"
#[1] 78
#[1] 85
#[1] "median size"
#[1] 6
#[1] 8
#[1] "median size per snp"
#[1] 0.04700855
#[1] 0.01038961
#[1] "total number markers"
#[1] 106263
#[1] 768502
#[1] "sizes"

dev.off()



#plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"-onlysignif..GIANTvsIMP")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC3_pvals_GIANTvsIMP,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qc3title,"-GIANT-GENO-overlappingOnly..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC3_GiantGenoOverlap,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"-GIANT-IMP-overlappingOnly..GIANTvsIMP")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3_GiantImpOverlap,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "Size-1 numbers"
#[1] 76
#[1] 123
#[1] "median size"
#[1] 6
#[1] 4
#[1] "median size per snp"
#[1] 0.04761905
#[1] 0.03921569
#[1] "total number markers"
#[1] 105421
#[1] 105421
#[1] "sizes"

dev.off()


plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"-GIANT-IMP-overlappingOnly..GIANTvsIM.prior0.02")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3_GiantImpOverlap_prior0.02,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "GIANT"
#[1] "IMP"

#[1] "Size-1 numbers"
#[1] 73
#[1] 121
#[1] "median size"
#[1] 6
#[1] 4
#[1] "median size per snp"
#[1] 0.05102041
#[1] 0.03985507
#[1] "total number markers"
#[1] 105421
#[1] 105421
#[1] "sizes"
dev.off()


plotName = paste0("GIANT-IMP-top",topXhits,"-",qc3title,"-GIANT-IMP-overlappingOnly..GIANTvsIMP.prior20")
pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

PlotCredibleSets(CredibleSets_QC3_GiantImpOverlap_prior20,theseRegions=theseRegionsGIANT_and_IMP,type1="GIANT",type2="IMP",pch=16)

#[1] "GIANT"
#[1] "IMP"

#[1] "Size-1 numbers"
#[1] 76
#[1] 122
#[1] "median size"
#[1] 6
#[1] 4
#[1] "median size per snp"
#[1] 0.04761905
#[1] 0.03921569
#[1] "total number markers"
#[1] 105421
#[1] 105421

dev.off()



# OVERLAPPING
DATA0.2 = extractCredibleData(CredibleSets_QC3_GiantImpOverlap)$sizes[[16]]
DATA20 = extractCredibleData(CredibleSets_QC3_GiantImpOverlap_prior20)$sizes[[16]]
DATA20.02 = extractCredibleData(CredibleSets_QC3_GiantImpOverlap_prior0.02)$sizes[[16]]

#prior=20  one region changed in credible set size in both studies, by one marker.
#prior=0.02 538 one region changed in credible set size in both studies, by one marker.

# ALL
DATAa0.2 = extractCredibleData(CredibleSets_QC3_GIANTvsIMP)$sizes[[16]]
DATAa20 = extractCredibleData(CredibleSets_QC3_prior20_GIANTvsIMP)$sizes[[16]]
DATAa0.02 = extractCredibleData(CredibleSets_QC3_prior0.02_GIANTvsIMP)$sizes[[16]]

change = (DATAa0.02[1,]-DATAa0.2[1,])/DATAa0.2[1,]
diff = (DATAa0.02[1,]-DATAa0.2[1,])/DATAa0.2[1,]

#max prop difference 0.00401186,

################################################################
#### TESTING LEGEND
pdf("plots/test.pdf")

#these = sample(1:nrow(TOPLOT_QC3),2000,replace=FALSE)

plot(log(TOPLOT_QC3$IMP..MyNewPosterior[these]),TOPLOT_QC3$IMP..INFO[these])
plot(TOPLOT_QC3$IMP..MyNewPosterior[these],TOPLOT_QC3$IMP..INFO[these])

oo = abind(sapply(1:length(my_breaks),function(w) rep(w,my_breaks[w])),along=1)
plot(rep(1,length(oo)),oo,xlab=NA,axes=FALSE,ylab=NA,col=add.alpha("blue",0.1),main="Count",pch=16,cex=8,xpd=NA)
text(x=rep(1,length(my_breaks)),y=1:length(my_breaks),labels=my_breaks,cex=2,xpd=NA,adj=0.5)

dev.off()

### After applying QC4

#plotName = paste0("GIANT-top",topXhits,"-",qc4title,"..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC4_GIANTvsGENO,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qc4title,"..GIANTvsIMP")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC4_GIANTvsIMP,theseRegions,type1="GIANT",type2="IMP",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qc4title,"-GIANT-GENO-overlappingOnly..GIANTvsGENO")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC4_GiantGenoOverlap,theseRegions,type1="GIANT",type2="GENO",pch=16)

#dev.off()


#plotName = paste0("GIANT-top",topXhits,"-",qc4title,"-GIANT-IMP-overlappingOnly..GIANTvsIMP")
#pdf(paste0(plotdir,"/plots/credible-sets.",plotName,".pdf"))

#PlotCredibleSets(CredibleSets_QC4_GiantImpOverlap,theseRegions,type1="GIANT",type2="IMP",pch=16)

#dev.off()




########################################################
# To plot for paper
########################################################






pdf("plots/GIANT-top-TEST.pdf")
         # barplot of the sizes
        specialBinsBarplot(x[1,],main=paste0("Credible set sizes in ",sourceNames[type1]," for ",nRegs," regions at level ",b))
        specialBinsBarplot(x[2,],main=paste0("Credible set sizes in ",sourceNames[type2]," for ",nRegs," regions at level ",b))
        
        specialBinsBarplot2(x,main=paste0("Size of ",b*100,"% credible set in ",sourceNames[type1]," and ",sourceNames[type2]," for ",nRegs," regions"),ylab="Number of regions",xlab=paste0("Number of markers in ",b*100,"% credible set"),col=colours,legend=sourceNames[c(type1,type2)])


dev.off()






#### BELOW NOT USED

pdf(paste0(plotdir,"/plots/posterior-qqplot.",plotName,".pdf"))
x = TOPLOT_NOQC_GiantImpOverlap$GIANT..MyNewPosteriors
y = TOPLOT_NOQC_GiantImpOverlap$IMP..MyNewPosteriors
keep = !is.na(y)& !is.na(x)

qqplot(x[keep],y[keep],xlab="GIANT",ylab="UKB Imputed")
abline(0,1,col="red",lty=3)

x = TOPLOT_NOQC$GIANT..MyNewPosteriors
y = TOPLOT_NOQC$GENO..MyNewPosteriors
keep = !is.na(y)& !is.na(x)

qqplot(x[keep],y[keep],xlab="GIANT",ylab="UKB Genotyped")
abline(0,1,col="red",lty=3)

dev.off()

#

#X[X$SNP2=="03:141105570_A_G",]
#save(MergedData,AllRegions,regionFileBase,dataFiles,file=paste0(outDir,"/",outName))



