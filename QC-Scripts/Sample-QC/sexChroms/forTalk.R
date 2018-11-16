source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/sexChroms/hmm-function.R')
source("/well/donnelly/ukbiobank_project_8874/clare/sex/scripts/commonScripts/plotGenesFunction.R")


#args = commandArgs(trailingOnly=T)

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-samples","b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-oddSamples.txt","-out","b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v3")

## b1__b11-b001__b095-sexchrom-sampleqc-sexCheck uses myForwardsBackwards.OLD, and has a fixed transition matrix.
## b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2 uses myForwardsBackwards in the hmm. This includes a better transition matrix, as well as options for updating transition probabilities.
## b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v3 as v2 but using fixed parameters across all snp as the starting positions + some noise for the means (variances equal). No updating of transition parameters...

print(args)
h = args[-c(which(args%in%c("-in","-out","-samples")),1+which(args%in%c("-in","-out","-samples")))]
for(helperScript in h){
    source(helperScript)
}

setwd("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/sexChroms")

outFile = args[which(args=="-out")+1]
sampleFile = args[which(args=="-samples")+1]

# get snp annotations (position etc)
UKBB = ukbiobank.ps2snp("sexchrom")
UKBL = ukbileve.ps2snp("sexchrom")
UKBL = tbl_df(UKBL)
UKBB = tbl_df(UKBB)

ps2snp = merge(UKBL, UKBB, all = TRUE)
ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]

# get sample info
otherInfo = read.multiple.batch.info(batchInfoFields)
otherInfo = tbl_df(otherInfo)

# read in cnv summary data
load(paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries.RData"),verbose=TRUE)
meanL2rY = outData[["log2ratio"]][,"means.Y"]
meanL2rX = outData[["log2ratio"]][,"means.X"]

# get list of SNPs used in cnv summaries (and to be used in HMM)
snpsToInclude = read.table(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc.bim"),header=FALSE)[,2]

# read in recombrates
recombrates = list()
for(i in c("X_nonPAR","X_PAR1","X_PAR2")){
    print(i)
    r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
    recombrates[[i]] = r
}
names(recombrates) = c("X","XPAR1","XPAR2")
recombRatesX = Reduce(rbind,recombrates)




########################################
# RUN HMM ON SET OF EXEMPLAR SAMPLES WITH EM TO LEARN PARAMETERS ()

# read in list of samples to run HMM on
oddSamples = read.table(sampleFile)[,1]
controlSamples = read.table("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-controlSamples.txt")[,1]

# read initial classifications for odd samples
load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-initialClassificationsOddSamples.Rdata",verbose=TRUE)

# read in raw cnv data for chromosome X & Y
snps = ps2snp$ProbeSetID[(ps2snp$AffySNPID%in%snpsToInclude) & (ps2snp$Chromosome %in% c(23,24,25)) ]

dataSubset = read.cnv.hdf5(inds=c(oddSamples,controlSamples),snps=snps,type="log2ratio",otherInfo=otherInfo)

log2RatiosAll = abind(dataSubset,along=2)  # careful: this only works with a limited number of samples

sampleOrderAll = colnames(log2RatiosAll)
snpOrderAll = rownames(log2RatiosAll)
snpChromAll = ps2snp$Chromosome[match(snpOrderAll,ps2snp$ProbeSetID)]
snpPosAll = ps2snp$Position[match(snpOrderAll,ps2snp$ProbeSetID)]


########################################
# LOAD GENES DATA

genes = load.genes(refGeneFile,fieldNames=fieldNames)
exons = get.exons(genes,chromosome = "X")


########################################
# LOAD HMM FOR ALL SAMPLES ON CHROM X.
    
#filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-155233115"
filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2-hmmResults-chrX-60425-to-155233115"
load(paste0(filename,".Rdata"),verbose=TRUE)

# read initial classifications for odd samples based on X-wide averages
load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-initialClassificationsOddSamples.Rdata",verbose=TRUE)

# read in cnv summary data
load(paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries.RData"),verbose=TRUE)
meanL2rY = outData[["log2ratio"]][,"means.Y"]
meanL2rX = outData[["log2ratio"]][,"means.X"]



################## Plot chrom X by chunks and order by cnv values
h5 = paste0(filename,"-byBatch.h5")
load(paste0(filename,"-hclustOrderPARoddSamples.RData"),verbose=TRUE) # get order of odd samples
load(paste0(filename,"-PAR1ExtendedSamples.RData"),verbose=TRUE) # get list of par extension samples


# just read in samples with odd values at 86-96 at start and end of X chromosome. samplesToPlot is from above.

samplesToExtract = unique(c(hclustOrderSamples,hclustOrderParExtendSamples))

theseSnps = snpOrder[snpPos <= (5*1e6) ]
dataHMMparExt = read.hmm.hdf5(h5,inds=samplesToExtract,snps=theseSnps,otherInfo=otherInfo)

                                        # find a good ordering
start = 0; chunkWidth = 5*1e6/1000000

W = transformHMMtoCalls(dataHMMparExt,threshold=0.98)
nNot0 = colSums(W!=0)
frac1 = colSums(W==1)/nNot0
frac2 = colSums(W==2)/nNot0
frac3 = colSums(W==3)/nNot0
frac0 = colSums(W==0)/nrow(W)
keepGood = colnames(W) %in% otherInfo$PIID # exclude references
#keep = (frac3 > 0.1) | ( (frac1 > 0.1 )&(frac1 < 0.9) ) | ( (frac2 > 0.1 )&(frac2 < 0.9) )
keepW = W[,keepGood]
samplesToPlotHClustOrder = hclust.order(keepW)
samplesToPlot = c(hclustOrderSamples,samplesToPlotHClustOrder[!samplesToPlotHClustOrder%in%hclustOrderSamples])
keepInds = colnames(W)[( frac0 < 0.5 )]

samplesToPlotKEEP = samplesToPlot


#plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    

plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtAndOddSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    


plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot[samplesToPlot%in%hclustOrderParExtendSamples],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    


start = 2.5; chunkWidth=0.5
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]

plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot[samplesToPlot%in%hclustOrderParExtendSamples],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2000)    


plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtAndOddSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2000,extraPlotting=FALSE)


plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtAndOddSamples-meanSorted-withGenes-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=TRUE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2500,addGenes=TRUE,cex=0.8)
abline(v=par1,lwd=5,lty=3)

plot.genes(chromosome = formatChrom("X"),region=c(start,start+chunkWidth)*1e6,local.genes=genes,exons=exons,label.cex=0.8,height_in_inches = 1,xaxt="n")
dev.off()


start = 2.65; chunkWidth=0.2
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtAndOddSamples-meanSorted-withGenes-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=TRUE,updateLegend=list("0"="uncertain"),pngRes=150,Width=4000,Height=2000,addGenes=TRUE,cex=0.1)

plot.genes(chromosome = formatChrom("X"),region=c(start,start+chunkWidth)*1e6,local.genes=genes,exons=exons,label.cex=0.8,height_in_inches = 1,xaxt="n")
dev.off()


############################################################
##### PICK A SELECTION OF SAMPLES FOR TALK (MAX two plots)!!

males = samplesToPlot[samplesToPlot%in%otherInfo$PIID[otherInfo$Inferred.Gender=="M"]]
females = samplesToPlot[samplesToPlot%in%otherInfo$PIID[otherInfo$Inferred.Gender=="F"]]
lowY = names(meanL2rY)[meanL2rY<(-0.25)&(meanL2rX < -0.2)]
weird = samplesToPlot[samplesToPlot%in%c(par1init,par3init,par3init2,x3init,y3init,lowY)]

ORDER = c(650:880,1100:1330,1600:1900,1365:1477)
ORDER = c(650:880,1100:1330,1600:1900,1380:1477)


start = 2.5; chunkWidth=0.5
#start = 0; chunkWidth=3
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
plotHMM(dataHMMparExt,samplesToPlot=samplesToPlot[ORDER][!samplesToPlot[ORDER]%in%weird],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtAndOddSamplesFORTALK-meanSorted-withGenes-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=TRUE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2500,addGenes=TRUE,cex=0.8)
abline(v=par1,lwd=5,lty=3)
abline(v=c(2694151,2694702,2808549,2809097),col="purple")

plot.genes(chromosome = formatChrom("X"),region=c(start,start+chunkWidth)*1e6,local.genes=genes,exons=exons,label.cex=0.8,height_in_inches = 1,xaxt="n")

dev.off()




# read in l2r for a set of individuals (male, female)

#examples = c(males[!males%in%weird][1:10],females[!females%in%weird][1:10])

examples = samplesToPlot[c(300,700,900,1200,1201,1400:1410,1500,1600)]


log2RatiosExamp = read.cnv.hdf5(inds=examples,snps=snpOrder[snpPos <= (5*1e6) ],type="log2ratio",otherInfo=otherInfo)
log2RatiosExamp = abind(log2RatiosExamp,along=2)


start = 2; chunkWidth=2
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
plotL2RIndividual(log2RatiosExamp,dataHMMparExt,examples[18],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))





##########################################
# Check for double deletion
##########################################
# frac0 etc. are from oddSamplesHMM in sex-hmm.R

#theseSnps = snpOrder[( (snpPos >= (par1 - 0.25*1e6) ) & ( snpPos <= (par1+0.25*1e6) ) )]
#oddSamplesHMM = read.hmm.hdf5(h5,inds=NULL,snps=theseSnps,otherInfo=otherInfo)
#W = transformHMMtoCalls(oddSamplesHMM,threshold=0.98)


toCheck = names(which((frac0b>0.3)&(frac0<0.05)))
dataHMMzeros = read.hmm.hdf5(h5,inds=toCheck,snps=snpOrder[snpPos <= (5*1e6) ],otherInfo=otherInfo)
sorted = names(sort(frac1+frac1b))
thiOrder = rev(sorted[sorted%in%toCheck])

start = 2.65; chunkWidth=0.2
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
plotHMM(dataHMMzeros,samplesToPlot=thiOrder,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parZeroSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2000,cex=0.4)


# or read in l2r for everyone around the PAR.
log2RatiosPAR = read.cnv.hdf5(inds=NULL,snps=snpOrder[(snpPos <= (2.85*1e6))&(snpPos >= (2.65*1e6)) ],type="log2ratio",otherInfo=otherInfo)
log2RatiosPAR = abind(log2RatiosPAR,along=2)


         
l2r1 = colMeans(log2RatiosPAR[rownames(log2RatiosPAR)%in%snpOrder[snpPos <= par1],])
l2r2 = colMeans(log2RatiosPAR[rownames(log2RatiosPAR)%in%snpOrder[(snpPos > par1)&(snpPos < 2.8*1e6)],])
l2r3 = colMeans(log2RatiosPAR[rownames(log2RatiosPAR)%in%snpOrder[(snpPos > 2.81*1e6)],])

png("test.png")
plot(l2r1,l2r2)
dev.off()

png("test2.png")
plot(l2r2,l2r3)
dev.off()

png("test3.png")
plot(l2r1,l2r3)
dev.off()

unusual = names(l2r2)[l2r2<(-0.75)]
unusual2 = names(l2r2)[l2r2<(-1.5)]
unusual1 = unusual[!unusual%in%unusual2]

log2RatiosLowPARExt = read.cnv.hdf5(inds=unusual,snps=snpOrder[snpPos <= (5*1e6) ],type="log2ratio",otherInfo=otherInfo)
log2RatiosLowPARExt = abind(log2RatiosLowPARExt,along=2)


#### Look for 'odd samples' based on PAR1 boundary as well.
theseSnps = snpOrder[( (snpPos >= (par1 - 0.25*1e6) ) & ( snpPos <= (par1+0.25*1e6) ) )]
oddSamplesHMM = read.hmm.hdf5(h5,inds=NULL,snps=theseSnps,otherInfo=otherInfo)


start = 2.65; chunkWidth=0.2
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]

for(i in 1:20){
    plotL2RIndividual(log2RatiosLowPARExt,oddSamplesHMM,unusual2[i],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    
}
    
start = par1/1e6 - 0.25; chunkWidth=0.5
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]

plotHMM(oddSamplesHMM,samplesToPlot=c(unusual1,unusual2,hclustOrderParExtendSamples[1710:1720]),snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-parExtZeroSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.999,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"),pngRes=150,Width=2000,Height=2000,cex=0.4)



########################################
# PLOT the l2r and inferred path for some strange individuals

# read in l2r for these individuals
log2RatiosKEEP = read.cnv.hdf5(inds=toCheck,snps=snpOrder[snpPos <= (5*1e6) ],type="log2ratio",otherInfo=otherInfo)
log2RatiosKEEP = abind(log2RatiosKEEP,along=2)

start = 2.65; chunkWidth=0.2
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]

for(i in 5:20){
    plotL2RIndividual(log2RatiosKEEP,dataHMMzeros,thiOrder[i],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    
}
    



##########################################
# QC failures by position
##########################################
fails = read.SNPQC.files(type="sexchrom")

theseSnps =  ps2snp$Chromosome%in%c(23,25)
pos = ps2snp$Position[theseSnps]
name = ps2snp$AffySNPID[theseSnps]

anyFail = name %in% unique(unlist(fails))
mfFail = name %in% fails$mfSNPs
hweFail = name %in% fails$batchHweSNPs


stepSize = 100000 # 
window = seq(0,max(pos),by=stepSize)
windowSize = 1000000
nSNPS = sapply(1:length(window), function(i){
    inWindow = pos%in%c(window[i]:(window[i]+windowSize))
    nSnps = sum(inWindow)
})

nHWE = sapply(1:length(window), function(i){
    inWindow = pos%in%c(window[i]:(window[i]+windowSize))
    nFails = sum(inWindow&hweFail)
    return(nFails)
})
nMF = sapply(1:length(window), function(i){
    inWindow = pos%in%c(window[i]:(window[i]+windowSize))
    nFails = sum(inWindow&mfFail)
    return(nFails)
})
nAny = sapply(1:length(window), function(i){
    inWindow = pos%in%c(window[i]:(window[i]+windowSize))
    nFails = sum(inWindow&anyFail)
    return(nFails)
})


mids = window + windowSize/2
    
png(paste0("plots/X-chromosome-qc-fail-position-density.png"),bg="transparent",width=1000,height=500)
plot(mids,nAny/nSNPS,pch=16,ylim=c(0,0.1))
points(mids,nHWE/nSNPS,pch=16,col="red",cex=0.1)
points(mids,nMF/nSNPS,pch=16,col="blue",cex=0.1)
lines(mids,nHWE/nSNPS,pch=16,col="red")
lines(mids,nMF/nSNPS,pch=16,col="blue")

dev.off()


mfPvals = lapply(c(ukbileve.batches(),ukbiobank.batches()),function(b){
    h = read.table(paste0("/well/ukbiobank/expt/V2_QCed.SNP-QC/data/",b,"/tests/V2_QCed.",b,".sexchrom.malefemale.txt"),stringsAsFactors=FALSE,header=TRUE)
})
mfPvalsMat = abind(mfPvals,along=1,force.array=FALSE)

minP = aggregate(by=list(mfPvalsMat$AffySNPID),x=mfPvalsMat$pFRQ_ceu,FUN="min")

hwePvals = lapply(c(ukbileve.batches(),ukbiobank.batches()),function(b){
    h = read.table(paste0("/well/ukbiobank/expt/V2_QCed.SNP-QC/data/",b,"/tests/V2_QCed.",b,".sexchrom.hwe.txt"),stringsAsFactors=FALSE,header=TRUE)
})
hwePvalsMat = abind(hwePvals,along=1,force.array=FALSE)

minPhwe = aggregate(by=list(hwePvalsMat$AffySNPID),x=hwePvalsMat$pFRQ_ceu,FUN="min")


png(paste0("plots/X-chromosome-qc-fail-position-pvalue.png"),bg="transparent",width=1000,height=500)
plot(ps2snp$Position[match(minP[,1],ps2snp$AffySNPID)],-log10(minP[,2]))
points(ps2snp$Position[match(minPhwe[,1],ps2snp$AffySNPID)],-log10(minPhwe[,2]),col="red")

these = which(-log10(minP[,2])>150)
points(ps2snp$Position[match(minP[these,1],ps2snp$AffySNPID)],rep(150,length(these)),col="blue")

dev.off()


##########################################
# Some cluster plots ==> QC failures and successes
##########################################

# SNPs are extracted from SNP-QC file.
# example of odd outside par1: Affx-34795806.

freqs = read.genotyped.maf( snpFrequencyData= "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr" )

toTest = freqs[freqs$CHR==23,]
pvalBig = minP[,1][which(-log10(minP[,2])<5)]

good = intersect(ps2snp$dbSNPRSID[match(pvalBig,ps2snp$AffySNPID)],toTest$SNP[toTest$MAF>0.1])

head(ps2snp[match(good,ps2snp$dbSNPRSID),])

# RELEASED DATA
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-34795806 plots/cluster_plots_gender -b Batch_b001,Batch_b002,Batch_b003 -released -onepage -nossp -sex -rNumber 102 &

# raw data snp which fails male vs female (un qc'd)

Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-35119527,Affx-37677012,Affx-37677016,Affx-35145751,Affx-37681251 plots/cluster_plots_gender -b Batch_b001,Batch_b002,Batch_b003 -onepage -nossp -sex -rNumber 103 &

                                        # good snp
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-34462427,Affx-34463634 plots/cluster_plots_gender -b Batch_b001,Batch_b002,Batch_b003 -onepage -nossp -sex -rNumber 104 &
    
