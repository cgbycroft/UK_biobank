source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/sexChroms/hmm-function.R')

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

#### NOTE: ABOVE HERE IS NECESSARY FOR LATER PARTS TOO ################

# run HMM based on an initialisation of the parameters.
# Parameters required are:
# mean and variance of l2ratios of each state in 0,1,2,3 copies; and for each SNP.


# Expectation: For each copy number, take a weighted average across all individuals, where weights are the probabilities of being in each state.

# Maximisation: Forwards/backwards algorithm using the computed emission probabilities.


# initialise with an obvious copy number state using sex (could be more sophisticated, but we'll see how we go)
############### What are we running it over?
log2Ratios = log2RatiosAll[snpChromAll%in%c(23,25) & (snpPosAll < 5*(10^6)),]  # X chromosome near PAR1
chrom="X"
############### 
log2Ratios = log2RatiosAll[snpChromAll%in%c(23,25) & (snpPosAll > (par2- 5*(10^6))),]  # X chromosome near PAR2
chrom="X"
###############
log2Ratios = log2RatiosAll[snpChromAll%in%c(23,25),]  # All of X!!!
chrom="X"
###############
log2Ratios = log2RatiosAll[snpChromAll%in%c(24),]  # Y chromosome
chrom="Y"
###############

sampleOrder = colnames(log2Ratios)
snpOrder = rownames(log2Ratios)
snpInfo = ps2snp[match(snpOrder,ps2snp$ProbeSetID),]
snpChrom = snpInfo$Chromosome
snpPos = snpInfo$Position

meanL2r = outData[["log2ratio"]][colnames(log2Ratios),]
# get male.female initial calls
males = otherInfo$Inferred.Gender[match(rownames(meanL2r),otherInfo$PIID)]=="M"
females = otherInfo$Inferred.Gender[match(rownames(meanL2r),otherInfo$PIID)]=="F"


## INITIALISE STATE PROBABILITIES (to get initial estimates of hmm parameters - i.e means and variances of log2ratios for each snp, and each copy number state)
stateProbabilities = array(0,dim=c(dim(log2Ratios)[1],dim(log2Ratios)[2],nStates=3))
# non-PAR X
stateProbabilities[snpChrom==23,sampleOrder%in%x1init,1] = 1  # one copy
stateProbabilities[snpChrom==23,sampleOrder%in%x2init,2] = 1  # two copies
stateProbabilities[snpChrom==23,sampleOrder%in%x3init,3] = 1  # three copies
# PAR1 X
stateProbabilities[(snpChrom==25) & (snpPos < 5e6),(!sampleOrder%in%par1init)&(!sampleOrder%in%par3init),2] = 1 # everyone has two copies
stateProbabilities[(snpChrom==25) & (snpPos < 5e6),sampleOrder%in%par1init,1] = 1 # except people who have 1 par
stateProbabilities[(snpChrom==25) & (snpPos < 5e6),sampleOrder%in%par3init,3] = 1 # three pars
# PAR2 X
stateProbabilities[(snpChrom==25) & (snpPos > 5e6),(!sampleOrder%in%par1init2)&(!sampleOrder%in%par3init2),2] = 1 # everyone has two copies
stateProbabilities[(snpChrom==25) & (snpPos > 5e6),sampleOrder%in%par1init2,1] = 1 # except people who have 1 par
stateProbabilities[(snpChrom==25) & (snpPos > 5e6),sampleOrder%in%par3init2,3] = 1 # three pars
# Y. NOTE: Y chrom copy number is one less than the value
stateProbabilities[snpChrom==24,(!sampleOrder%in%y2init)&(!sampleOrder%in%y3init),1] = 1  # zero copies (females)
stateProbabilities[snpChrom==24,sampleOrder%in%y2init,2] = 1  # two copies
stateProbabilities[snpChrom==24,sampleOrder%in%y3init,3] = 1  # three copies


# check they all sum to 1
sumTest = apply(stateProbabilities,2,rowSums)
sum(sumTest!=1)

set.seed(123456) # to account for the random noise added at initialisation. actually not used.

# set up transition matrix (it is later scaled by bp distance)
stayProb = 0.999
trans=diag(stayProb,dim(stateProbabilities)[3])  # make it hard to change state
trans[1,2:3] = (1-stayProb)*c(0.5,0.5)
trans[2,c(1,3)] = (1-stayProb)*c(0.5,0.5)
trans[3,c(1,2)] = (1-stayProb)*c(0.5,0.5)

# set filename
filename = paste0(outFile,"-hmmResults-chr",chrom,"-",min(snpPos),"-to-",max(snpPos))

# plot initial parameter estimates
ExpectationInit = getSNPposteriors(l2r=log2Ratios,stateProbabilities)       
png(paste0("plots/",filename,"-InitialParametersCheck-hist.png"),width=1000,height=1000,res=150)
plotParametersHist(hmmToPlot=NULL,modelParameters=ExpectationInit,snpsToPlot=NULL)
dev.off()

# all snps the same, but with noise   
png(paste0("plots/",filename,"-InitialParametersCheck-histTEST.png"),width=1000,height=1000,res=150)
plotParametersHist(hmmToPlot=NULL,modelParameters=Expectation)
dev.off()


# Run HMM!!!!!!
# for versions v1,v2
#hmmResults = runHMM.em(log2Ratios,initialStateProbabilities=stateProbabilities,maxIterations=20,trans=trans,outIterations=2,tvdTol=1e-4,saveParameters=TRUE)
# for version v3 (don't update transitions)
hmmResults = runHMM.em(log2Ratios,initialStateProbabilities=stateProbabilities,maxIterations=20,trans=trans,outIterations=2,tvdTol=1e-4,saveParameters=TRUE,computeTransitions=FALSE)


# save a rounded version   
hmmResultsRound = hmmResults; hmmResultsRound$HMM[[1]] = round(hmmResultsRound$HMM[[1]],3)
save(hmmResultsRound,file=paste0(filename,".Rdata"))


###################
# plot results

load(paste0(filename,".Rdata"),verbose=TRUE)
#load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-4993450.Rdata",verbose=TRUE)
#load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-155233115.Rdata",verbose=TRUE)
#load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2-hmmResults-chrX-60425-to-155233115.Rdata",verbose=TRUE)

hmmToPlot = hmmResultsRound


# convergence
tvd = hmmToPlot$convStats[,1]
cert = hmmToPlot$convStats[,2]
iters = length(cert)
    
png(paste0("plots/",filename,"-convergenceStats-",iters,"-iterations.png"),width=1000,height=2000,res=150)
par(mfrow=c(2,1))
plot(1:length(tvd),tvd,xlab="iteration",ylab="mean tvd")
lines(1:length(tvd),tvd)
plot(1:length(cert),cert,xlab="iteration",ylab="mean certainty")
lines(1:length(cert),cert)
dev.off()

# parameter estimates by position
# Plot estimates of parameters after final iteration (actually this is the parameters estimated after the final iteration, not the set that were actually USED in the final iteration. This is true, as of August 1st 2016)
modelParameters = hmmToPlot$parameters
W = transformHMMtoCalls(hmmToPlot,threshold=0.98)
snpPositions = ps2snp$Position[match(rownames(W),ps2snp$ProbeSetID)]
t = getXY(W[,1:20],snpPositions)

chunkWidth=10
snpOrder = hmmToPlot$orderOfHMMSNPs
snpPos = ps2snp$Position[match(snpOrder,ps2snp$ProbeSetID)]
for( start in  seq(0,round(max(snpPos)/(10^6)),chunkWidth)){
    #start = 86
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
    png(paste0("plots/",filename,"-LearnedParameters-",start,"-",(start+chunkWidth),"MB.png"),width=2000,height=1000,res=150)
    plotParameters(hmmToPlot,snpsToPlot=theseSnps,colors=t$colLabels)
    dev.off()
}


# histograms of parameter estimates (all SNPs)
png(paste0("plots/",filename,"-LearnedParameters-hist.png"),width=1000,height=1000,res=150)
plotParametersHist(hmmResultsRound,snpsToPlot=NULL,colors=t$colLabels)
dev.off()



# plot hmm calls - basic - all samples
plotHMM(hmmToPlot,samplesToPlot=NULL,filename=paste0("plots/",filename),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# plot samples in oddSamples based on which category
# one x
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,x1init),filename=paste0("plots/",filename,"-putative1X"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# two x
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,x2init),filename=paste0("plots/",filename,"-putative2X"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# three x
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,x3init),filename=paste0("plots/",filename,"-putative3X"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 1 par1
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,par1init),filename=paste0("plots/",filename,"-putative1PAR1"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 3 par1
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,par3init),filename=paste0("plots/",filename,"-putative3PAR1"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 1 par2
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,par1init2),filename=paste0("plots/",filename,"-putative1PAR2"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 3 par2
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,par3init2),filename=paste0("plots/",filename,"-putative3PAR2"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 1 y
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,y2init),filename=paste0("plots/",filename,"-putative1Y"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# 2 y
plotHMM(hmmToPlot,samplesToPlot=intersect(oddSamples,y3init),filename=paste0("plots/",filename,"-putative2Y"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))


# all odd samples
plotHMM(hmmToPlot,samplesToPlot=oddSamples,filename=paste0("plots/",filename,"-allOddSamples"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# all control samples
plotHMM(hmmToPlot,samplesToPlot=controlSamples,filename=paste0("plots/",filename,"-controlSamples"),threshold=0.98,recombRates=recombRatesX,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))


# plot 5MB chunks of chromosome X
# get a sensible sample order for the plot, that puts similar individuals together, at least in the PAR1 and PAR2 regions
theseSnps = snpOrder[( (snpPos <= (5)*10^6) ) | ( snpPos > (par2-1e6) )]
# For v2 onwards, order by viterbi path...
#W = transformHMMtoCalls(hmmToPlot,samplesToPlot=oddSamples,snpsToPlot=theseSnps,threshold=0.98)
W = transformHMMtoCalls(hmmToPlot,samplesToPlot=oddSamples,snpsToPlot=theseSnps,threshold="viterbi")

Wna = W; Wna[Wna==0] = NA
d = dist(t(Wna)) # euclidian distance, ignoring uncertain values coded as NAs. takes a few seconds on 1500 x 1500
j = hclust(d)
hclustOrderSamples = colnames(W)[j$order]

# or load from version 2 (this is what I've done)
load(file=paste0("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2-hmmResults-chrX-60425-to-155233115-hclustOrderPARoddSamples.RData"),verbose=TRUE)


chunkWidth=10  # with of chunks in MB
for( start in seq(0,round(max(snpPos)/(10^6)),chunkWidth)){
    print(start)
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
    # threshold on snp posteriors
    plotHMM(hmmToPlot,samplesToPlot=hclustOrderSamples,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allOddSamples-PARsorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))
    # viterbi path
    #plotHMM(hmmToPlot,samplesToPlot=hclustOrderSamples,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allOddSamples-PARsorted-",start,"-",(start+chunkWidth),"MB"),threshold="viterbi",recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))  
}

# special for odd region around 90MB
start = 86
theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
plotHMM(hmmToPlot,samplesToPlot=hclustOrderSamples,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allOddSamples-PARsorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))  
plotHMM(hmmToPlot,samplesToPlot=hclustOrderSamples,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allOddSamples-PARsorted-",start,"-",(start+chunkWidth),"MB"),threshold="viterbi",recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))  


# print certain individuals
# number in hclustOrderSamples
indivIndex = c(100:120,350,500,700,900,1100)
for(start in c(0,86)){
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]    
    for( ind in indivIndex ){
        print(ind)
        plotL2RIndividual(log2Ratios,hmmToPlot,hclustOrderSamples[ind],snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allOddSamples-PARsorted-",start,"-",(start+chunkWidth),"MB-sampleIndex-",ind),threshold=0.98,title=ind,colors=colors,recombRates=recombRatesX)

    }
}


# what is the relationship between individual's l2r and the estimated parameters?
# exclude PAR and high homology region
in86 = hmmToPlot$orderOfHMMSNPs%in%ps2snp$ProbeSetID[ (ps2snp$Position > 8.8e7) & (ps2snp$Position < 9.25e7)] 
Xonly = hmmToPlot$orderOfHMMSNPs%in%ps2snp$ProbeSetID[ps2snp$Chromosome==23] 
keepSnps = Xonly&!in86
system("mkdir plots/mixtureResults/")

ind = 100;
for( ind in indivIndex ){
    print(ind)
    
    individual = hclustOrderSamples[ind]
    Y = log2Ratios[hmmToPlot$orderOfHMMSNPs[keepSnps],individual]
    x1 = hmmToPlot$parameters[[1]][keepSnps,1,1]
    x2 = hmmToPlot$parameters[[1]][keepSnps,2,1]

    Ws = mixtureTest(Y,x1,x2)
    w = Ws[[1]]; w2 = Ws[[2]]; w3 = Ws[[3]]

    residQP = Y - (w3$solution[1]*x1 + w3$solution[2]*x2)
    compareW1 = sum(residQP^2)/sum((Y - x1)^2) # if w=1
    compareW0 = sum(residQP^2)/sum((Y - x2)^2) # if w=0
         
    # linear regression
    png(paste0("plots/mixtureResults/",filename,"-allOddSamples-PARsorted-sampleIndex-",ind,"-",individual,"-1-2copiesMixtureLM.png"),width=1000,height=1000,res=150)
    plot(x1 - x2,Y - x2,ylab="indl2r - X2mean",xlab="X1mean - X2mean",pch=20,main=paste0("w = ",w$coefficients),xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),col=add.alpha("black",0.01))
    abline(0,1,lty=3,col="gray")
    abline(0,w$coefficients,col="red")
    dev.off()

    # quad programming 
    png(paste0("plots/mixtureResults/",filename,"-allOddSamples-PARsorted-sampleIndex-",ind,"-",individual,"-1-2copiesMixtureQP.png"),width=1000,height=1000,res=150)
    plot(x1 - x2,Y - x2,ylab="indl2r - X2mean",xlab="X1mean - X2mean",pch=20,main=paste0("w = ",w3$solution[1],"\nw1Fit = ",compareW1,";  w0Fit = ",compareW0),xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),col=add.alpha("black",0.01))
    abline(0,1,lty=3,col="gray")
    abline(0,w3$solution[1],col="red")
    dev.off()
    
    # plot disribution of residuals
    png(paste0("plots/mixtureResults/",filename,"-allOddSamples-PARsorted-sampleIndex-",ind,"-",individual,"-1-2copiesMixtureQP-resid.png"),width=1000,height=1000,res=150)
    hist(residQP,breaks=100)
    abline(v=0)
    dev.off()
}

# apply to all samples with low X (and low Y)
w300 = sapply(c(1:300,700:750,1100:1150), function(ind) {
    print(ind)
    individual = hclustOrderSamples[ind]
    Y = log2Ratios[hmmToPlot$orderOfHMMSNPs[keepSnps],individual]
    x1 = hmmToPlot$parameters[[1]][keepSnps,1,1]
    x2 = hmmToPlot$parameters[[1]][keepSnps,2,1]
    Ws = mixtureTest(Y,x1,x2)
    w = Ws[[1]]; w2 = Ws[[2]]; w3 = Ws[[3]]
    return(w3$solution[1])
})


inds = hclustOrderSamples[c(1:300,700:750,1100:1150)]
y = meanL2rY[inds]
colors = rep("green",length(y))
colors[inds%in%otherInfo$PIID[otherInfo$Inferred.Gender=="F"]] = "green"
colors[inds%in%otherInfo$PIID[otherInfo$Inferred.Gender=="M"]] = "purple"

png(paste0("plots/",filename,"-allOddSamples-PARsorted-1-2copiesMixtureQP-mixedSamplesByYl2r.png"),width=1000,height=1000,res=150)
plot(w300,y,xlab="w estimate",ylab="mean log2ratio Y",col=colors,pch=4)
dev.off()

y = otherInfo$Y.intensity/otherInfo$X.intensity
y = y[match(inds,otherInfo$PIID)]
png(paste0("plots/",filename,"-allOddSamples-PARsorted-1-2copiesMixtureQP-mixedSamplesByYXintensityAffy.png"),width=1000,height=1000,res=150)
plot(w300,y,xlab="w estimate",ylab="Y:X intensity ratio (Affy)",col=colors,pch=4)
dev.off()

h = inds[(y > 0.5)&(colors=="green") & (w300<0.4)]
for(i in quantInfoFields){
    print(i)
    y = otherInfo[[i]][match(inds,otherInfo$PIID)]
    coln = i
    png(paste0("plots/",filename,"-allOddSamples-PARsorted-1-2copiesMixtureQP-mixedSamplesBy-",coln,".png"),width=1000,height=1000,res=150)
    plot(w300,y,xlab="w estimate",ylab=coln,col=colors,pch=4)
    dev.off()
}


h = inds[(y > 0.5)&(colors=="green") & (w300<0.4)]
x = meanL2rY[inds][colors=="green"]
for(i in quantInfoFields){
    print(i)
    y = otherInfo[[i]][match(inds[colors=="green"],otherInfo$PIID)]
    coln = i
    png(paste0("plots/",filename,"-allOddSamples-PARsorted-1-2copiesMixtureQP-Yl2rBy-",coln,"-females.png"),width=1000,height=1000,res=150)
    plot(x,y,xlab="mean log2ratio Y",ylab=coln,col=colors[colors=="green"],pch=4)
    dev.off()
}

coln = "Standing.height"
inds = names(meanL2rY)
y = otherInfo[[coln]][match(inds,otherInfo$PIID)]
x = meanL2rY
colors = rep("green",length(y))
colors[inds%in%otherInfo$PIID[otherInfo$Inferred.Gender=="F"]] = "green"
colors[inds%in%otherInfo$PIID[otherInfo$Inferred.Gender=="M"]] = "purple"
png(paste0("plots/",filename,"-allOddSamples-PARsorted-1-2copiesMixtureQP-Yl2rBy-",coln,"-allSamples.png"),width=1000,height=1000,res=150)
plot(x,y,xlab="mean log2ratio Y",ylab=coln,col=colors,pch=4)
dev.off()






################################################################
# RUN ALL SAMPLES CHROM Y (BY BATCH) USING EM HMM AS STATE EMISSIONS. SAVE RDATA FILE.

filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrY-2734854-to-28767492"
load(paste0(filename,".Rdata"),verbose=TRUE)

log2Ratios = log2RatiosAll[snpChromAll%in%c(24),]  # Y chromosome
chrom="Y"

# get model parameters by using EM hmm from 1300 samples above (is this enough???)
modelParameters = getSNPposteriors(l2r=log2Ratios[hmmResultsRound$orderOfHMMSNPs,],stateProbabilities=hmmResultsRound$HMM[[2]])

# run one iteration of hmm for each batch. Will take some time.
snps = rownames(log2Ratios)

###### Y chromosome can fit in one file

hmmResultsBatches = list()
for( batch in all.batches()){
    print(batch)
    print(date())
    l2r = read.cnv.hdf5(inds=NULL,snps=snps,type="log2ratio",batch=batch)    
    results = runHMM.em(l2r[[1]],initialStateProbabilities=NULL,trans=trans,outIterations=2,tvdTol=1e-4,Expectation=modelParameters,ExpectationOrderOfSNPs=hmmResultsRound$orderOfHMMSNPs)
    resultsRound = results; resultsRound$HMM[[1]] = round(resultsRound$HMM[[1]],3)
    save(hmmResultsBatches,trans,file=paste0(filename,"-byBatch-",batch,".RData"))
    
    hmmResultsBatches[[batch]] = resultsRound
}
save(hmmResultsBatches,trans,file=paste0(filename,"-byBatch.RData"))

                                        # plot each batch separately
for( batch in all.batches()){
    print(batch)
    results = hmmResultsBatches[[batch]]
    plotHMM(results,samplesToPlot=NULL,filename=paste0("plots/",filename,"-",batch),threshold=0.98,plotRecomb=FALSE,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))
}

# sort by mean l2r on Y

# just pick 'males'
keep = (meanL2rX < -0.2) | (meanL2rY > -1)
# plot 5% quantiles
quantile(meanL2rY[keep],seq(0,1,0.05))
breaks = seq(-2,0.6,0.05)
breaks = c(-2,-0.5,-0.3,-0.25)

q1000 = quantile(meanL2rY[keep],seq(0,1,1000/sum(keep)))
breaks1 = c(0,q1000[1:10])
breaks2 = c(q1000[(length(q1000) - 9):length(q1000)],1)

for( breaks in list(breaks1,breaks2) ){
    print(breaks)
    for(b in 2:length(breaks)){
        ranges = breaks[(b-1):b]
        if(ranges[1]==0) next
        inds = names(meanL2rY)[(meanL2rY >= ranges[1]) & (meanL2rY < ranges[2]) & keep]
        print(length(inds))
        indsBatch=otherInfo$Batch[match(inds,otherInfo$PIID)]
        # get data from batch output
        bs = unique(indsBatch)
        hmm = sapply(bs,function(batch) {
            bInds = inds[indsBatch==batch]
            thisHMM = hmmResultsBatches[[batch]]
            out = thisHMM$HMM[[1]][,match(bInds,thisHMM$orderOfHMMSamples),]
        })
        allHmm = abind(hmm,along=2)
        hmmSampleOrder = unlist(sapply(bs,function(batch) {
            bInds = inds[indsBatch==batch]
        },USE.NAMES=FALSE))
        hmmToPlot = list("HMM"=list(allHmm),"orderOfHMMSamples"=hmmSampleOrder,"orderOfHMMSNPs"=hmmResultsBatches[[1]]$orderOfHMMSNPs)
        r = round(ranges,3)
        # sort by meanl2r:
        sOrder = hmmSampleOrder[order(meanL2rY[hmmSampleOrder])]
        plotHMM(hmmToPlot,samplesToPlot=sOrder,filename=paste0("plots/",filename,"-meanL2RY-btwn_",r[1],"_",r[2]),threshold=0.98,plotRecomb=FALSE,sortInds=TRUE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))
    }
}




########################################
# RUN ALL SAMPLES ON CHROM X. SAVE AS HDF5. THIS USES PARAMETER ESTIMATES FROM EXEMPLAR SUBSET.
    
#filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-155233115"
filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2-hmmResults-chrX-60425-to-155233115"
load(paste0(filename,".Rdata"),verbose=TRUE)

log2Ratios = log2RatiosAll[snpChromAll%in%c(23,25),]  # X chromosome.  This is used just for the SNP information!
chrom="X"

snpOrder = rownames(log2Ratios)
snpInfo = ps2snp[match(snpOrder,ps2snp$ProbeSetID),]
snpChrom = snpInfo$Chromosome
snpPos = snpInfo$Position

# get model parameters by using EM hmm from 1300 samples above (is this enough???)
#modelParameters = getSNPposteriors(l2r=log2Ratios[hmmResultsRound$orderOfHMMSNPs,],stateProbabilities=hmmResultsRound$HMM[[2]])
modelParameters = hmmResultsRound$parameters

# get list of snps to use in this version    
snps = hmmResultsRound$orderOfHMMSNPs

# run one iteration of hmm for each batch. Will take some time.
###### Save separately in batches - hdf5 files
outh5 = paste0(filename,"-byBatch.h5")
h5createFile(outh5)

for( batch in all.batches() ){
    
    print(batch)
    print(date())
    l2r = read.cnv.hdf5(inds=NULL,snps=snps,type="log2ratio",batch=batch)    
    results = runHMM.em(l2r[[1]],initialStateProbabilities=NULL,trans=modelParameters[[2]],outIterations=2,tvdTol=1e-4,Expectation=modelParameters[[1]],ExpectationOrderOfSNPs=hmmResultsRound$orderOfHMMSNPs,saveParameters=FALSE,computeTransitions=FALSE)
    resultsRound = results; resultsRound$HMM[[1]] = round(resultsRound$HMM[[1]],3)

    h5write(resultsRound$HMM[[1]],outh5,paste0(batch,".HMM"))
    h5write(resultsRound$orderOfHMMSNPs,outh5,paste0(batch,".orderOfHMMSNPs"))
    h5write(resultsRound$orderOfHMMSamples,outh5,paste0(batch,".orderOfHMMSamples"))
    # no need to store convergence stats, as these are easily derivable from hmm output.
}

h5write(modelParameters,outh5,"parameters")        
H5close()


################## Plot chrom X by chunks and order by cnv values
h5 = paste0(filename,"-byBatch.h5")
load(paste0(filename,"-hclustOrderPARoddSamples.RData")) # get order of odd samples


chunkWidth=10  # with of chunks in MB
for( start in seq(0,round(max(snpPos)/(10^6)),chunkWidth)){

    #start = 86
    print(start)
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
    dataHMM = read.hmm.hdf5(h5,inds=NULL,snps=theseSnps,otherInfo=otherInfo)
    W = transformHMMtoCalls(dataHMM,threshold=0.98)
    nNot0 = colSums(W!=0)
    frac1 = colSums(W==1)/nNot0
    frac2 = colSums(W==2)/nNot0
    frac3 = colSums(W==3)/nNot0
    frac0 = colSums(W==0)/nrow(W)
    keepGood = colnames(W) %in% otherInfo$PIID # exclude references
    keep = (frac3 > 0.1) | ( (frac1 > 0.1 )&(frac1 < 0.9) ) | ( (frac2 > 0.1 )&(frac2 < 0.9) )
    keepInds = colnames(W)[keep & keepGood & ( frac0 < 0.5 )]
    keepW = W[,keepInds]
    samplesToPlot = hclust.order(keepW)
    samplesToPlot = c(hclustOrderSamples,samplesToPlot[!samplesToPlot%in%hclustOrderSamples])
    samplesToPlotKEEP = samplesToPlot
    
    plotHMM(dataHMM,samplesToPlot=samplesToPlotKEEP,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))
    
}


# just read in samples with odd values at 86-96 at start and end of X chromosome. samplesToPlot is from above.
theseSnps = snpOrder[( (snpPos <= (5)*10^6) ) | ( snpPos > (par2-1e6) )]
dataHMM86 = read.hmm.hdf5(h5,inds=samplesToPlot,snps=theseSnps,otherInfo=otherInfo)
for( start in c(0,150)){    
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
    plotHMM(dataHMM86,samplesToPlot=samplesToPlot,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-allSamples-meanSorted-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    
}

### get order for 'odd samples' based on PAR sorting
theseSnps = snpOrder[( (snpPos <= (5)*10^6) ) | ( snpPos > (par2-1e6) )]
oddSamplesHMM = read.hmm.hdf5(h5,inds=oddSamples,snps=theseSnps,otherInfo=otherInfo)
W = transformHMMtoCalls(oddSamplesHMM,threshold=0.98)
Wna = W; Wna[Wna==0] = NA
d = dist(t(Wna)) # euclidian distance, ignoring uncertain values coded as NAs. takes a few seconds on 1500 x 1500
j = hclust(d)
hclustOrderSamples = colnames(W)[j$order]

# save this!
#save(hclustOrderSamples,file=paste0(filename,"-hclustOrderPARoddSamples.RData"))



######
#### Look for 'odd samples' based on PAR1 boundary as well.
theseSnps = snpOrder[( (snpPos >= (par1 - 0.25*1e6) ) & ( snpPos <= (par1+0.25*1e6) ) )]
oddSamplesHMM = read.hmm.hdf5(h5,inds=NULL,snps=theseSnps,otherInfo=otherInfo)
W = transformHMMtoCalls(oddSamplesHMM,threshold=0.98)
Wna = W; Wna[Wna==0] = NA
# count fractions of 0/1/2/3 inside and outside PAR
inPAR = rownames(W)%in%snpOrder[snpPos<par1]
nNot0 = colSums(W[inPAR,]!=0)
frac1 = colSums(W[inPAR,]==1)/nNot0
frac2 = colSums(W[inPAR,]==2)/nNot0
frac3 = colSums(W[inPAR,]==3)/nNot0
frac0 = colSums(W[inPAR,]==0)/nrow(W[inPAR,])

nNot0b = colSums(W[!inPAR,]!=0)
frac1b = colSums(W[!inPAR,]==1)/nNot0b
frac2b = colSums(W[!inPAR,]==2)/nNot0b
frac3b = colSums(W[!inPAR,]==3)/nNot0b
frac0b = colSums(W[!inPAR,]==0)/nrow(W[!inPAR,])

X = cbind(frac1,frac2,frac1b,frac2b)
X[is.na(X)]=0
pcaFracs = princomp(X)
colours = rep(NA,nrow(X)) # <=== Only plotting the samples in otherInfo
males = rownames(X)%in%otherInfo$PIID[otherInfo$Inferred.Gender=="M"]
females = rownames(X)%in%otherInfo$PIID[otherInfo$Inferred.Gender=="F"]
colours[males] = sexCols["M"]
colours[females] = sexCols["F"]

png(paste0("plots/",filename,"-PAR1-boundary-fraction1or2.png"),bg='transparent')
#plot(pcaFracs$scores[,1],pcaFracs$scores[,2],col=colours,pch=16)
#plot(pcaFracs$scores[,3],pcaFracs$scores[,4],col=colours,pch=16)
#plot(frac2,frac2b,col=colours,pch=16)
hist(frac2b[females],breaks=100,ylim=c(0,1000),xlim=c(0,1),col=add.alpha(sexCols["F"],0.5),border=NA,xlab="Fraction 2 copies outside PAR1 boundary")
hist(frac2b[males],breaks=100,col=add.alpha(sexCols["M"],0.5),border=NA,add=TRUE)
dev.off()

parExtendSamplesList = colnames(W)[(frac2b>0.4)&(frac2b<0.65)&(colnames(W)%in%otherInfo$PIID)]
d = dist(t(Wna[,parExtendSamplesList])) # euclidian distance, ignoring uncertain values coded as NAs. takes a few seconds on 1500 x 1500
j = hclust(d)
hclustOrderParExtendSamples = colnames(Wna[,parExtendSamplesList])[j$order]

# save this!
save(hclustOrderParExtendSamples,file=paste0(filename,"-PAR1ExtendedSamples.RData"))




########################################
# PLOT the l2r and inferred path for some strange individuals

# read in l2r for these individuals
log2RatiosKEEP = read.cnv.hdf5(inds=samplesToPlotKEEP,snps=snpOrder,type="log2ratio",otherInfo=otherInfo)
log2RatiosKEEP = abind(log2RatiosKEEP,along=2)  

for( start in c(0,150)){    
    theseSnps = snpOrder[(snpPos > start*10^6) & (snpPos <= (start+chunkWidth)*10^6)]
    plotL2RIndividual(log2RatiosKEEP,dataHMM86,individual,snpsToPlot=theseSnps,filename=paste0("plots/",filename,"-",start,"-",(start+chunkWidth),"MB"),threshold=0.98,recombRates=recombRatesX,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))    
}
    

########################################
# CHECK CONSISTENCY OF PARAMETER ESTIMATES.

# NOTE: This version of the hmm uses different means and variances for each SNP. The question is: is this necessary? Learning much fewer parameters would be helpful, but there is so much data that it might not be necessary to have such a restriction if there is sufficient variation in log2Ratios across SNPs.

# Read in results from exemplar set across whole X chromosome
#filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-155233115"
filename="b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-v2-hmmResults-chrX-60425-to-155233115"
load(paste0(filename,".Rdata"),verbose=TRUE)





########################################
# RUN SAMPLES WITH LOW Y ON FIRST CHUNK OF X CHROMOSOME, USING EM HMM AS STATE EMISSIONS.
# ACTUALLY THIS IS REDUNDANT AFTER RUNNING ALL SAMPLES ON X
    
# get list of low Y samples to do HMM on chunk of X chromosome
###############
ranges = breaks1[2:3]
inds = names(meanL2rY)[(meanL2rY >= ranges[1]) & (meanL2rY < ranges[2]) & keep]
log2Ratios = log2RatiosAll[snpChromAll%in%c(23,25) & (snpPosAll < 5*(10^6)),]  # X chromosome near PAR1 ===> samples with low Y meanYl2r
chrom="X"
############### 

sampleOrder = colnames(log2Ratios)
snpOrder = rownames(log2Ratios)
snpInfo = ps2snp[match(snpOrder,ps2snp$ProbeSetID),]
snpChrom = snpInfo$Chromosome
snpPos = snpInfo$Position

load("b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-155233115.Rdata",verbose=TRUE)

# get modelParameters for these snps
theseSNPs = hmmResultsRound$orderOfHMMSNPs%in%snpOrder
modelParameters = getSNPposteriors(l2r=log2Ratios[hmmResultsRound$orderOfHMMSNPs[theseSNPs],],stateProbabilities=hmmResultsRound$HMM[[2]][theseSNPs,,])

l2r = read.cnv.hdf5(inds=inds,snps=snpOrder,type="log2ratio",otherInfo=otherInfo)
l2r2 = abind(l2r,along=2)

results = runHMM.em(l2r2,initialStateProbabilities=NULL,trans=trans,outIterations=2,tvdTol=1e-4,Expectation=modelParameters,ExpectationOrderOfSNPs=hmmResultsRound$orderOfHMMSNPs[theseSNPs])

filename = paste0(outFile,"-hmmResults-chr",chrom,"-",min(snpPos),"-to-",max(snpPos),"-lowYinds")

hmmResultsRound = results; hmmResultsRound$HMM[[1]] = round(hmmResultsRound$HMM[[1]],3)
save(hmmResultsRound,trans,file=paste0(filename,".Rdata"))


#### Plot alongside X chromosome results for these individuals
fnameY = "b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrY-2734854-to-28767492-byBatch"
load(paste0(fnameY,".RData"),verbose=TRUE)
fnameX = "b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-hmmResults-chrX-60425-to-4993450-lowYinds"
load(paste0(fnameX,".Rdata"),verbose=TRUE)
toPlotX = hmmResultsRound


# get data from batch output
indsBatch=otherInfo$Batch[match(inds,otherInfo$PIID)]
bs = unique(indsBatch)
hmm = sapply(bs,function(batch) {
    bInds = inds[indsBatch==batch]
    thisHMM = hmmResultsBatches[[batch]]
    out = thisHMM$HMM[[1]][,match(bInds,thisHMM$orderOfHMMSamples),]
})
allHmm = abind(hmm,along=2)
hmmSampleOrder = unlist(sapply(bs,function(batch) {
    bInds = inds[indsBatch==batch]
},USE.NAMES=FALSE))
hmmToPlot = list("HMM"=list(allHmm),"orderOfHMMSamples"=hmmSampleOrder,"orderOfHMMSNPs"=hmmResultsBatches[[1]]$orderOfHMMSNPs)

# exclude (from plots) samples with lots of red in non-PAR X (likely females)
snpPos = ps2snp$Position[match(toPlotX$orderOfHMMSNPs,ps2snp$ProbeSetID)]
xData = toPlotX$HMM[[1]][snpPos>par1,,]
fracRed = colSums(xData[,,2] >=0.98)/colSums(!is.na(xData[,,2]))
g=hist(fracRed,plot=FALSE)
indsM = inds[otherInfo$Submitted.Gender[match(inds,otherInfo$PIID)]=="M"]
inds1X = toPlotX$orderOfHMMSamples[fracRed<0.1]
indsMnomismatch = indsM[otherInfo$Inferred.Gender[match(indsM,otherInfo$PIID)]=="M"]
r = round(max(meanL2rY[inds]),3)

# get order of samples - male submitted
Order = order(meanL2rY[indsM]) # sort by mean l2r
# plot the Y chromosome
plotHMM(hmmToPlot,samplesToPlot=indsM[Order],filename=paste0("plots/",fnameY,"-lowYinds-MaleSubmitted-lt_",r),threshold=0.98,plotRecomb=FALSE,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain","1"="0","2"="1","3"="2"))

# plot the X chromosome in the same order
plotHMM(toPlotX,samplesToPlot=indsM[Order],filename=paste0("plots/",fnameX,"-lowYinds-MaleSubmitted-lt_",r),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))

# get order of samples - 1X
Order = order(meanL2rY[inds1X]) # sort by mean l2r
# plot the Y chromosome
plotHMM(hmmToPlot,samplesToPlot=inds1X[Order],filename=paste0("plots/",fnameY,"-lowYinds-1X-lt_",r),threshold=0.98,plotRecomb=FALSE,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain","1"="0","2"="1","3"="2"))

# plot the X chromosome in the same order
plotHMM(toPlotX,samplesToPlot=inds1X[Order],filename=paste0("plots/",fnameX,"-lowYinds-1X-lt_",r),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))


# get order of samples - male submitted and inferred
Order = order(meanL2rY[indsMnomismatch]) # sort by mean l2r
# plot the Y chromosome
plotHMM(hmmToPlot,samplesToPlot=indsMnomismatch[Order],filename=paste0("plots/",fnameY,"-lowYinds-MaleSubmittedAndInferred-lt_",r),threshold=0.98,plotRecomb=FALSE,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain","1"="0","2"="1","3"="2"))

# plot the X chromosome in the same order
plotHMM(toPlotX,samplesToPlot=indsMnomismatch[Order],filename=paste0("plots/",fnameX,"-lowYinds-MaleSubmittedAndInferred-lt_",r),threshold=0.98,recombRates=recombRatesX,sortInds=FALSE,extraPlotting=FALSE,updateLegend=list("0"="uncertain"))
