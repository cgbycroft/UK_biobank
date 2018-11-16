args = commandArgs(trailingOnly=TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

# snp annotations in new release
BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
BB.ps2snpSex = ukbiobank.ps2snp("sexchrom")
BL.ps2snpSex = ukbileve.ps2snp("sexchrom")
ps2snpAll=rbind(BB.ps2snp,BL.ps2snp,BB.ps2snpSex,BL.ps2snpSex)
ps2snpAll=unique(ps2snpAll)

# load the final sample table
load(paste0(baseSampleQCDir,"/data/ForRelease/b1__b11-b001__b095-sampleTable_v4_allColumns.RData"),verbose=TRUE) 

# get list of R files
rFiles = list.files(path=paste0(baseSampleQCDir,"/data/Combined/interimReleaseSamples/Comparison"),pattern="compareToNew.RData")

print(paste0("Plotting ",length(rFiles)," files."))

# loaded in batches x chromosomes (takes a few minutes to load)
CS = sapply(rFiles,function(x){
    load(paste0(baseSampleQCDir,"/data/Combined/interimReleaseSamples/Comparison/",x),verbose=TRUE)
    return(cs)
},simplify=FALSE)

nSNPsSamples = sapply(rFiles,function(x){
    load(paste0(baseSampleQCDir,"/data/Combined/interimReleaseSamples/Comparison/",x),verbose=TRUE)
    return(c(length(intersectSnps),length(intersectSamples)))
},simplify=FALSE)


# sum counts over chromosomes
batchPrefixes = unique(sapply(str_split(rFiles,"\\.chr"),function(x) x[1]))
batches = all.batches()[1:33]

CSsum = sapply(batches,function(b){
    bs = grep(paste0("*",b,".chr"),names(CS),value=TRUE)
    #print(b)
    #print(bs[1])
    if(length(bs)<26) print(paste0(b," only has ",length(bs)," chromosomes done."))
    allChrom = Reduce('+', CS[bs])  # sum over all the chromosomes
    return(allChrom)
},simplify=FALSE,USE.NAMES=FALSE)


nSNPsSamplessum = sapply(batches,function(b){
    bs = grep(paste0("*",b,".chr"),names(CS),value=TRUE)
    #print(b)
    #print(bs[1])
    if(length(bs)<26) print(paste0(b," only has ",length(bs)," chromosomes done."))
    allChrom = sum(sapply(nSNPsSamples[bs],function(x) x[1]))  # sum over all the chromosomes
    allSamples = nSNPsSamples[bs][[1]][2]
    if(sum(sapply(nSNPsSamples[bs],function(x) x[2])!=allSamples)>0) print("WARNGING!")
     # all chroms should have the same value

    return(c(allChrom,allSamples))
},simplify=TRUE,USE.NAMES=TRUE)



# compute various stats per sample
fracNonMissing = unlist( lapply(CSsum,function(x) x[,"nCallChanges"]/x[,"nonMissingBoth"]) )
missToNonmiss = unlist( lapply(CSsum,function(x) x[,"mToNm"]/rowSums(x[,c("mToNm","NmTom")])))
nonmissToMiss = unlist( lapply(CSsum,function(x) x[,"NmTom"]/rowSums(x[,c("mToNm","NmTom")])))
nonMissingBoth = unlist( lapply(CSsum,function(x) x[,"nonMissingBoth"]) )

MissingInOne = unlist( lapply(CSsum,function(x) rowSums(x[,c("mToNm","NmTom")])) )


MissingRateInterim = unlist( lapply(CSsum,function(x) rowSums(x[,c("mToNm","NmTom")])) )
MissingRate = unlist( lapply(CSsum,function(x) rowSums(x[,c("mToNm","NmTom")])) )


nSamples = sum(!is.na(fracNonMissing))
nSNPS = nSNPsSamplessum[1,1]
    
print( paste0( nSamples, " samples in overlap.") )

outTableIndex = match(names(fracNonMissing),outTable$PIID)
thisBatch = outTable$Batch[outTableIndex]
thisBatch = factor(thisBatch,levels=batches)


#### Save output
save(CSsum,nSNPsSamplessum,thisBatch,file="Interim-final-comparison-Totals.RData")



####################################
# Investigate low values in nonmissToMiss

toCheck = names(nonmissToMiss)[(nonmissToMiss>0.4)&(!is.na(nonmissToMiss))]

ch = outTable[outTable$PIID%in%toCheck,]

imageArte = read.table('/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-proj40-dim50-pc17-19-samplesInCluster.txt',stringsAsFactors=FALSE)[,1]

sum(!toCheck%in%imageArte) # <===== ALL ARE THE IMAGE ARTEFACT INDIVIDUALS!
    
#write.table(cbind(toCheck,toCheck),file="odd-samples-to-check.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)





####################################
# Plot histogram of various measures

png(paste0("plots/InterimComparison-histograms-%02d.png"),width=1000,height=1000,res=150 )

hist(fracNonMissing,breaks=100,xlab="Fraction of non-missing markers discordant",main=paste0(nSamples," overlapping samples.\n",nSNPS," overlapping markers" ))

hist(nonmissToMiss,breaks=100,xlab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples.\n",nSNPS," overlapping markers" ))

hist(missToNonmiss,breaks=100,xlab="Fraction of one missing calls in interim data",main=paste0(nSamples," overlapping samples.\n",nSNPS," overlapping markers" ))

hist(log10(fracNonMissing),breaks=100,xlab="Fraction of non-missing markers discordant (log10)",main=paste0(nSamples," overlapping samples.\n",nSNPS," overlapping markers" ))
#axis(1,at=log10(c(1:7)*(10^-4)),labels=as.expression(paste0(c(1:7),"e-4")) )

hist(log10(1-fracNonMissing),breaks=100,xlab="Fraction of non-missing markers concordant (log10)",main=paste0(nSamples," overlapping samples.\n",nSNPS," overlapping markers" ))

dev.off()


####################################
# Plot boxplot by batch

png(paste0("plots/InterimComparison-batch-boxplots-%02d.png"),width=1000,height=1000,res=150 )

boxplot(fracNonMissing~thisBatch,ylab="Fraction of non-missing markers discordant",main=paste0(nSamples," overlapping samples." ),las=2)
boxplot(nonmissToMiss~thisBatch,ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),las=2)
boxplot(missToNonmiss~thisBatch,ylab="Fraction of one missing calls in interim data",main=paste0(nSamples," overlapping samples." ),las=2)

dev.off()



####################################
# Plot scatter

png(paste0("plots/InterimComparison-scatter-%02d.png"),width=1000,height=1000,res=150 )

# colour by het-miss outliers
color=rep("black",length(fracNonMissing))
color[outTable$het.missing.outliers[outTableIndex]==1] = "red"

plot(fracNonMissing,nonmissToMiss,xlab="Fraction of non-missing markers discordant",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(nonMissingBoth,nonmissToMiss,xlab="Number of non-missing markers in both",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(MissingInOne,nonmissToMiss,xlab="Number of markers missing in one",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(MissingInOne,nonMissingBoth,xlab="Number of markers missing in one",ylab="Number of non-missing markers in both",main=paste0(nSamples," overlapping samples." ),col=color)


plot(outTable$dQC[outTableIndex],nonmissToMiss,xlab="dish QC",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(outTable$QC.CR[outTableIndex],nonmissToMiss,xlab="QC.CR",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(outTable$Cluster.CR[outTableIndex],nonmissToMiss,xlab="Cluster.CR",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

plot(outTable$sample.qc.missing.rate[outTableIndex],nonmissToMiss,xlab="Sample QC missing rate",ylab="Fraction of one missing calls in final data",main=paste0(nSamples," overlapping samples." ),col=color)

dev.off()

