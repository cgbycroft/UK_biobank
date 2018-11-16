## Script to summarise cnv data for sex chromosomes.
## Output is for each sample, the mean, median, and sd of baf and l2r for each part of the sex chromosomes - PAR1, PAR2, X, Y, MT separately
#NOTE: R library rhdf5 must be installed for this to work. This is a part of bioConductor. To install:
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-out","b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries")


print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}

outFile = args[which(args=="-out")+1]

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

outData=list()
for(type in c("log2ratio","baf")){
    
    #compute summary statistics by chromosome. Also only read in samples in fam file (i.e excludes references, but not duplicates(?)).
    # out2 file is nSamples x 3*nChromosomes
    # takes about 10secs per batch, if reading X chromosome. Faster for others. Probably takes about 1.5hrs to do the whole thing with both measures.
    print(type)
    
    out = sapply(chromosomes,function(chr){
        print(chr)        
        snps = ps2snp$ProbeSetID[(ps2snp$AffySNPID%in%snpsToInclude) & (ps2snp$Chromosome2==chr) ]
        Summaries = sapply(batches,function(batch){
            dataSubset = read.cnv.hdf5(inds=fam,snps=snps,type=type,otherInfo=otherInfo,batch=batch)
            means = colMeans(dataSubset[[1]],na.rm=TRUE)
            sds = apply(dataSubset[[1]],2,sd,na.rm=TRUE)
            miss = colSums( is.na(dataSubset[[1]]) )
            summaries = cbind(means,sds,miss)
            colnames(summaries) = paste0(c("means","sds","miss"),".",chr)
            return(summaries)
        })
        Summaries = Reduce(rbind,Summaries)
        
        return(Summaries)
        
    },simplify=FALSE)

    out2 = Reduce(cbind,out)
    outData[[type]] = out2
    
}


# write to file
save(outData,file=paste0(baseSampleQCDir,"/data/CNV/",outFile,".RData"))

#[c(52235,52397)]


## compute over samples (i.e a value for each SNP)
outData2=list()
for(type in c("log2ratio","baf")){
    
    #compute summary statistics by chromosome. Also only read in samples in fam file (i.e excludes references, but not duplicates(?)).
    # out2 file is nSamples x 3*nChromosomes
    # takes about 10secs per batch, if reading X chromosome. Faster for others. Probably takes about 1.5hrs to do the whole thing with both measures.
    print(type)
    
    snps = ps2snp$ProbeSetID[(ps2snp$AffySNPID%in%snpsToInclude)]
    out = sapply(batches,function(batch){
        dataSubset = read.cnv.hdf5(inds=fam,snps=snps,type=type,otherInfo=otherInfo,batch=batch)
        means = rowMeans(dataSubset[[1]],na.rm=TRUE)
        sds = apply(dataSubset[[1]],1,sd,na.rm=TRUE)
        miss = rowSums( is.na(dataSubset[[1]]) )
        summaries = cbind(means,sds,miss)
        colnames(summaries) = paste0(c("means","sds","miss"))
        return(summaries)
    })
    
    out2 = Reduce(rbind,out)   
    outData2[[type]] = out2    
}


# write to file
save(outData2,file=paste0(baseSampleQCDir,"/data/CNV/",outFile,"-bySNP.RData"))
