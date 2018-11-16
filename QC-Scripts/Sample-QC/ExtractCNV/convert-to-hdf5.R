# script to read in and convert cnv data to hdf5 format
#NOTE: R library rhdf5 must be installed for this to work. This is a part of bioConductor. To install:
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")


args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

for(helperScript in args){
    source(helperScript)
}

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/ExtractCNV"))
batches= all.batches()

# get snp annotations (position etc)
UKBB = ukbiobank.ps2snp("sexchrom")
UKBL = ukbileve.ps2snp("sexchrom")
UKBL = tbl_df(UKBL)
UKBB = tbl_df(UKBB)

ps2snp = merge(UKBL, UKBB, all = TRUE)
ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]


# extract, combine and convert. Takes about 2mins per batch = 
#batch ="UKBiLEVEAX_b2"
for(type in c("log2ratio","baf")){
    
    outh5 = paste0(baseSampleQCDir,"/data/CNV/b1__b11-b001__b095.V2_All.",type,".sexchrom.h5")
    h5createFile(outh5)
    
    for(batch in batches){
        
        print(batch)
        print(date())
        
        data2 = extract.cnv(inds=NULL,batch,type,snpsToInclude=NULL,ps2snp=ps2snp)    
        snpOrder = data2[,1:3]
        sampleOrder = colnames(data2)[-(1:3)]
        data = as.matrix(data2[,-(1:3)])
        # data is a matrix with dimensions nSNPs x nSamples
        h5write(data, outh5,paste0(batch,".data"))
        h5write(snpOrder, outh5,paste0(batch,".snpOrder"))
        h5write(sampleOrder, outh5,paste0(batch,".sampleOrder"))
    }
    H5close()
}


print("DONE CONVERSION.")
print(outh5)


## testing
#otherInfo = read.multiple.batch.info(batchInfoFields)

#indsToRead = otherInfo$PIID[otherInfo$Y.intensity<100]
#snpsToRead = NULL

#dataSubset = read.cnv.hdf5(indsToRead,snpsToRead,type=type,otherInfo=otherInfo)
