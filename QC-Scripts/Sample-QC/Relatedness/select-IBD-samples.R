#####
# Script to select set of samples to run plink IBD
#####


args = commandArgs(trailingOnly=TRUE)
print(args)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-wb","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt","-rel","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0")

print(args)
h = args[-c(which(args%in%c("-wb","-rel")),1+which(args%in%c("-wb","-rel")))]
for(helperScript in h){
    source(helperScript)
}

kingFile = args[which(args=="-rel")+1]

kin = read.table(kingFile,stringsAsFactors=FALSE,header=TRUE)
wb = read.table(args[which(args=="-wb")+1],stringsAsFactors=FALSE,header=FALSE)[,1]


                                        # select white british samples in kinship table, where both of the pair are white british

toKeepKin = kin[(kin$ID1%in%wb)&(kin$ID2%in%wb),]

print("N pairs to keep, i.e where both of the pair are white british.")
print(dim(toKeepKin)[1])

toKeep = unique(c(toKeepKin$ID1,toKeepKin$ID2))
print("Makes N individuals to keep:")
print(length(toKeep))

                                        # write a kinship table with just the white british filtered
write.table(toKeepKin,quote=FALSE,col.names=TRUE,row.names=FALSE,file=gsub(".kin0","-White_British.kin0",kingFile))

# write the list of samples to keep for ibd analysis
write.table(cbind(toKeep,toKeep),quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/plink-IBD-samples.txt"))

print("List for plink saved in:")
print(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/plink-IBD-samples.txt"))


print("WhiteBritish only King kinship file save as:")
print( gsub(".kin0","-White_British.kin0",kingFile) )

