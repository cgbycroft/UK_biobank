######
# R script to write a list of snp ids in the format provided in the release data, from a list in AFFY ID format (e.g from sample QC).
#####
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
print(args)

#args = c('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.bim','/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.bim')

AffySNPs = read.table(args[1],stringsAsFactors=FALSE,header=FALSE)
genoSNPs = read.table(args[2],stringsAsFactors=FALSE,header=FALSE)

a = tbl_df(AffySNPs); a$chromPos = paste0(a$V1,"_",a$V4)
g = tbl_df(genoSNPs); g$chromPos = paste0(g$V1,"_",g$V4)

check = length(unique(a$chromPos)) == dim(a)[1] # check that snps are all unique on chromosome and position
print("arg1 data is unique on chromosome and position")
print(check)

check = length(unique(g$chromPos)) == dim(g)[1] # check that snps are all unique on
print("arg2 data is unique on chromosome and position")
print(check)

joined = left_join(a,g,by=c("chromPos"="chromPos"))
dim(joined)[1]==dim(a)[1]

                                        # are there any snps in args1 which aren't found in args2?
print("number of snps not found in second file:")
print( sum(is.na(joined$V2.y)) )

write.table(joined$V2.y[!is.na(joined$V2.y)],file=args[3],quote=FALSE,col.names=FALSE,row.names=FALSE)
