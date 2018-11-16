##########
# Script to define a set of regions for comparing results across studies, using just the length of the genome (and no GWAS results).
##########
# Output is a regions file: chrom start end

#####################
# Preliminaries

args = commandArgs(TRUE)

h=c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


# args = c("-ws","1")

##################################################
# Basic: just define without using GWAS results.

windowSize = as.numeric(args[which(args=="-ws")+1])  # in units of Mb

chroms=1:22

# downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes on 10/8/2018
chromosomeLengths = read.table(chromLengthsFile)


# Divide each chromosome into regions

outFile = paste0(baseSampleQCDir,"/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-",windowSize,"Mb.txt")

write.table(t(c("chr","start","end")),file=outFile,quote=FALSE,row.names=FALSE,col.names=FALSE)

nRegions = 0
for( i in chroms){
    print(paste0("chr:",i))
    
    chr = paste0("chr",i)
    l = chromosomeLengths[chromosomeLengths[,1]==chr,2]
    starts = seq(1,l,by=windowSize*1e6)
    ends = starts + windowSize*1e6 - 1
    pos = cbind(i,starts,ends)

    test = (max(starts)<l) & (max(ends) > 1) # check that last bin contains the end.
    if(!test) print("WARNING: end region does not contain end of chromosome!")
    posF = format(pos,scientific=FALSE,trim=TRUE)

    check = sum(as.numeric(posF)!=pos)==0

    if(!check) print("WARNING: something wrong with integer formatting!")
    
    write.table(posF,file=outFile,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)

    nRegions = nrow(pos) + nRegions
    print(nrow(pos))   
}


print("DONE! See:")
print(outFile)
print(paste0("Should be ",nRegions," regions in total, of size ",windowSize," Mb."))
