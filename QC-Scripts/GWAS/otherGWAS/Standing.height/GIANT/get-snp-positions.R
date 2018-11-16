########
# Script to add snp positions to the giant file
########



giantFile ='/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz'    


giant = read.table(gzfile(giantFile),stringsAsFactors=FALSE,header=TRUE)

# get positions from ucsc
snpDBFile = '/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/ucsc.snps.annot/snp147.subsetFields.txt'
snpDF = read.table(snpDBFile,stringsAsFactors=FALSE,header=FALSE,colClasses=c("character","numeric","numeric","character","character","character"))


# include position in giant file
giantIndex = match(giant$MarkerName,snpDF$V4)
giant$BP = snpDF$V3[giantIndex] # end position is the right one!
colnames(giant)[c(1,7)] = c("SNP","P2")
giant$P = giant$P2
giant$CHRraw = snpDF$V1[giantIndex]
giant$CHR = as.numeric(gsub("chr","",giant$CHRraw))
giant$CHR[is.na(giant$CHR)] = 0
giant$MAF = as.numeric(giant$Freq.Allele1.HapMapCEU)
giant$MAF[(giant$MAF>0.5)&(!is.na(giant$MAF))] = 1-giant$MAF[(giant$MAF>0.5)&(!is.na(giant$MAF))]
# exclude snps not found in uscs data (we can't plot them anyway)
#giant2 = giant[!is.na(giantIndex),]

#giantData = list("DF"=giant2,"Pvalset"="GIANT.chrgenome")


outFile = '/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz'

gz1 = gzfile(outFile,"w")
write.table(giant, gz1);
close(gz1)
