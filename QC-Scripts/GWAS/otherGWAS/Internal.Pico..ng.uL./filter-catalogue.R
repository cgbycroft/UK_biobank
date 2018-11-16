h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

library(dplyr)
library(qqman)
library(stringr)


################################
phenotype = "Internal.Pico..ng.uL."
################################

outFile = paste0(baseSampleQCDir,"/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-",phenotype)


######## Get GWAS catalogue information

catFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog.txt.gz"

print(paste0("Using catalogue file: ",catFile))
catRaw = read.delim(gzfile(catFile),header=FALSE,stringsAsFactors=FALSE,sep="\t")

# select the appropriate variants associated with this phenotype
queryString = "White blood|cell|Lymphocyte|Neutro"
hits = grep(queryString,catRaw$V11,ignore.case=TRUE)
catPheno = catRaw[hits,]
print( table(catPheno$V11) )

# Then manually check these hits - are they appropriate??
catPheno = catPheno[!grepl("cancer|carcinoma|lymphoma|tumor|chemotherapy|pressure|anemia|F-cell|beta-cell|leukemia",catPheno$V11,ignore.case=TRUE),]
print( table(catPheno$V11) )

print( paste0("Apparently ",dim(catPheno)[1]," raw hits.") )

# Check for European ancestry
catPheno$Ancestry = gsub("[[:digit:]]|,|individuals","",catPheno$V12,ignore.case=TRUE)
catPheno$Ancestry = str_trim(catPheno$Ancestry)

print( "Raw ancestry information: ") 
print( table(catPheno$Ancestry) )

print( "After subsetting for European: ") 
Euro = grep("European",catPheno$Ancestry)
print( table(catPheno$Ancestry[Euro]) )

catPhenoEur = catPheno[Euro,]
print( paste0("Apparently ",dim(catPhenoEur)[1]," hits for Europeans only.") )


# any other filters??

# write output files
# European ancestry
write.table(catPhenoEur,file=paste0(outFile,"-Eur.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

# Any ancestry
write.table(catPheno,file=paste0(outFile,".txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

save(catPheno,catPhenoEur,file=paste0(outFile,".RData"))

print( "files written to: ")
print( paste0(outFile,"*") )
