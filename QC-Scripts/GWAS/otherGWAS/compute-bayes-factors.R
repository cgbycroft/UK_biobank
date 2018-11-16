#####################
# Script to computing bayes factors for a given GWAS at a given set of regions.
#####################

#####################
# Preliminaries

args = commandArgs(TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

                                        # args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out.signif","-regions","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts","-chr","2")

# args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz","-prior","0.2")

#args = c("-gwas","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v20/Standing.height-BOLT-LMM-v20-chr%%.out","-regions","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts","-chr","5","-prior","0.2")



print(args)

dataFile = args[which(args=="-gwas")+1]


priors = as.numeric(str_split(args[which(args=="-prior")+1],",")[[1]])
    

if(grepl("%%",dataFile)) {
    
    chrom = as.numeric(args[which(args=="-chr")+1])
    dataFile = gsub("%%",chrom,dataFile)
} else {
    chrom="genome"
}

outFile = paste0(dataFile,".",priors[1],".bf")

print(outFile)
#quit()

## Read in the GWAS data
if(grepl("GIANT",dataFile)) {
    theseColumns = get.colClasses(file=dataFile,c("b","SE"),colTypes=c("numeric","numeric"))
    theseColumns = c(NA,theseColumns) # R reads in a rownames column too
} else {
    theseColumns = get.colClasses(file=dataFile,c("BETA","SE"),colTypes=c("numeric","numeric"))
}

if(grepl(".gz$",dataFile)) DF = tryCatch(read.table(gzfile(dataFile),sep="",header=TRUE,stringsAsFactors=FALSE,colClasses=theseColumns), error=function(e) NULL) else DF = tryCatch(read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE,colClasses=theseColumns), error=function(e) NULL)


if(is.null(DF)){
    print("No GWAS data available for this chromosome.")
    quit()
}

colnames(DF)[colnames(DF)=="b"] = "BETA"
print(head(DF))

print(paste0("Computing bayes factors for ",nrow(DF)," markers on chromosome: ",chrom))


#######
# Compute the Bayes Factors

#test = calculate.BF(DFsubset$BETA,DFsubset$SE,g=0.2)

for(prior in priors){
    print(prior)
    BF = calculate.BF(DF$BETA,DF$SE,g=prior)
    colnames(BF) = paste0(colnames(BF),".",prior)
    if(prior != priors[1]) outBF = cbind(outBF,BF) else outBF=BF
}


####### print the data

print("Writing out to file (in order of original GWAS out data).")
write.table(outBF,file=outFile,quote=FALSE,col.names=TRUE,row.names=FALSE)


print("DONE! see:")
print(outFile)
