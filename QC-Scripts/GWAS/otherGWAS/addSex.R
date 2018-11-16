#### Script to add sex column to bgen *sample file

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

args = commandArgs(trailingOnly=TRUE)
    
#args = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.sample"

input = args[1]
out = paste0(input,".sex")
    
samples = read.table(input,header=TRUE,skip=1,stringsAsFactors=FALSE)

sex = read.multiple.batch.info("Inferred.Gender")

if((sum(duplicated(sex$PIID))>0)|(sum(duplicated(samples[,2]))>0)) print("ERROR: input file id columns are not unique!")

matchSex = sex$Inferred.Gender[match(samples[,2],sex$PIID)]

samples$D=NA
samples$D[matchSex=="F"] = 2
samples$D[matchSex=="M"] = 1

print(table(samples$D,useNA="ifany"))
print(table(matchSex,useNA="ifany"))

myHeader=cbind(read.table(input,header=FALSE,nrow=2,stringsAsFactors=FALSE),c("sex","D"))

write.table(myHeader,file=out,quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(samples,file=out,quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE)

print("DONE!")
print(paste0("File save to ",out))
