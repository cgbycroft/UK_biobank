#######
# This script creates a list of samples to exclude in the pca-UKBiobank analysis, based on relatedness and missingness

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-k","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0","-ibs0","0.002","-out","pca-UKBio-sample-exclusions-init.txt")

h = args[-c(which(args%in%c("-k","-out","-ibs0","-exclude")),1+which(args%in%c("-k","-out","-ibs0","-exclude")))]
for(helperScript in h){
    source(helperScript)
}

outputFile = args[which(args%in%c("-out")) + 1]
relatednessFile = args[which(args%in%c("-k")) + 1]
ibs0Threshold = args[which(args%in%c("-ibs0")) + 1]
missingRatesFile = paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc.imiss.gz")

miss = read.table(gzfile(missingRatesFile),stringsAsFactors=FALSE,header=T)
gender = read.multiple.batch.info(c("Submitted.Gender","Inferred.Gender"))


#### Interim release documentation:
#We also removed samples who were related to multiple other samples (to the 1st, 2nd or 3rd degree), one sample from each remaining related pair (chosen randomly), as well as removing all twins and gender mismatches and samples with a high missing rate
## ---> This version we just exclude 1st and 2nd-degree relatives (based on first round of KING) to see if it works!

#### add samples to pca exclude list

# REFS
referenceList="../../referenceList.txt"
references = read.table(referenceList,header=FALSE,stringsAsFactors=FALSE)[,1]

# Gender mismatches
genderMismatch = gender$PIID[gender$Submitted.Gender!=gender$Inferred.Gender]

# High missing rates
highMissing = miss$IID[miss$F_MISS > 0.02]

# other exclusions (e.g het/missing list) if asked for.
exclusions = c()
if("-exclude"%in%args) {
    exFiles = args[which(args=="-exclude")+1]
    for(file in exFiles){
        print( paste0("Reading: ",file))
        e = read.table(file,stringsAsFactors=FALSE)[,1]
        print( paste0( length(e)," samples.") )
        exclusions = c(exclusions,e)
    }

    exclusions = unique(exclusions)
    print(paste0(length(exclusions)," samples in the extra exclusions list(s)"))
}
    
# Relatives
#first exclude the above samples from kinship table
badSamples = unique(c(references,genderMismatch,highMissing,exclusions))
print( paste0(length(badSamples)," excluded for QC reasons.") )
write.table(badSamples,file="badSamples.tmp",quote=FALSE,row.names=FALSE,col.names=FALSE)


# find unrelated samples : after excluding the above individuals
args <- c(h,"-in",relatednessFile,"-ibs0",ibs0Threshold,"-exclude","badSamples.tmp","-outdir",".")
source( paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/postking/get-unrelated.R")) 


#### Combine all filters
# NOTE: this list has taken into account the earlier filters.
unrelated = read.table(paste0(gsub(".kin0","",basename(relatednessFile)),"-unrelated.txt" ),stringsAsFactors=FALSE,header=FALSE)[,1]

toExcludePCA = otherInfo$PIID[!otherInfo$PIID%in%unrelated]


print( paste0(length(toExcludePCA)," samples in pca exclusion list") )
print( paste0(length(references)," references" ) )
print( paste0(length(genderMismatch)," gender mismatches" ) )
print( paste0(length(highMissing)," missing rate > 0.02" ) )


# check that there are no more related pairs in what's left over, after excluding poor quality samples from the relatedness table
allSamples = gender$PIID
left = allSamples[!allSamples%in%toExcludePCA]

print( paste0( "There should be ",length(left), " samples in the pca" ) )


write.table(cbind(toExcludePCA,toExcludePCA),file=outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE)

print( paste0("DONE! List of samples to exclude in ",outputFile) )
