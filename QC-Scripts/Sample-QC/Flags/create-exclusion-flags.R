args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-hm","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt","-rel","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-excess-relatives.txt","-dup", "/well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt")

print(args)

h = args[-c(which(args%in%c("-hm","-rel1","-rel2","-sex","-dup")),1+which(args%in%c("-hm","-rel1","-rel2","-sex","-dup")))]
for(helperScript in h){
    source(helperScript)
}


otherInfo = read.multiple.batch.info(batchInfoFields)


# Het/missing outliers
hmFile = args[which(args=="-hm")+1]
hm = read.table(hmFile,header=FALSE)[,1]


# Excess related
relFile1 = args[which(args=="-rel1")+1]
rel1 = read.table(relFile1,header=FALSE)[,1]

# Related exclusions
relFile2 = args[which(args=="-rel2")+1]
rel2 = read.table(relFile2,header=FALSE)[,1]


# Dubious sex calling (not necessarily sex mismatches)
dubSexFile = args[which(args=="-sex")+1]
dubSex = read.table(dubSexFile,header=FALSE,stringsAsFactors=FALSE)[,1]

# sex mismatch
sexMismatch = otherInfo$PIID[otherInfo$Inferred.Gender!=otherInfo$Submitted.Gender]


# Duplicates (to exclude from any list)
dupFile = args[which(args=="-dup")+1]
dup = read.table(dupFile,header=FALSE)[,1]


# create table
allSamples = unique(c(hm,rel1,rel2,dubSex,sexMismatch))

flagTable = as.data.frame(allSamples,stringsAsFactors=FALSE)

flagTable$het.missing.outliers = 0
flagTable$het.missing.outliers[flagTable$allSamples%in%hm] = 1

flagTable$excess.relatives = 0
flagTable$excess.relatives[flagTable$allSamples%in%rel1] = 1 # anyone with more than 10 3rd-degree relatives in kinship table (188 samples)

flagTable$excluded.from.kinship = 0
flagTable$excluded.from.kinship[flagTable$allSamples%in%rel2] = 1

flagTable$possible.sex.aneuploidy = 0
flagTable$possible.sex.aneuploidy[flagTable$allSamples%in%dubSex] = 1

flagTable$sex.mismatch = 0
flagTable$sex.mismatch[flagTable$allSamples%in%sexMismatch] = 1

# exclude duplicates
flagTable2 = flagTable[!flagTable$allSamples%in%dup,]

colnames(flagTable)[1] = colnames(flagTable2)[1] = "PIID"


tabulate <- function(data){
    out = diag(colSums(data))
    for( i in which(diag(out)!=0) ){
        for( j in which(diag(out)!=0) ){
            if(i==j) next
            tab = table(data[,c(i,j)])
            out[i,j] = tab["1","1"]
        }
    }
    colnames(out) = rownames(out) = colnames(data)
    return(out)
}

print("with dupes")
print( tabulate(flagTable[,-1]) )

print("exclude dupes")
print( tabulate(flagTable2[,-1]) )

print( paste0("Total non-duplicate samples flagged: ",nrow(flagTable2) ) )


# write the table
outName = str_split(basename(hmFile),"-autosome-sampleqc-hetcorrected-")[[1]][1]

write.table(flagTable2,file=paste0(outName,"-sampleQC-flags.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(flagTable,file=paste0(outName,"-sampleQC-flags-with-dupes.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)

print("Files saved here:")
print( paste0(outName,"-sampleQC-flags*") )

print( "See Flags/Logs for table of numbers.")
