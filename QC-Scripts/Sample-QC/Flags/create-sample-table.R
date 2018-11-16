#################
# This script is IMPORTANT. It creates the sample table for release to researchers.
#################

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-qcflags","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/Flags/b1__b11-b001__b095-sampleQC-flags.txt","-fam", "/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.fam","-pcs","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-pca-UKbio-round2","-wb","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt","-het","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-imiss.RData","-rel", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0","-out","b1__b11-b001__b095-sampleTable_v1")

print(args)

options = c("-fam","-wb","-qcflags","-pcs","-het","-rel","-out")
h = args[-c(which(args%in%options),1+which(args%in%options))]
for(helperScript in h){
    source(helperScript)
}


otherInfo = read.multiple.batch.info(batchInfoFields,printBatches=FALSE)
otherInfo = tbl_df(otherInfo)

outName = args[which(args=="-out")+1]

#############
# 1. QC FLAGS

# Read in flags table (this is without duplicates).
flagsFile = args[which(args=="-qcflags")+1]
flags = read.table(flagsFile,header=TRUE)


#############
# 2.  PCA
pcFile = paste0( args[which(args=="-pcs")+1],".RData" )
load(pcFile,verbose=TRUE)


#############
# 3.  White British
wbFile = args[which(args=="-wb")+1]
wb = read.table(wbFile,header=FALSE)[,1]


#############
# 4.  Samples in genotype data (from Colin)
famFile = args[which(args=="-fam")+1]
fam = read.table(famFile,header=FALSE)
colnames(fam)[1:2]=c("BPW","PIID")


#############
# 5.  Het corrected and missing
hetFile =  args[which(args=="-het")+1]
load(hetFile,verbose=TRUE)


#############
# 6.  Kinship table
relFile = args[which(args=="-rel")+1]
kinship = read.table(relFile,header=TRUE)


#############
# 7.  Who was filtered out of the kinship table??
relFileExcl = gsub(".kin0","-sample-effective_exclusions.txt",relFile) # exclusions from kinship (see filter-king.R)
kinshipExcl = read.table(relFileExcl,header=FALSE)[,1]



#############
# 8.  Who was used in the PC-calculation??
pcFileRaw = paste0( gsub("-pca-UKbio-","-autosome-sampleqc-fastpca-",args[which(args=="-pcs")+1]),"-highquality-pruned.fam" )
print("PCA pruned file: ")
print(pcFileRaw)
pcaSamples = read.table(pcFileRaw,header=TRUE,stringsAsFactors=FALSE)[,1]



print( paste0( "nInds in flag file: ",dim(flags)[1]) )
print( paste0( "nInds in flag file but not in fam file: ", sum(!flags[,1]%in%fam[,2]) ) )
print( paste0( "nInds wb file: ",length(wb) ) )
print( paste0( "nInds wb file but not in fam file: ",sum(!wb%in%fam[,2]) ) )
print( paste0( "nInds used for pcs: ",length(pcaSamples) ) )
print( paste0( "nInds in pc file but not in fam file: ", sum(!pcaSamples%in%fam[,2]) ) )
print( paste0( "nInds in fam file: ",dim(fam)[1]) )




############
# Create giant table! (base it on OtherInfo)

# fold in set of samples in genotype file (so they're in the same order)
outTable = left_join(fam[,1:2],otherInfo,by=c("PIID"="PIID")) 

# check batchPlateWell ids match up (they won't now, because colin has changed the .fam file to have just PIID PIID in the first two columns)
bpw = paste0(outTable$Batch,"-",outTable$Processed.Plate,"-",outTable$Processed.Well)
print( sum(bpw!=outTable$BPW) )


# fold in pcs
outTable = left_join(outTable,PCs[,!colnames(PCs)%in%c("Colors","Pops","Chars")],by=c("PIID"="PIID"))

# fold in who was used in PCA
outTable[["used.in.pca.calculation"]] = 0
outTable[["used.in.pca.calculation"]][ outTable[["PIID"]]%in%pcaSamples ] = 1

# fold in flags
outTable = left_join(outTable,flags,by=c("PIID"="PIID"))
for(i in colnames(flags)[-1]){
    # set to zero anything that wasn't in the flag file
    outTable[[i]][is.na(outTable[[i]])] = 0
}

# fold in white british
outTable[["white.british"]] = 0 
outTable[["white.british"]][ outTable[["PIID"]]%in%wb ] = 1


# fold in heterozygosity and missing rates
outTable = left_join(outTable,Table[,c("IID","miss","het","het.corrected")],by=c("PIID"="IID"))


# create and 'array' variable
outTable$genotyping.array = "UKBB"
outTable$genotyping.array[outTable$Batch%in%ukbileve.batches()] = "UKBL"



# create Kinship table for release
outKin = kinship[( kinship$ID1%in%outTable$PIID )&( kinship$ID2%in%outTable$PIID ),c("ID1","ID2","HetHet","IBS0","Kinship") ]
allKin = unique(c(outKin$ID1,outKin$ID2))

# create and 'in.kinship.table' variable
outTable$in.kinship.table = 0
outTable$in.kinship.table[outTable$PIID %in% allKin] = 1

# create and 'excluded.from.kinship.inference' variable
outTable$excluded.from.kinship.inference = 0
outTable$excluded.from.kinship.inference[outTable$PIID %in% kinshipExcl] = 1
    


############
# Rename some columns
colnames(outTable)[colnames(outTable)=="white.british"] = "in.white.British.ancestry.subset"
colnames(outTable)[colnames(outTable)=="miss"] = "sample.qc.missing.rate"  # this is the missing rate that was used to compute het/missing outliers. Computed on set of snps used in sample QC.
colnames(outTable)[colnames(outTable)=="het.corrected"] = "heterozygosity.pc.corrected"
colnames(outTable)[colnames(outTable)=="het"] = "heterozygosity"
colnames(outTable)[colnames(outTable)=="BPW"] = "batch.plate.well"
colnames(outTable)[colnames(outTable)=="Processed.Plate.Name"] = "Plate.Name" # change to match discription of batch/plate effects etc.
colnames(outTable)[colnames(outTable)=="Processed.Well"] = "Well"
colnames(outTable)[colnames(outTable)=="possible.sex.aneuploidy"] = "putative.sex.chromosome.aneuploidy"
print(colnames(outTable))


############
# Check that the order matches the fam file
print('checking matches with .fam file...') 
sum( fam[,1]!= outTable$BPW )
sum( fam[,2]!= outTable$Best.Array )



############
# Save this table with all columns - in text and R format (along with names of files used to put into the table!). Note: order of columns may not match those in the file for release.

outFile = paste0( baseSampleQCDir,"/data/ForRelease/" ,outName )

save(outTable, outKin, relFile, hetFile, famFile, flagsFile, pcFile, wbFile, file=paste0(outFile,"_allColumns.RData") )
# write a comma-separated file as there are spaces in some fields <-- which ones?
write.csv(outTable, paste0(outFile,"_allColumns.csv"),quote=FALSE,row.names=FALSE)

# write kinship table for output
write.csv(outKin, paste0(outFile,"_Kinship.csv"),quote=FALSE,row.names=FALSE)
# space-separated file too
write.table(outKin, paste0(outFile,"_Kinship.txt"),sep=" ",quote=FALSE,row.names=FALSE)



############
# Save this table with selected columns - in text and R format

columnsToRelease = c("Best.Array",batchInfoFields[5],"genotyping.array","Batch","Plate.Name","Well",batchInfoFields[c(19,17,7,23,24,28,29,1,2)],
    "sample.qc.missing.rate",
    "heterozygosity",
    "heterozygosity.pc.corrected",
    "het.missing.outliers",
    "putative.sex.chromosome.aneuploidy",
    "in.kinship.table",
    "excluded.from.kinship.inference",
    "excess.relatives",
    "in.white.British.ancestry.subset",
    "used.in.pca.calculation",
    grep("^PC",colnames(outTable),value=TRUE))

outTable2 = outTable[,columnsToRelease]

# write a comma-separated file.
write.csv(outTable2, paste0(outFile,"_releaseColumns.csv"),quote=FALSE,row.names=FALSE)

# write a space-separated file.
write.table(outTable2, paste0(outFile,"_releaseColumns.txt"),sep=" ",quote=FALSE,row.names=FALSE)


# print table with column descriptions
source( paste0(baseSampleQCDir,'/QC-Scripts/Sample-QC/Flags/SampleQCFieldDescriptions.R') )
desc = t( sapply(colnames(outTable2), function(co){
    if( grepl("PC",co) ){
        pc = gsub("PC","",co)
        d = get( "PC1" )
        description = gsub("1",pc,d)
    } else {
        description = get( co )
    }
}) )

write.table(desc[colnames(outTable2),],paste0(outFile,"_releaseColumns_FieldDescriptions.txt"),sep=" ",quote=TRUE,row.names=TRUE,col.names=FALSE)



desc = t( sapply(colnames(outKin), function(co){
        description = get( co )
    }) )

write.table(desc[colnames(outKin),],paste0(outFile,"_Kinship_FieldDescriptions.txt"),sep=" ",quote=TRUE,row.names=TRUE,col.names=FALSE)


# can it be read properly? ==> YES, but heterozygosity and hetpccorrected are different by machine precision (e-16). This is because the original files with this information were RData files, rather than text files... probably doesn't matter.
#test = read.table( paste0(outFile,"_releaseColumns.txt"),header=TRUE)
#j = sapply(1:ncol(test),function(x) sum(test[,x]!=outTable2[,x]))




############
# Print some summaries of the fields

print( paste0( "n samples in output table: ",dim(outTable)[1]) ) 

# 488,410 samples in expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.fam



for(i in colnames(outTable2) ){
    print(i)
    d = outTable2[[i]]
    if( ( length(unique(d)) > 20 ) & class(d) == "numeric" ) {
        print( paste0( "range: " ,range(d ,na.rm=TRUE) ) )
        print( paste0( sum(is.na(d))," NAs") )
    }
    if( ( length(unique(d)) > 20 ) & class(d) != "numeric" ) {
        print( head(d) )
        print( paste0( sum(is.na(d))," NAs") )
    }
    if( length(unique(d)) <= 20 ) print( table(d,useNA="ifany") )        
}



print(paste0("DONE!! Tables saved as: ",outFile,"*"))


