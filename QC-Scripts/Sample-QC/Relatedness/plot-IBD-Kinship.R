## Script to plot output from Plink IBD, using same set of snps as for final KING run + only samples that ended up in kinship table.


args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-rel","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0" ,"-ibd" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome","-ibs0","0.0012")

print(args)
h = args[-c(which(args%in%c("-rel","-ibd","-ibs0")),1+which(args%in%c("-rel","-ibd","-ibs0")))]
for(helperScript in h){
    source(helperScript)
}

KingOutFile = args[which(args=="-rel")+1]
PlinkOutFile = args[which(args=="-ibd")+1]
OutputFile = basename(PlinkOutFile)
OutDir = "plots"
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])


## get sample info
otherInfo = read.multiple.batch.info(c("Ethnic.background","Place.of.birth.in.UK...north.co.ordinate","Pops","Chars","Colors"))


## get kinship data
print( paste0("Reading kinship data... ",KingOutFile) )
kinKing = read.table(KingOutFile,header=TRUE,stringsAsFactors=FALSE)

## get kinship classes
classes = get.kin.classes(kinKing,ibs0Threshold=ibs0Threshold)

## IBD Plink results
print( paste0("Reading plink ibd data... ",PlinkOutFile) )
files = list.files(path=dirname(PlinkOutFile),pattern=basename(PlinkOutFile))
header = read.table(paste0(PlinkOutFile,".1"),header=FALSE,nrow=1,stringsAsFactors=FALSE)
columnClasses = c(rep("character",6),rep("numeric",13))
columnClasses2 = columnClasses; columnClasses2[c(5:9,11:14,18,19)] = "NULL"

print(paste0("reading in ",length(files)," Plink IBD files..."))
                                        # NOTE: each pair of individuals will turn up twice! Need to find and remove duplicates.
# Outfile will look like the KING output.

outFileCombined = paste0(PlinkOutFile,".combined")
#write.table(header[,columnClasses2!="NULL"],file=outFileCombined,col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(t(c("FID1", "ID1", "FID2", "ID2", "N_SNP","HetHet", "IBS0", "Kinship")),file=outFileCombined,col.names=FALSE,row.names=FALSE,quote=FALSE)


myPairs=c()
for(file in files){
    print(file)
    if(grepl("genome.1$",file)) ibd = read.table(paste0(dirname(PlinkOutFile),"/",file),header=TRUE,stringsAsFactors=FALSE,colClasses=columnClasses2) else ibd = read.table(paste0(dirname(PlinkOutFile),"/",file),header=FALSE,stringsAsFactors=FALSE,colClasses=columnClasses2)    
                                        # colnames(ibd) = header
                                        # throw away any pair that's alread in our list
    AlreadyRead = (paste0(ibd[,1],"__",ibd[,3])%in%myPairs)|(paste0(ibd[,3],"__",ibd[,1])%in%myPairs)
    sum(!AlreadyRead)
    
    newPairs = paste0(ibd[!AlreadyRead,1],"__",ibd[!AlreadyRead,3])
    myPairs = c(myPairs,newPairs)
    kin = ibd[!AlreadyRead,]
    N_SNP = rowSums(kin[,6:8])
    IBS0 = kin[,6]/N_SNP # compute IBS0 from counts
    HetHet = kin[,7]/N_SNP
    Kinship = kin[,5]/2
    kin$N_SNP = N_SNP
    kin$HetHet = HetHet
    kin$IBS0 = IBS0
    kin$Kinship=Kinship
    
    write.table(kin[,c("V1","V2","V3","V4","N_SNP","HetHet","IBS0", "Kinship")],file=outFileCombined,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)

}


kinPlink = read.table(outFileCombined,header=TRUE,stringsAsFactors=FALSE,colClasses=columnClasses2)

