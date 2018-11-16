## Script to process output from Plink IBD, using same set of snps as for final KING run + only samples that ended up in kinship table. Output is a file that looks like the KING output, but has plink results in it.


args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-ibd" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome","-hm","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt")

print(args)
h = args[-c(which(args%in%c("-rel","-ibd","-ibs0","-hm")),1+which(args%in%c("-rel","-ibd","-ibs0","-hm")))]
for(helperScript in h){
    source(helperScript)
}

PlinkOutFile = args[which(args=="-ibd")+1]
OutputFile = basename(PlinkOutFile)

# het miss outliers
#hetmissoutliers = read.table(args[which(args=="-hm")+1],stringsAsFactors=FALSE,header=FALSE)[,1]

## IBD Plink results
print( paste0("Reading plink ibd data... ",PlinkOutFile) )
files = list.files(path=dirname(PlinkOutFile),pattern=basename(PlinkOutFile))
files = files[!grepl("combined",files)]

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
    if(grepl("genome.1$",file)) ibd = read.table(paste0(dirname(PlinkOutFile),"/",file),header=FALSE,skip=1,stringsAsFactors=FALSE,colClasses=columnClasses2) else ibd = read.table(paste0(dirname(PlinkOutFile),"/",file),header=FALSE,stringsAsFactors=FALSE,colClasses=columnClasses2)    
                                        # colnames(ibd) = header
    # exclude any pair that involves a het-miss outlier
    #hetMiss = (ibd[,1]%in%hetmissoutliers)|(ibd[,3]%in%hetmissoutliers)
    #ibd=ibd[!hetMiss,]
    # throw away any pair that's alread in our list
    AlreadyRead = (paste0(ibd[,1],"__",ibd[,3])%in%myPairs)|(paste0(ibd[,3],"__",ibd[,1])%in%myPairs)
    print( sum(!AlreadyRead) )
    print( sum(AlreadyRead) )
    
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
    #print(head(kin))
    
    write.table(kin[,c("V1","V2","V3","V4","N_SNP","HetHet","IBS0", "Kinship")],file=outFileCombined,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)

}

#myInd = "A550465-4196550-091314-229_C08"

#kinPlink = read.table(outFileCombined,header=TRUE,stringsAsFactors=FALSE,colClasses=columnClasses2)

print("File output to:")
print(outFileCombined)
