h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

args = commandArgs(TRUE)

print(args)

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""
useLmmInf=TRUE
if("-lreg"%in%args){
    useLmmInf=FALSE # use linear regression stats instead of lmm results
    extraTitle = paste0(extraTitle,"-lreg")
}

# do we fix the y-axis?
Ymax = NULL
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

# which chroms?
if(chroms!="genome") {
    if(chroms=="all") {
        
        chroms = 1:22
        if( "-sex" %in% args ) chroms = 23:26
        
    } else {
        
        chroms = parse.range.string(chroms)        
    }
}


# do we apply qc thresholds?
if("-qc" %in% args){
                                        # QC THRESHOLDS 
        minmaf = 0.001
        mininfo = 0.3
        maxmiss = 0.05  # maximum 5% missing data
            
    } else {
        minmaf=0
        mininfo=0
        maxmiss=1
    }


if("-minmaf" %in% args) minmaf = as.numeric(args[which(args=="-minmaf")+1])
if("-mininfo" %in% args) mininfo = as.numeric(args[which(args=="-mininfo")+1])
if("-maxmiss" %in% args) maxmiss = as.numeric(args[which(args=="-maxmiss")+1])




################ Actual plotting: run through chromosomes

print("You asked for chromosomes")
print( chroms )

if(!grepl("%%",dataFile)){
    # read in the whole set (could be multiple chromosomes)
    Data = read.gwas(dataFile,chrom="genome",minmaf,mininfo,maxmiss,Ymax,QCexclusions,extraTitle,useLmmInf=useLmmInf)
    rawPvalset=Data$Pvalset
    # are all the chromosomes in this dataset?
    chromsData = unique(Data$DF$CHR)
    if(chroms[1]!="genome") chroms = chroms[chroms%in%chromsData]
}

print("printing the following chromosomes")
print( chroms )

for(chrom in chroms){

    print(chrom)
        
    if(grepl("%%",dataFile)) {
        Data = read.gwas(gsub("%%",chrom,dataFile),chrom,minmaf,mininfo,maxmiss,Ymax,QCexclusions,extraTitle,useLmmInf=useLmmInf)
    } 
        
    plot.NULL.qqplot(Data,chrom,plotOutDir,Ymax,NULLdata=)
    
}


print(paste0("SEHR GUT! Plots saved in files with this kind of name: ",Data$Pvalset))

