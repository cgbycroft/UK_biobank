
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)




######## Which chromosome and what data are we plotting?

args = commandArgs(TRUE)

#args = c("test.out","1","plots", "-ymax","50","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-Euro-hits")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v16/Internal.Pico..ng.uL.-BOLT-LMM-v16.out","genome","plots", "-ymax","50","-qc","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Internal.Pico..ng.uL..RData", "-title", "-Euro-hits")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s-23.outt","23","plots", "-ymax","50","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-raw-Euro-hits")

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

######## Get GWAS catalogue information
# NOTE: field descriptons are here: http://genome.ucsc.edu/cgi-bin/hgTables
columns = c("bin","chrom",
"chromStart",
"chromEnd",
"name",
"pubMedID",
"author",
"pubDate",
"journal",
"title",
"trait",
"initSample",
"replSample",
"region",
"genes",
"riskAllele",
"riskAlFreq",
"pValue",
"pValueDesc",
"orOrBeta",
"ci95",
"platform",
"cnv")

catFile = NULL
if("-hits"%in%args) {
    # NOTE: this overrides the -hi for QC colours
    catInfo = args[which(args=="-hits")+1]
    load(catInfo,verbose=TRUE)
    
    catFile = catPheno  # any ancestry
    colnames(catFile)[ncol(catFile)] = "Ancestry"
        
    catFile$Pvalue = catFile$V18
    catFile$SNP = catFile$V5
    catFile$BP = catFile$V4 # this is the chromEnd field
    catFile$CHR = gsub("chr","",catFile$V2)
    catFile$CHR[catFile$CHR%in%names(sexChroms)] = sexChroms[catFile$CHR[catFile$CHR%in%names(sexChroms)]]
    catFile$CHR = as.numeric(catFile$CHR)

    # adjust any Xchrom hits for PAR regions
    parHits = ( catFile$CHR == 23 ) & ( ( catFile$BP < par1 ) | ( catFile$BP > par2 ) )
    print( paste0( sum(parHits)," hits in PAR X") )
    catFile$CHR[parHits] = 25
    
    print( head(catFile) )
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


# Should we highlight some SNPs for QC reasons?
highlightSNPs=NULL
highlightCols=NULL
QCexclusions=c()

#if(( "-hi"%in%args )|("-qc" %in%args )){

#    print("Reading QC snps lists...")

#    if( "-sex" %in% args ) type = "sexchrom"  else type = "autosome"
#    QCSNPList = read.SNPQC.files(justSNPs=TRUE,type=type)

                                        #    QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$imageSNPs,QCSNPList$concordanceSNPs))    
                                        #QCexclusions = unique(c(QCSNPList$arraySNPs,QCSNPList$concordanceSNPs))    # this is only required if using -hi or clare's versions of plink genotype files.
#    QCexclusions = unique(QCSNPList$concordanceSNPs)
    
#    if( "-hi"%in%args ){
#        highlightSNPs = unique(unlist(QCSNPList))
#        colors = rep("black",length(highlightSNPs))

#        colors[highlightSNPs%in%QCSNPList$batchHweSNPs] = "green"   # HWE (apply first)
#        colors[highlightSNPs%in%c(QCSNPList$plateSNPs,QCSNPList$batchSNPs)] = "purple"  # BATCH/PLATE
#        colors[highlightSNPs%in%c(QCSNPList$imageSNPs)] = "orange" # IMAGE ARTEFACT
#        colors[highlightSNPs%in%c(QCSNPList$arraySNPs)] = "red" # ARRAY
#        colors[highlightSNPs%in%c(QCSNPList$concordanceSNPs)] = "blue" # CONCORDANCE
#        highlightCols = colors
        
#        print(table(highlightCols))
#    }
    
#}


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
        
    plot.BOLT.pvalues(Data,chrom,plotOutDir,plotQQ=TRUE,catFile,Ymax,highlight=highlightSNPs,highlightCols=highlightCols)
    
}

#dev.off()

print(paste0("SEHR GUT! Plots saved in files with this kind of name: ",Data$Pvalset))



# testing
# Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v5/Standing.height-BOLT-LMM-v5.out 19,21-22 plots -ymax 50 -title -test

# Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v5/Standing.height-BOLT-LMM-v5.out genome plots -ymax 50 -title -test

# Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v7/Standing.height-BOLT-LMM-v7.out 21-22 plots -ymax 50 -hits /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData -title -test-raw-hits

