h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
QCexclude=c()

args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v5/Standing.height-BOLT-LMM-v5.out","2","plots","-ymax","50","-qc","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-recombRegions","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v6a/Standing.height-BOLT-LMM-v6a-chr2.out","-sample","/well/ukbiobank/phasing/final/phased_chunks/chr2.test2.sample","-dontComputeLD")


print(args)

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

plotRegions=TRUE
computeLD=TRUE
if("-dontPlotRegions"%in%args) plotRegions=FALSE
if("-dontComputeLD"%in%args) computeLD=FALSE

chroms = parse.range.string(chroms)

# any specified bgen or genotypes file for computing LD?
bgenData=args[which(args=="-bgenFile")+1]

# should we compute ld stats data and save it?
saveData=FALSE
if(computeLD){
    samplesFile = args[which(args=="-phenoFile")+1]
    system(paste0("cut -f1,2 -d' ' ",samplesFile," | tail -n +2 > ",samplesFile,".samples.tmp"))
    saveData=TRUE
}

system( paste0("mkdir ",plotOutDir,"/Region.plots") )



######## Get GWAS catalogue information
# NOTE: field descriptons are here: http://genome.ucsc.edu/cgi-bin/hgTables
catFile = NULL
if("-hits"%in%args) {
    # NOTE: this overrides the -hi for QC colours
    catInfo = args[which(args=="-hits")+1]
    load(catInfo,verbose=TRUE)
    
    catFile = catPheno  # just europeans
    colnames(catFile)[ncol(catFile)] = "Ancestry"
        
    catFile$Pvalue = catFile$V18
    catFile$SNP = catFile$V5
    catFile$BP = catFile$V4 # this is the chromEnd field
    catFile$CHR = gsub("chr","",catFile$V2)
    catFile$CHR[catFile$CHR%in%names(sexChroms)] = sexChroms[catFile$CHR[catFile$CHR%in%names(sexChroms)]]
    catFile$CHR = as.numeric(catFile$CHR)
    print( head(catFile) )
}


# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

if("-lreg" %in% args ) {
    useLmmInf = FALSE
    extraTitle = paste0(extraTitle,"-lreg")
} else {
    useLmmInf=TRUE
}

### read in GWAS results
if("-qc" %in% args){
                                        # QC THRESHOLDS  (used in plots made before 6-07-2016)
    minmaf = 0.001
    mininfo = 0.3
    maxmiss = 0.05  # maximum 5% missing data    
} else {
    minmaf=0
    mininfo=0
    maxmiss=1
}


resultsGENO = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF1 = resultsGENO$DF

resultsIMP = read.gwas(bgenData,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF2 = resultsIMP$DF


### read in recombination rates
print('Reading recombination rate files from /well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr*')

recombrates = list()
for(i in c(1:22,"X_nonPAR","X_PAR1","X_PAR2")){    
    print(i)
    r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
    recombrates[[i]] = r
}


###################################
# for each chromosome, get set of 'hits' based on recombination distances

# does this data exist already, if we haven't asked to computeLD? If we have asked to computeLD (default) then we run hits and LD functions.

if(!computeLD){
    hitsDataFile = paste0(resultsGENO$Pvalset,"-hitsAndLDcalc.RData")
    if("-ldRData"%in%args) hitsDataFile = args[which(args=="-ldRData")+1]
    if(hitsDataFile%in%list.files()) load(hitsDataFile,verbose=TRUE)
}





###################################                 
# plotting regions separately

                                        # get gene information
    source("/well/donnelly/ukbiobank_project_8874/clare/sex/scripts/commonScripts/plotGenesFunction.R")
    genes = load.genes(refGeneFile,fieldNames=fieldNames)


for (chrom in chroms){

    print(chrom)
                                        # get exon information
    exons = get.exons(genes,chromosome = chrom)

    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    HITS = HITSchrom[[as.character(chrom)]]
    ld = LDchrom[[as.character(chrom)]]
    recomb = recombrates[[as.character(chrom)]]
    
    print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )
#### plot the regions
    if(length(HITS$topVariants) > 0){
        sapply(1:length(HITS$topVariants),FUN=function(s) {
                                        # s = 17;
            #for(s in h) source(s)
     #      s=49; 

            regionToPlot = c(HITS$regs[s,1],HITS$regs[s,2])
                
            plot.gwas.region.v2(DF1,HITS$topVariants[s],catFile,ld=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG1"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=NULL)
            dev.off()

            plot.gwas.region.v2(DF1,HITS$topVariants[s],catFile,ld=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG2"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=DF2)
            dev.off()
            

        })
    }
}



###################################
# Read in GIANT results

#giantFile ='/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt'

giantFile ='/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz'    

giant = read.table(gzfile(giantFile),stringsAsFactors=FALSE,header=TRUE)

# get positions from ucsc
snpDBFile = '/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/ucsc.snps.annot/snp147.subsetFields.txt'
snpDF = read.table(snpDBFile,stringsAsFactors=FALSE,header=FALSE,colClasses=c("character","numeric","numeric","character","character","character"))


# include position in giant file
giantIndex = match(giant$MarkerName,snpDF$V4)
giant$BP = snpDF$V3[giantIndex] # end position is the right one!
colnames(giant)[c(1,7)] = c("SNP","P2")
giant$P = giant$P2
giant$CHRraw = snpDF$V1[giantIndex]
giant$CHR = as.numeric(gsub("chr","",giant$CHRraw))
giant$CHR[is.na(giant$CHR)] = 0
giant$MAF = as.numeric(giant$Freq.Allele1.HapMapCEU)
# exclude snps not found in uscs data
giant2 = giant[!is.na(giantIndex),]

giantData = list("DF"=giant2,"Pvalset"="GIANT.chrgenome")

###################################
# Plot GIANT results under genotype results

#NOTE: ASHG2 are plots with imputed data in gray as well.

for (chrom in chroms){

    print(chrom)
                                        # get exon information
    exons = get.exons(genes,chromosome = chrom)

    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    HITS = HITSchrom[[chrom]]
    ld = LDchrom[[chrom]]
    recomb = recombrates[[as.character(chrom)]]
    
    print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )
    
#### plot the regions
    if(length(HITS$topVariants) > 0){
        sapply(1:length(HITS$topVariants),FUN=function(s) {
                                        # s = 17;
            #for(s in h) source(s)
     #      s=49; 

            regionToPlot = c(HITS$regs[s,1],HITS$regs[s,2])
            
            plot.gwas.region.v2(resultsGENO$DF,HITS$topVariants[s],catFile=NULL,ld=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG2.withGIANT"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=resultsIMP$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantColour)
            dev.off()
            

        })
    }


###################################
# Plot GIANT results manhattan
sum(catFile$V6 == 25282103) # latest GIANT study
giantColour = "dodgerblue2"
    
plot.BOLT.pvalues(Data=giantData,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,extraFilename=".withCatalogueHits",col=giantColour)
    
plot.BOLT.pvalues(Data=giantData,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,col=giantColour)


    
###################################
# Plot GIANT results manhattan with ours

    # genotypes
plot.BOLT.pvalues(Data=resultsGENO,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=".ASHG.withGIANT")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT","UKBiobank genotyped"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()
plot.BOLT.pvalues(Data=resultsGENO,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename="-EuroHits.ASHG.withGIANT")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT","UKBiobank genotyped"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()


    # Imputed
plot.BOLT.pvalues(Data=resultsIMP,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=".ASHG.withGIANT",col="black")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT","UKBiobank imputed"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()
plot.BOLT.pvalues(Data=resultsIMP,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename="-EuroHits.ASHG.withGIANT",col="black")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT","UKBiobank imputed"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()


    
###################################
# Plot GIANT results manhattan with ours zoomed-in, with imputed snps too (actually done this for all of the hits, above.)

    
snp = "Affx-19070092"
s = which(HITS$topVariants==snp)
    
    regionToPlot = c(HITS$regs[s,1]-(10^6)/2,HITS$regs[s,2]+(10^6)/2)

                plot.gwas.region.v2(resultsGENO$DF,HITS$topVariants[s],catFile=NULL,ld=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG2.withGIANTwide"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=NULL,
                                basicCol="black",
                                imputedSNPs=resultsIMP$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantColour)
            dev.off()


# find stuff in gene
gene = "ACP1"
gRange = c(exons$txStart[exons$name2=="ACP1"][1],
exons$txEnd[exons$name2=="ACP1"][1])
    
thisGiant = giantData$DF[(giantData$DF$BP <= gRange[2])&(giantData$DF$BP >= gRange[1]) &(giantData$DF$CHR==chrom),]
thisImp = resultsIMP$DF[(resultsIMP$DF$BP <= gRange[2])&(resultsIMP$DF$BP >= gRange[1]) &(resultsIMP$DF$CHR==chrom),]
thisGeno = resultsGENO$DF[(resultsGENO$DF$BP <= gRange[2])&(resultsGENO$DF$BP >= gRange[1]) &(resultsGENO$DF$CHR==chrom),]

write.table(thisGiant,file=paste0(giantData$Pvalset,"-gene-",gene,"-stats.txt"),col.names=TRUE,quote=FALSE,row.names=FALSE)

write.table(thisImp,file=paste0(resultsIMP$Pvalset,"-gene-",gene,"-stats.txt"),col.names=TRUE,quote=FALSE,row.names=FALSE)

write.table(thisGeno,file=paste0(resultsGENO$Pvalset,"-gene-",gene,"-stats.txt"),col.names=TRUE,quote=FALSE,row.names=FALSE)    
