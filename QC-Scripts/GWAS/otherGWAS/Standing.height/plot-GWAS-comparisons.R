h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
QCexclude=c()

args = commandArgs(trailingOnly=TRUE)

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v14/Standing.height-BOLT-LMM-v14.out","2","plots","-ymax","50","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v13/Standing.height-BOLT-LMM-v13-chr2.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr2.hrc+uk10k.I4.v1.1.sample","-ldRData","Standing.height-BOLT-LMM-v14.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp","Standing.height-BOLT-LMM-v13-chr2.out.chr2.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData")

args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out","2","plots","-ymax","50","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15/Standing.height-BOLT-LMM-v15-chr2.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr2.hrc+uk10k.I4.v1.1.sample","-ldRData","Standing.height-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp","Standing.height-BOLT-LMM-v15-chr2.out.chr2.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData")


print(args)
print(getwd())

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chroms = parse.range.string(chroms)

# any specified bgen or genotypes file for computing LD?
bgenData=args[which(args=="-bgenFile")+1]

# phased file?
phaseFile = "/well/ukbiobank/expt/V2_QCed.export/data/imputation_pipeline_input/v4/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.phasing_QC.chr"

if("-phased"%in%args){
    phaseFile = args[which(args=="-phase")+1]
} 

# should we compute ld stats data and save it?
#saveData=FALSE
#if(computeLD){
#    samplesFile = args[which(args=="-phenoFile")+1]
#    system(paste0("cut -f1,2 -d' ' ",samplesFile," | tail -n +2 > ",samplesFile,".samples.tmp"))
#    saveData=TRUE
#}

system( paste0("mkdir ",plotOutDir,"/Results.comparison.plots") )

                                        # get gene information
source("/well/donnelly/ukbiobank_project_8874/clare/sex/scripts/commonScripts/plotGenesFunction.R")
genes = load.genes(refGeneFile,fieldNames=fieldNames)


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

###  Set QC thresholds!
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


if("-minmaf" %in% args) minmaf = as.numeric(args[which(args=="-minmaf")+1])
if("-mininfo" %in% args) mininfo = as.numeric(args[which(args=="-mininfo")+1])
if("-maxmiss" %in% args) maxmiss = as.numeric(args[which(args=="-maxmiss")+1])


###################################
# read GWAS results!

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


hitsDataFile = paste0(resultsGENO$Pvalset,"-hitsAndLDcalc.RData")
if("-ldRData"%in%args) hitsDataFile = args[which(args=="-ldRData")+1]
if(hitsDataFile%in%list.files()) {
    load(hitsDataFile,verbose=TRUE)
    HITSchromGENO = HITSchrom
    LDchromGENO = LDchrom
}


if("-ldRDataImp"%in%args) hitsDataFileImputed = args[which(args=="-ldRDataImp")+1]
if(hitsDataFileImputed%in%list.files()) {
    load(hitsDataFileImputed,verbose=TRUE)
    HITSchromIMP = HITSchrom
    LDchromIMP = LDchrom
}



###################################
# Get giant data

giantFile ='/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz'

giant = read.table(gzfile(giantFile),stringsAsFactors=FALSE,header=TRUE)

giant2 = giant[!is.na(giant$BP),]

                                        # cap the p-values
if(is.numeric(Ymax)){
    
    giant2$P2 = giant2$P
    giant2$P2[-log10(giant2$P)>Ymax] = 10^-Ymax
}

giantData = list("DF"=giant2,"Pvalset"="GIANT.chrgenome")





###################################
# Plot GIANT results manhattan

giantColour = "dodgerblue2"

for( chrom in 1:22 ){
    sum(catFile$V6 == 25282103) # latest GIANT study in NCBI
    
    plot.BOLT.pvalues(Data=giantData,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=".withCatalogueHits",col=giantColour)
    legend("topright",legend=c("GIANT 2014","GWAS catalogue hits 2016"),col=c(giantColour,"red"),bty="n",pch=c(16,16),cex=3,pt.lwd=c(1,1))
    dev.off()
    
    plot.BOLT.pvalues(Data=giantData,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=FALSE,col=giantColour)
    
}



###################################
# Plot GIANT results under genotype results

#NOTE: ASHG2 are plots with imputed data in gray as well.

for (chrom in chroms){

    print(chrom)
                                        # get exon information
    exons = get.exons(genes,chromosome = chrom)
   
    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    HITS = HITSchromGENO[[as.character(chrom)]]
    ld = LDchromGENO[[as.character(chrom)]]
    ldIMP = LDchromIMP[[as.character(chrom)]]

    recomb = recombrates[[as.character(chrom)]]
    
    print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )
    
#### plot the regions
    if(length(HITS$topVariants) > 0){
        sapply(1:length(HITS$topVariants),FUN=function(s) {
                                        # s = 17;
            #for(s in h) source(s)
     #      s=49; 
            # s=51
            regionToPlot = c(HITS$regs[s,1],HITS$regs[s,2])
            
            plot.gwas.region.v2(resultsGENO$DF,HITS$topVariants[s],catFile=NULL,ld=ld,ldIMP=ldIMP,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=NULL,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG2.withGIANT"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=resultsIMP$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantColour)
            dev.off()
            

        })
    }

}





###################################
# Plot GIANT results manhattan with ours

                                        # genotypes
plot.BOLT.pvalues(Data=resultsGENO,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=".ASHG.withGIANT")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT 2014","UKBiobank genotyped"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()
plot.BOLT.pvalues(Data=resultsGENO,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename="-EuroHits.ASHG.withGIANT")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT 2014","UKBiobank genotyped"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()


                                        # Imputed
plot.BOLT.pvalues(Data=resultsIMP,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=".ASHG.withGIANT",col="black")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT 2014","UKBiobank imputed"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()
plot.BOLT.pvalues(Data=resultsIMP,chrom,plotOutDir="plots",plotQQ=FALSE,catFile,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename="-EuroHits.ASHG.withGIANT",col="black")
points(giantData$DF$BP[giantData$DF$CHR==chrom],-log10(giantData$DF$P2[giantData$DF$CHR==chrom]),col=giantColour,pch=16)
legend("top",legend=c("GIANT 2014","UKBiobank imputed"),col=c(giantColour,"black"),bty="n",pch=16,cex=3)
dev.off()



###################################
# Plot GIANT results manhattan with ours zoomed-in on specific region where GIANT has significant p-values but we don't.

# SNP in GWAS catalogue: rs7567288-T at 134434823 by Lango Allen H 2010; European ancestry...

myRegion = as.data.frame(t(c(134434824 - 2.5e5,134434824 + 2.5e5))) # go half a megabase wide
myRegion2 = as.data.frame(t(c(134434824 - 2e6,134434824 + 2e6))) # go half a megabase wide
variant = resultsIMP$DF$SNP2[resultsIMP$DF$BP==134434824]

# first compute LD (or load it)

#load(paste0(resultsIMP$Pvalset,"-",variant,"-hitsAndLDcalc.RData"),verbose=TRUE)
load("Standing.height-BOLT-LMM-v13-chr2.out.chr2.maf0.info0.pruned-raw-lreg-2:134434824_T_C-hitsAndLDcalc.RData",verbose=TRUE)
ld = LDchrom[["2"]]

ld = compute.ld.around.snp.bgen(variants=variant,bgenFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr2.hrc+uk10k.I4.v1.1.bgen",bgenSampleFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr2.hrc+uk10k.I4.v1.1.sample",chrom=2,samplesFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt.samples.tmp",regions=myRegion,DF=resultsIMP$DF)

LDchrom=list()
LDchrom[["2"]] = ld
save(LDchrom,file=paste0(resultsIMP$Pvalset,"-",variant,"-hitsAndLDcalc.RData"))


# where are the genotyped snps?
GenoSNPsHere = resultsGENO$DF[(resultsGENO$DF$CHR==2)&(resultsGENO$DF$BP>=myRegion[1,1])&(resultsGENO$DF$BP<=myRegion[1,2]),]
ld = LDchrom[["2"]]
    
plot.gwas.region.v2(resultsGENO$DF,GenoSNPsHere$SNP[GenoSNPsHere$P==min(GenoSNPsHere$P)],catFile,ld=NULL,ldIMP=ld,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=c(myRegion2[1,1],myRegion2[1,2]),Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(resultsIMP$Pvalset,"-forASHG2.withGIANTwide"),addGenes=list("genes"=genes,"exons"=exons),thisChrom=2,
                    ldIMPvariant = variant,
                    showRegion=NULL,
                    basicCol="black",
                    imputedSNPs=resultsIMP$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantColour)
dev.off()

                                        # "Affx-17862683" = "rs7567288"

# plot with the 'north coordinate' gwas!
resultsGENOnorth = read.gwas(gsub("Standing.height","Place.of.birth.in.UK...north.co.ordinate",dataFile),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF1north = resultsGENOnorth$DF

resultsIMPnorth = read.gwas(gsub("Standing.height","Place.of.birth.in.UK...north.co.ordinate",bgenData),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF2north = resultsIMPnorth$DF


plot.gwas.region.v2(resultsGENOnorth$DF,GenoSNPsHere$SNP[GenoSNPsHere$P==min(GenoSNPsHere$P)],catFile,ld=NULL,ldIMP=ld,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=c(myRegion2[1,1],myRegion2[1,2]),Ymax=20,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(resultsIMPnorth$Pvalset,"-forASHG2.withGIANTwide"),addGenes=list("genes"=genes,"exons"=exons),thisChrom=2,
                    ldIMPvariant = variant,
                    showRegion=NULL,
                    basicCol="black",
                    imputedSNPs=resultsIMPnorth$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantColour)
dev.off()


## Read in pc-loads for white british
SnpLoads = '/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British.snpload.map'
loads = dplyr::tbl_df(read.table(SnpLoads,header=FALSE,stringsAsFactors=FALSE))

# plot chromosome 2 in specific region
toCheckRegion=c(myRegion2[1,1]-3e6,myRegion2[1,2]+3e6)
myLoads = loads[(loads$V1==chrom)&(loads$V4>=toCheckRegion[1])&(loads$V4<=toCheckRegion[2]),]

png(paste0("plots/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British.snploads-TEST.png"),width=41,height=12,units="in",res=150)

par(mar=c(5,5,4.1,5),cex.axis=1.5,cex.lab=1.5)

plot(NULL,xlim=toCheckRegion,ylim=c(-max(abs(myLoads[,8:27])),max(abs(myLoads[,8:27]))))

for( i in 1:20){
    points(myLoads$V4,myLoads[[i+7]],col=colors()[i],pch=16)    
}
abline(v=c(myRegion2[1,1],myRegion2[1,2]),col="red",lwd=2)
dev.off()


###################################
# Plot GIANT results manhattan with ours zoomed-in, with imputed snps too (actually done this for all of the hits, above.) THIS IS FOR THE ACP1 region


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


###################################
# Plot GIANT results p-values with ours
giantData$DF$SNP2a = altSNPID(giantData$DF$CHR,giantData$DF$BP,giantData$DF$Allele1,giantData$DF$Allele2)
giantData$DF$SNP2b = altSNPID(giantData$DF$CHR,giantData$DF$BP,giantData$DF$Allele2,giantData$DF$Allele1)

# genotypes
genoIndexA = match(DF1$SNP2,giantData$DF$SNP2a)
genoIndexB = match(DF1$SNP2,giantData$DF$SNP2b)
genoIndex = genoIndexA; genoIndex[is.na(genoIndexA)] = genoIndexB[is.na(genoIndexA)]
    
giantGeno = giantData$DF[genoIndex,] # Giant at the same genotypes. 
                                        # Ignore anything that does not match on one of the alleles
print( paste0( sum(!is.na(giantGeno$SNP))," snps in genotype data match with giant snps.") )
print( paste0( dim(DF1)[1]," snps in genotype data.") )

# imputed
impIndexA = match(DF2$SNP2,giantData$DF$SNP2a)
impIndexB = match(DF2$SNP2,giantData$DF$SNP2b)
impIndex = impIndexA; impIndex[is.na(impIndexA)] = impIndexB[is.na(impIndexA)]

giantImp = giantData$DF[impIndex,] # Giant at the same genotypes

print( paste0( sum(!is.na(giantImp$SNP))," snps in genotype data match with giant snps.") )
print( paste0( dim(DF2)[1]," snps in imputed data.") )


# with north coordinate
genoIndex2 = match(DF1$SNP,resultsGENOnorth$DF$SNP)
northGeno = resultsGENOnorth$DF[genoIndex2,] # Giant at the same genotypes
impIndex2 = match(DF2$SNP,resultsIMPnorth$DF$SNP)
northImp = resultsIMPnorth$DF[impIndex2,] # Giant at the same genotypes


###### 
phasedColour = "red"
rareColour = "blue"
infoColour = "orange"
rareThreshold = 0.001

#chrom=chroms
phased = read.phased.snps(chrom)[[1]]
bgenSNPs = read.imputed.snps(chrom)[[1]]  # get alternative ID from bgen files. It's a bit slow (~2mins) for large chromosomes. About as slow as reading the gwas output. SNPs should be in the same order.
DF2$alternate_ids = bgenSNPs$alternate_ids[match(DF2$SNP2,bgenSNPs$SNP2)]

a=match(DF2$SNP2,phased$SNP2) # more are matched this way...
b=match(DF2$SNP,phased$V2) # miss-matches due to multi-allelics
c=match(DF2$SNP2,phased$SNP3) # nothing matches with allele flips.
d = grep("2:",bgenSNPs$alternate_ids,invert=TRUE) # how many imputed snps have rsids or Affyids
e = match(bgenSNPs$alternate_ids,phased$V2) # how many imputed snps have rsids or Affyids
f = match(DF2$alternate_ids,phased$V2) # how many imputed snps have rsids or Affyids



rareGeno = DF1$MAF < rareThreshold
inPhasedGeno = DF1$SNP%in%DF2$alternate_ids
poorInfoGeno = DF1$F_MISS >= 0.05

rareImp = DF2$MAF < rareThreshold
inPhasedImp = DF2$alternate_ids%in%DF1$SNP
poorInfoImp = DF2$INFO < 0.3

colorsGeno = rep("black",dim(DF1)[1]) 
colorsGeno[inPhasedGeno]=phasedColour # in phased
colorsImp = rep("black",dim(DF2)[1]) 
colorsImp[inPhasedImp]=phasedColour # in phased


# scatter plots of pvalues
png(paste0( "plots/",resultsGENO$Pvalset,"-withGIANT-pvalues-scatter.png" ),width=1000,height=1000)

x = -log10(giantGeno$P2)
y = -log10(DF1$P2)

plot(x,y,ylab="-log10(p) UK Biobank genotype data",xlab="-log10(p) GIANT data",xlim=c(0,Ymax),ylim=c(0,Ymax),main=paste0(sum(!is.na(giantGeno$SNP))," overlapping markers."),pch=16,col=add.alpha(colorsGeno,0.5))
# phased
points(x[inPhasedGeno],y[inPhasedGeno],col=phasedColour,pch=16)
# rare
points(x[rareGeno],y[rareGeno],col=rareColour,pch=1)
# low info
points(x[poorInfoGeno],y[poorInfoGeno],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

dev.off()


png(paste0( "plots/",resultsGENO$Pvalset,"-withGIANT-pvaluesUncapped-scatter.png" ),width=1000,height=1000)
x = -log10(giantGeno$P)
y = -log10(DF1$P)
plot(x,y,ylab="-log10(p) UK Biobank genotype data",xlab="-log10(p) GIANT data",main=paste0(sum(!is.na(giantGeno$SNP))," overlapping markers."),pch=16,col=add.alpha(colorsGeno,0.5))
abline(0,1,col="red",lty=3)
# phased
points(x[inPhasedGeno],y[inPhasedGeno],col=phasedColour,pch=16)
# rare
points(x[rareGeno],y[rareGeno],col=rareColour,pch=1)
# low info
points(x[poorInfoGeno],y[poorInfoGeno],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

dev.off()

# just one chromosome
png(paste0( "plots/",gsub("genome",chrom,resultsGENO$Pvalset),"-withGIANT-pvalues-scatter.png" ),width=1000,height=1000)
x = -log10(giantGeno$P2)
y = -log10(DF1$P2)
plot(x[giantGeno$CHR==2],y[giantGeno$CHR==2],ylab="-log10(p) UK Biobank genotype data",xlab="-log10(p) GIANT data",xlim=c(0,Ymax),ylim=c(0,Ymax),main=paste0(sum(!is.na(giantGeno$SNP[giantGeno$CHR==2]))," overlapping markers."),pch=16,col=add.alpha(colorsGeno[giantGeno$CHR==2],0.5))

# phased
points(x[(inPhasedGeno)&(giantGeno$CHR==2)],y[(inPhasedGeno)&(giantGeno$CHR==2)],col=phasedColour,pch=16)
# rare
points(x[(rareGeno)&(giantGeno$CHR==2)],y[(rareGeno)&(giantGeno$CHR==2)],col=rareColour,pch=1)
# low info
points(x[(poorInfoGeno)&(giantGeno$CHR==2)],y[(poorInfoGeno)&(giantGeno$CHR==2)],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

abline(0,1,col="red",lty=3)
dev.off()


png(paste0( "plots/",gsub("genome",chrom,resultsGENO$Pvalset),"-withGIANT-pvaluesUncapped-scatter.png" ),width=1000,height=1000)

x = -log10(giantGeno$P)
y = -log10(DF1$P)

plot(x[giantGeno$CHR==2],y[giantGeno$CHR==2],ylab="-log10(p) UK Biobank genotype data",xlab="-log10(p) GIANT data",main=paste0(sum(!is.na(giantGeno$SNP[giantGeno$CHR==2]))," overlapping markers."),pch=16,col=add.alpha(colorsGeno[giantGeno$CHR==2],0.5))

# phased
points(x[(inPhasedGeno)&(giantGeno$CHR==2)],y[(inPhasedGeno)&(giantGeno$CHR==2)],col=phasedColour,pch=16)
# rare
points(x[(rareGeno)&(giantGeno$CHR==2)],y[(rareGeno)&(giantGeno$CHR==2)],col=rareColour,pch=1)
# low info
points(x[(poorInfoGeno)&(giantGeno$CHR==2)],y[(poorInfoGeno)&(giantGeno$CHR==2)],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

abline(0,1,col="red",lty=3)
dev.off()


# Imputed data
png(paste0( "plots/",resultsIMP$Pvalset,"-withGIANT-pvalues-scatter.png" ),width=1000,height=1000)
x = -log10(giantImp$P2)
y = -log10(DF2$P2)
plot(x,y,ylab="-log10(p) UK Biobank imputed data",xlab="-log10(p) GIANT data",xlim=c(0,Ymax),ylim=c(0,Ymax),main=paste0(sum(!is.na(giantImp$SNP))," overlapping markers."),pch=16,col="black")

# phased
points(x[inPhasedImp],y[inPhasedImp],col=phasedColour,pch=16)
# rare
points(x[rareImp],y[rareImp],col=rareColour,pch=1)
# low info
points(x[poorInfoImp],y[poorInfoImp],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)
abline(h=-log10(5e-8),col=add.alpha("lightgray",0.8),lty=3)
abline(v=-log10(5e-8),col=add.alpha("lightgray",0.8),lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("INFO < 0.3")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

abline(0,1,col="red",lty=3)
dev.off()

png(paste0( "plots/",resultsIMP$Pvalset,"-withGIANT-pvaluesUncapped-scatter.png" ),width=1000,height=1000)

x = -log10(giantImp$P)
y = -log10(DF2$P)

plot(x,y,ylab="-log10(p) UK Biobank imputed data",xlab="-log10(p) GIANT data",main=paste0(sum(!is.na(giantImp$SNP))," overlapping markers."),pch=16,col="black")

# phased
points(x[inPhasedImp],y[inPhasedImp],col=phasedColour,pch=16)
# rare
points(x[rareImp],y[rareImp],col=rareColour,pch=1)
# low info
points(x[poorInfoImp],y[poorInfoImp],col=infoColour,pch=1)

abline(0,1,col="red",lty=3)
abline(h=-log10(5e-8),col=add.alpha("lightgray",0.8),lty=3)
abline(v=-log10(5e-8),col=add.alpha("lightgray",0.8),lty=3)

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (UK Biobank) < ",rareThreshold),paste0("INFO < 0.3")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))

dev.off()



##### where Giant does better

to=which((-log10(giantImp$P)>(-log10(5e-8)))&(-log10(DF2$P)<(-log10(1e-5)))&(DF2$CHR==2))
length(to)
toCheckGiant=giantImp[to,]
toCheckImp=DF2[to,]
                                        # exclude anything with poor info or major differences in allele freq
toExclude = ( abs(toCheckGiant$MAF-toCheckImp$MAF)>0.1 )|( toCheckImp$INFO<0.3 )|( toCheckImp$MAF < rareThreshold )|( !toCheckImp$alternate_ids %in% DF1$SNP )

toCheckImp2 = toCheckImp[!toExclude,]
toCheckGiant2 = toCheckGiant[!toExclude,]

for(r in 1:nrow(toCheckImp2)){
    reg = c(toCheckImp2$BP[r]-0.5e6,toCheckImp2$BP[r]+0.5e6)
    
    print(r)
    
    GenoSNPsHere = resultsGENO$DF[(resultsGENO$DF$CHR==chrom)&(resultsGENO$DF$BP>=reg[1])&(resultsGENO$DF$BP<=reg[2]),]
    if(dim(GenoSNPsHere)[1]==0) {
        print()
        next
    }

    giantCols = rep(add.alpha(giantColour,0.6),sum(giantData$DF$CHR==2))
    giantCols[giantData$DF[giantData$DF$CHR==2,"SNP"]==toCheckGiant2$SNP[r]]="green"
    
    plot.gwas.region.v2(resultsGENO$DF,GenoSNPsHere$SNP[GenoSNPsHere$P==min(GenoSNPsHere$P)],catFile,ld=NULL,ldIMP=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=reg,Ymax=Ymax,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(resultsIMP$Pvalset,"-forASHG2.withGIANTbetter-",toCheckImp2$SNP[r]),addGenes=list("genes"=genes,"exons"=exons),thisChrom=chrom,
                        ldIMPvariant = toCheckImp2$SNP[r],
                        showRegion=NULL,
                        basicCol=add.alpha("black",0.8),
                        imputedSNPs=resultsIMP$DF,imputeCol="lightgray",extraPoints=giantData$DF[giantData$DF$CHR==2,],extraPointsCol=giantCols)
    dev.off()
    
}



##### Highlight the lactase snp.

lactasePos = 136608650
lac = (DF2$BP==lactasePos)&(DF2$CHR==2)

# the nearest in the genotype data is:
lactasePos2 = 136608646
lac2 = (DF1$BP==lactasePos2)&(DF1$CHR==2)

png(paste0( "plots/",resultsGENO$Pvalset,"-with-Place.of.birth.in.UK...north.co.ordinate-pvaluesUncapped-scatter.png" ),width=1000,height=1000)
toPlot=(northGeno$P<5e-8)|(DF1$P<5e-8)
x = -log10(northGeno$P)[toPlot]
y = -log10(DF1$P)[toPlot]
plot(x,y,ylab="-log10(p) Genotype data Standing Height",xlab="-log10(p) Genotype data North Coordinate",main=paste0(sum(!is.na(northGeno$SNP))," overlapping markers."),pch=16,col=add.alpha("black",0.5))

abline(0,1,col="red",lty=3)
points(-log10(northGeno$P)[lac2],-log10(DF1$P)[lac2],col="red",pch=16)
text(-log10(northGeno$P)[lac2],-log10(DF1$P)[lac2],paste0(DF1$SNP[lac2]," LCT"),adj=c(0.01,0.01))
dev.off()


png(paste0( "plots/",resultsIMP$Pvalset,"-with-Place.of.birth.in.UK...north.co.ordinate-pvaluesUncapped-scatter.png" ),width=1000,height=1000)
toPlot=(northImp$P<5e-8)|(DF2$P<5e-8)
x = -log10(northImp$P)[toPlot]
y = -log10(DF2$P)[toPlot]
plot(x,y,ylab="-log10(p) Imputed data Standing Height",xlab="-log10(p) Imputed data North Coordinate",main=paste0(sum(!is.na(northImp$SNP))," overlapping markers."),pch=16,col=add.alpha("black",0.5))

abline(0,1,col="red",lty=3)
points(-log10(northImp$P)[lac],-log10(DF2$P)[lac],col="red",pch=16)
text(-log10(northImp$P)[lac],-log10(DF2$P)[lac],paste0(DF2$SNP[lac]," LCT"),adj=c(0.01,0.01))

dev.off()



###################################
# What happens at LACTASE? ===> no association!!

pos = 136608650

pos = 136616754
IMPvariant = (DF2$BP==pos)&(DF2$CHR==2)


###################################
# Any differences between v14 and v16? Only 33 inds removed.

resultsGENOv16 = read.gwas(gsub("v14","v16",dataFile),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)

png(paste0( "plots/",resultsGENO$Pvalset,"-raw-compared-with-v16.png" ),width=1000,height=1000)

x = -log10(resultsGENO$DFraw$P_LINREG)[resultsGENO$DFraw$SNP%in%resultsGENOv16$DFraw$SNP]
y = -log10(resultsGENOv16$DFraw$P_LINREG)
plot(x,y,ylab="-log10(p) Genotypes v14",xlab="-log10(p) Genotypes v16")
abline(0,1,col="red",lty=3)

dev.off()


###################################
# Any differences between v15 and v13? Only 33 inds removed.

#resultsIMPv15 = read.gwas(gsub("v13","v15",bgenData),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)

resultsIMPv13 = read.gwas(gsub("v15","v13",bgenData),chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)

#resultsIMPv15=resultsIMP
overlap = intersect(resultsIMPv13$DF$SNP2,resultsIMPv15$DF$SNP2)

snpOrder13 = match(overlap,resultsIMPv13$DF$SNP2)
snpOrder15 = match(overlap,resultsIMPv15$DF$SNP2)

png(paste0( "plots/",resultsIMPv15$Pvalset,"-compared-with-v13.png" ),width=1000,height=1000,res=150)

y = -log10(resultsIMPv13$DF$P2)[snpOrder13]
x = -log10(resultsIMPv15$DF$P2)[snpOrder15]

plot(x,y,ylab="-log10(p) Imputed genotypes v13",xlab="-log10(p) Imputed genotypes v15",col=colorsImp[snpOrder15],pch=16,cex=0.5,main=paste0(length(y)," markers"))
# phased
points(x[colorsImp[snpOrder15]==phasedColour],y[colorsImp[snpOrder15]==phasedColour],col=phasedColour,pch=16,cex=0.5)
# rare
points(x[resultsIMPv15$DF$MAF[snpOrder15]<rareThreshold],y[resultsIMPv15$DF$MAF[snpOrder15]<rareThreshold],col=rareColour,pch=1,cex=0.5)

abline(0,1,col="red",lty=1,main=paste0(length(overlap)," overlapping snps (by chr:pos:alleles)"),main=paste0(length(overlap)," overlapping snps (by chr:pos:alleles)"))

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (v15) < ",rareThreshold)),col=c("black",phasedColour,rareColour),bty="n",pch=c(16,16,1),cex=1,pt.lwd=c(1,1,2))       

dev.off()


png(paste0( "plots/",resultsIMPv15$Pvalset,"-compared-with-v13-uncapped.png" ),width=1000,height=1000)

y = -log10(resultsIMPv13$DF$P)[snpOrder13]
x = -log10(resultsIMPv15$DF$P)[snpOrder15]

plot(x,y,ylab="-log10(p) Imputed genotypes v13",xlab="-log10(p) Imputed genotypes v15",col=colorsImp[snpOrder15],pch=16,cex=1,main=paste0(length(y)," markers"))
# phased
points(x[colorsImp[snpOrder15]==phasedColour],y[colorsImp[snpOrder15]==phasedColour],col=phasedColour,pch=16,cex=1)
# rare
points(x[resultsIMPv15$DF$MAF[snpOrder15]<rareThreshold],y[resultsIMPv15$DF$MAF[snpOrder15]<rareThreshold],col=rareColour,pch=1,cex=1)

abline(0,1,col="red",lty=1,main=paste0(length(overlap)," overlapping snps (by chr:pos:alleles)"),main=paste0(length(overlap)," overlapping snps (by chr:pos:alleles)"))

legend("topleft",legend=c("Not phased","Phased",paste0("MAF (v15) < ",rareThreshold)),col=c("black",phasedColour,rareColour),bty="n",pch=c(16,16,1),cex=1,pt.lwd=c(1,1,2))       

dev.off()

# any big differences?

check15=resultsIMPv15$DF[snpOrder15,][abs(y-x)>1,]
check13=resultsIMPv13$DF[snpOrder13,][abs(y-x)>1,]


