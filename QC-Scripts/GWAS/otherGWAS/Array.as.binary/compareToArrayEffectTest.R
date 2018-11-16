h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")
BB.ps2snp.sex = ukbiobank.ps2snp("sexchrom")
BL.ps2snp.sex = ukbileve.ps2snp("sexchrom")
ps2snp = unique(rbind(BB.ps2snp,BL.ps2snp,BB.ps2snp.sex,BL.ps2snp.sex ))

QCexclude=c()

args = commandArgs(trailingOnly=TRUE)

args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Array.as.binary/BOLTLMM.v16/Array.as.binary-BOLT-LMM-v16.out","/well/ukbiobank/expt/V2_QCed.SNP-QC/data/Aggregate_Batches/V2_QCed.UKBiLEVEAX_b1-b11.Batch_b001-b095.autosome.array_test.txt","all","plots","-ymax","50","-lreg")


print(args)
print(getwd())

dataFile = args[1]
chroms = args[3]
plotOutDir = args[4]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

if(chroms%in%c("all","genome")) {
    chrom="genome"
    chroms=1:22
} else {
    chroms = parse.range.string(chroms)
}

if("-lreg" %in% args ) {
    useLmmInf1 = useLmmInf2 = FALSE
} else {
    useLmmInf1 = useLmmInf2 = TRUE
    if("-lreg1" %in% args ) useLmmInf1 = FALSE
    if("-lreg2" %in% args ) useLmmInf2 = FALSE
}

# Array test results
arrayTestFile=args[2]


source("/well/donnelly/ukbiobank_project_8874/clare/sex/scripts/commonScripts/plotGenesFunction.R")
genes = load.genes(refGeneFile,fieldNames=fieldNames)



### read in recombination rates
print('Reading recombination rate files from /well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr*')

if(args[3]%in%c("all","genome")) recomChroms=c(1:22,"X_nonPAR","X_PAR1","X_PAR2") else recomChroms=chroms

recombrates = list()
for(i in recomChroms){    
    print(i)
    r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
    recombrates[[i]] = r
}


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

resultsGENO = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf1)
DF1 = resultsGENO$DF


###################################
# read array effect results!

                                        # frequency calculations (from compute-frequencies.sh)
snpFrequencyData = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr"


DFraw = read.table(arrayTestFile,sep="",header=TRUE,stringsAsFactors=FALSE)

minP = unlist(apply(DFraw[,-1],1,function(x) x[x==min(x,na.rm=TRUE)][1]))

DFraw$P=minP    # This is the results from the minimum of the three p-values

if(!Ymax) Ymax = ceiling(max(-log10(DFraw$P[DFraw$P!=0]),na.rm=T)) + 10

DFraw$P2=DFraw$P
DFraw$P2[DFraw$P<(10^-Ymax)] = 10^-Ymax
DFraw$P_LINREG=DFraw$P

index = match(DFraw$AffySNPID,ps2snp$AffySNPID)

DFraw$SNP=ps2snp$dbSNPRSID[index]
DFraw$BP=ps2snp$Position[index]
DFraw$CHR=ps2snp$Chromosome[index]
DFraw$ALLELE0 = ps2snp$AlleleA[index]
DFraw$ALLELE1 = ps2snp$AlleleB[index]
DFraw$SNP2  = altSNPID(DFraw$CHR,DFraw$BP,DFraw$ALLELE0,DFraw$ALLELE1)
    
freqs = read.genotyped.maf(mafBins = c(0,1/1000,1/100,5/100,Inf))
freqs$AffySNPID =ps2snp$AffySNPID[match(freqs$SNP,ps2snp$dbSNPRSID)]
freqs$F_MISS = freqs$C.MISSING./(freqs$C.MISSING.+rowSums(freqs[,5:9]))

indexFreqs = match(DFraw$AffySNPID,freqs$AffySNPID)
# replace any missing miss or freq values with 0 (this is artificial anyway)
DFraw$ingenotypes=1
DFraw$ingenotypes[is.na(indexFreqs)]=0

DFraw$F_MISS = freqs$F_MISS[indexFreqs]
DFraw$MAF = freqs$MAF[indexFreqs]
DFraw$MAF[is.na(indexFreqs)]=0
DFraw$F_MISS[is.na(indexFreqs)]=0
DFraw$A1FREQ = DFraw$MAF # NOTE: this is not necessarily true! Just there so read.gwas function works.

DF = dplyr::tbl_df(DFraw)

resultsTEST = list("DFraw"=DFraw,"DF"=DF,"Pvalset"=gsub(".txt","-chrgenome",basename(arrayTestFile)))
DF2 = resultsTEST$DF

dummyResultsFile = paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/Array.as.binary/",gsub(".txt",".out",basename(arrayTestFile)))

write.table(DF2,file=dummyResultsFile,quote=FALSE,col.names=TRUE,row.names=FALSE)




###################################
# index for combined data
index2a = match(DF2$SNP2,DF1$SNP2)
index2b = match(altSNPID(DFraw$CHR,DFraw$BP,DFraw$ALLELE1,DFraw$ALLELE0),DF1$SNP2)
index2 =  index2a
index2[is.na(index2a)] = index2b[is.na(index2a)]
    
DF1a = DF1[index2,]


###################################
# Filename

vers1 = gsub("\\.out","",str_split(dataFile,"-v")[[1]][2])
vers2 = gsub("\\.txt","",basename(arrayTestFile))
    
pvalset1=resultsGENO$Pvalset

filename=paste0(pvalset1,"-with-",vers2)
if(!useLmmInf1) filename=paste0(pvalset1,"-lreg-with-",vers2)


###################################
# Get  colours for plots

phasedColour = "red"
rareColour = "blue"
infoColour = "orange"
rareThreshold = 0.001
excludedCol = "darkgreen"


imputedData=FALSE
if( grepl("chr", dataFile) ){ # are we comparing results from imputation?
    print("Data are imputed.")
    imputedData=TRUE
    chrom = as.numeric( gsub("\\.out","",str_split(basename(dataFile1),"chr")[[1]][2]) )
    phased = read.phased.snps(chrom)[[1]]
                                        # all based on DF1 values
    inPhased = DF1a$alternate_ids%in%phased$V2
    poorInfo = DF1a$INFO < 0.3
    
} else {
    
    p = read.phased.snps(chroms)
    
    phased = rbind_all(p)
    
    inPhased = DF1a$SNP%in%phased$V2
    
    poorInfo = DF1a$F_MISS >= 0.05

    excluded = is.na(DF1a$SNP)

}

rare = DF1a$MAF < rareThreshold

colorsGeno = rep("black",dim(DF1a)[1]) 
colorsGeno[inPhased]=phasedColour # in phased
colorsGeno[excluded]=excludedCol # excluded from released genotypes!



###################################
# Plot manhattans on top of each other

for(chrom in chroms){
    
    print(chrom)
    # Array test
    plot.BOLT.pvalues(resultsTEST,chrom,plotOutDir,plotQQ=TRUE,catFile = NULL,Ymax=Ymax,extraData=NULL,extraPointsCol="lightgray",moreToPlot=TRUE,extraFilename=paste0(".colourPhasedandRareLT",rareThreshold),qqylim=NULL)

    inchrom=DF2$CHR==chrom
    points(DF2$BP[(inPhased)&(inchrom)],-log10(DF2$P2[(inPhased)&(inchrom)]),col=phasedColour,pch=16)
    points(DF2$BP[(excluded)&(inchrom)],-log10(DF2$P2[(excluded)&(inchrom)]),col=excludedCol,pch=16)
    points(DF2$BP[(rare)&(inchrom)],-log10(DF2$P2[(rare)&(inchrom)]),col=rareColour,pch=1,lwd=2)    
    points(DF2$BP[(poorInfo)&(inchrom)],-log10(DF2$P2[(poorInfo)&(inchrom)]),col=infoColour,pch=1,lwd=2)

    if(!imputedData) legend("topleft",legend=c("Not phased","Phased",paste0("MAF v",vers1," < ",rareThreshold),paste0("MISS >= 0.05"),"excluded from genotypes"),col=c("black",phasedColour,rareColour,infoColour,excludedCol),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2))  

    dev.off()
    

}

############################################
# combine with GWAS results into pdf


for(chrom in chroms){

    print("combining plots into one pdf")
    
    impplot = list.files(path="./plots",pattern=paste0(gsub("genome",chrom,resultsTEST$Pvalset),"\\.*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

    genoplot = list.files(path="./plots",pattern=paste0(gsub("genome",chrom,resultsGENO$Pvalset),".*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

    filename = paste0("plots/", str_split(genoplot,".pruned")[[1]][1],"_COMBINED_",impplot)
    filename=gsub(".png",".pdf",filename)

                                        # IOP makes filename too long!!
    if(grepl("Intra.ocular.pressure",filename)) filename=gsub("Intra.ocular.pressure","IOP",filename)

    print(impplot)
    print(genoplot)
    
    pdf(filename,width=61,height=24)
    
    par(mar=c(0,0,0,0))
    plot.new()
    
    img1 <- readPNG(paste0("plots/",impplot))
    img2 <- readPNG(paste0("plots/",genoplot))
                                        # plot top row, then 2 x 1 matrix below
    rasterImage(img1,0,1/2,1,1)
    rasterImage(img2,0,0,1,1/2)
    
    dev.off()
    
}

############################################
# Plot specific regions!!!


myRegions = c("15 7.88e7 7.895e7","15 78633683 7.88e7")
    
write.table(myRegions,file="regions.to.check.array2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


forSystem = paste0("Rscript ../get.hit.regions.R ",dummyResultsFile," 15 plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw.for.array.test -regions regions.to.check.array2.txt -fixRegions")
system(forSystem)

r=read.gwas(dummyResultsFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)

hitsDataFile = paste0(basename(dummyResultsFile),".chrgenome.maf0.miss1.pruned-raw.for.array.test-lreg-hitsAndLDcalc.RData")
load(hitsDataFile,verbose=TRUE)
HITSchromGENO = HITSchrom
LDchromGENO = LDchrom
    
for (chrom in chroms){

    print(chrom)
                                        # get exon information
    inchrom=DF2$CHR==chrom
    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    HITS = HITSchromGENO[[as.character(chrom)]]
    ld = LDchromGENO[[as.character(chrom)]]
    recomb = recombrates[[as.character(chrom)]]

    print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )

#### plot the regions
    ldIMP=NULL
    if( "LDchromIMP"%in%ls() ){
        ldIMP = LDchromIMP[[as.character(chrom)]]
        if(!is.null(ldIMP)) ldIMP = as.data.frame(ldIMP)
    }
    
    if(length(HITS$topVariants) > 0){

        exons = get.exons(genes,chromosome = chrom)
            
        sapply(1:length(HITS$topVariants),FUN=function(s) {
                                        # s = 17;
                                        #for(s in h) source(s)
                                        #      s=49;
                                        # s=51 (chromosome 2, this is the ACP1 region)
            
            regionToPlotWide = c(HITS$regs[s,1],HITS$regs[s,2])

############################################
            regionToPlot = c(78633683,79061002) # this is what I told Colin to exclude.
############################################
            
            regionToPlotWide = c(regionToPlot[1]-1e6,regionToPlot[2]+1e6)
            regionToPlotWide = c(77874000,79868000) # includes the range of ld calculation for the main snp.

            Legend = list("topright","legend"=c("In genotype data release","Excluded from genotype data release"),"pch"=c(23,23),"col"=c("darkgray",excludedCol),"pt.bg"=c("darkgray",excludedCol),"bty"="n","cex"=2,"pt.lwd"=c(1,1))
            
            plot.gwas.region.v2(DF2,HITS$topVariants[s],catFile,ld=ld,ldIMP=NULL,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlotWide,Ymax=Ymax+5,Ymin=0,width=NULL,moreToPlot=FALSE,Pvalset= gsub("genome",chrom,resultsTEST$Pvalset),printTopSNP=TRUE,
                                showRegion=regionToPlot,addGenes=list("exons"=exons,"genes"=genes),
                                basicCol="black",Legend=Legend,main="",
                                extraPoints=DF2[(inchrom)&(excluded),],extraPointsCol=excludedCol)


            
        })
        
        
    }
    
}

inRegion=(DF2$BP>regionToPlot[1])&(DF2$BP<regionToPlot[2])&(inchrom)


############
# Cluster plots for snps
############

forSystem = paste0( "Rscript ../cluster-region-plots.R Array.as.binary -vers ",vers2," -outdir clusterplots_snpsToCheck -snps ",ps2snp$AffySNPID[match(HITSchromGENO[["15"]]$topVariants,ps2snp$dbSNPRSID)]," -rsids -prefix ",vers2,".\\*")
system(forSystem)


############################################
# Plot the p-values where they overlap

x =  DF1a$P2 # GWAS
y =  DF2$P2 # array test

plotDir = plotOutDir

excl = (x>10^-3)&(y>10^-3)


png(paste0(plotDir,"/",filename,"-scatter-pvals.png"),width=1000,height=1000,res=150)      

plot(-log10(x)[!excl],-log10(y)[!excl],xlab=paste0("-log10(p) version v",vers1),ylab=paste0("-log10(p) version v",vers2),pch=16,main=paste0(sum(!is.na(y))," overlapping snps."),col="black",xlim=c(0,Ymax),ylim=c(0,Ymax),cex=0.8)
# phased
points(-log10(x)[(!excl)&(inPhased)],-log10(y)[(!excl)&(inPhased)],col=phasedColour,pch=16,cex=0.8)
# rare
points(-log10(x)[(!excl)&(rare)],-log10(y)[(!excl)&(rare)],col=rareColour,pch=1,cex=0.8)
# low info
points(-log10(x)[(!excl)&(poorInfo)],-log10(y)[(!excl)&(poorInfo)],col=infoColour,pch=1,cex=0.8)

abline(0,1,col="red",lty=3)

if(imputedData) legend("topleft",legend=c("Not phased","Phased",paste0("MAF ",vers1," < ",rareThreshold),paste0("INFO ",vers1," < 0.3")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2))       

if(!imputedData) legend("topleft",legend=c("Not phased","Phased",paste0("MAF ",vers1," < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2))       

dev.off()


png(paste0(plotDir,"/",filename,"-scatter-pvalsUncapped.png"),width=1000,height=1000,res=150)

x =  DF1a$P # GWAS
y =  DF2$P # array test


plot(-log10(x)[!excl],-log10(y)[!excl],xlab=paste0("-log10(p) version v",vers1),ylab=paste0("-log10(p) version v",vers2),pch=16,main=paste0(sum(!is.na(y))," overlapping snps."),col="black",cex=0.8)
# phased
points(-log10(x)[(!excl)&(inPhased)],-log10(y)[(!excl)&(inPhased)],col=phasedColour,pch=16,cex=0.8)
# rare
points(-log10(x)[(!excl)&(rare)],-log10(y)[(!excl)&(rare)],col=rareColour,pch=1,cex=0.8)
# low info
points(-log10(x)[(!excl)&(poorInfo)],-log10(y)[(!excl)&(poorInfo)],col=infoColour,pch=1,cex=0.8)

abline(0,1,col="red",lty=3)

if(imputedData) legend("topleft",legend=c("Not phased","Phased",paste0("MAF ",vers1," < ",rareThreshold),paste0("INFO ",vers1," < 0.3")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

if(!imputedData) legend("topleft",legend=c("Not phased","Phased",paste0("MAF ",vers1," < ",rareThreshold),paste0("MISS >= 0.05")),col=c("black",phasedColour,rareColour,infoColour),bty="n",pch=c(16,16,1,1),cex=1,pt.lwd=c(1,1,2,2))       

dev.off()
