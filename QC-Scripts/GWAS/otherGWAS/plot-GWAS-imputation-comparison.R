h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = rbind(ukbiobank.ps2snp("autosome"),ukbiobank.ps2snp("sexchrom"))
BL.ps2snp = rbind(ukbileve.ps2snp("autosome"),ukbileve.ps2snp("sexchrom"))
QCexclude=c()

args = commandArgs(trailingOnly=TRUE)

                                        #args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v14/Standing.height-BOLT-LMM-v14.out","2","plots","-ymax","50","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v13/Standing.height-BOLT-LMM-v13-chr2.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr2.hrc+uk10k.I4.v1.1.sample","-dontPlotManhattans","-ldRData","Standing.height-BOLT-LMM-v14.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp","Standing.height-BOLT-LMM-v13-chr2.out.chr2.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v16s/Internal.Pico..ng.uL.-BOLT-LMM-v16s.out","23","plots","-ymax","50","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v15s/Internal.Pico..ng.uL.-BOLT-LMM-v15s-chr23.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr2.hrc+uk10k.I4.v1.1.sample","-dontPlotManhattans","-ldRData","Internal.Pico..ng.uL.-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp","Internal.Pico..ng.uL.-BOLT-LMM-v15s-chr2.out.chr2.maf0.001.info0.3.pruned-QCFiltered-lreg-hitsAndLDcalc.RData")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v16s/Internal.Pico..ng.uL.-BOLT-LMM-v16s.out","23","plots","-ymax","50","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgenFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v15s/Internal.Pico..ng.uL.-BOLT-LMM-v15s-chr23.out","-sample","/well/ukbiobank/imputation/final/full/bgen/chr2.hrc+uk10k.I4.v1.1.sample","-dontComputeLD")

                                        #args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s.out", "23","plots","-ymax", "50", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg", "-title", "-raw","-bgenFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15s/Standing.height-BOLT-LMM-v15s-chr23.out","-sample", "/well/ukbiobank/imputation/final/full/bgen/chrX.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans","-ldRData","Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp","Standing.height-BOLT-LMM-v15s-chr23.out.chrX.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData")

## v2 of bgen data
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16/Standing.height-BOLT-LMM-v16.out", "2","plots","-ymax", "50", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg", "-title", "-raw","-bgenFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19/Standing.height-BOLT-LMM-v19-chr2.out","-sample", "/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample","-dontComputeLD")

#PAR
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s.out" ,"25" ,"plots", "-ymax" ,"50" ,"-phenoFile" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-rawQChits.v19s", "-bgenFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr25.out" ,"-sample" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrPAR.v2.sample" ,"-ldRData" ,"Standing.height-BOLT-LMM-v16s-25.out.chr25.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp" ,"Standing.height-BOLT-LMM-v19s-chr23.out.chrX.maf0.001.info0.3.pruned-QCFiltered-lreg-ct0.1-hitsAndLDcalc.RData", "-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-dontPlotManhattans")

# just manhattans
#args=c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s.out","25" ,"plots", "-ymax" ,"50" ,"-qc","-phenoFile" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-title", "-QCFiltered", "-bgenFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr25.out", "-sample", "/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample/v2/chrPAR1.v2.sample", "-dontComputeLD", "-par", "2")

###### X chrom
# V19s
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s.out" ,"23" ,"plots", "-ymax" ,"50" ,"-phenoFile" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-rawQChits.v19s", "-bgenFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr23.out" ,"-sample" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.v2.sample" ,"-ldRData" ,"Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-ldRDataImp" ,"Standing.height-BOLT-LMM-v19s-chr23.out.chrX.maf0.001.info0.3.pruned-QCFiltered-lreg-ct0.1-hitsAndLDcalc.RData", "-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-dontPlotManhattans")


print(args)
print(getwd())

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

plotRegions=TRUE
plotManhattans=TRUE
if("-dontPlotRegions"%in%args) plotRegions=FALSE
if("-dontPlotManhattans"%in%args) plotManhattans=FALSE

chroms = parse.range.string(chroms)
if(length(chroms)==1){
    if(chroms%in%c(23:26,"X","XY","Y","MT")) dataFile = paste0(dirname(dataFile),"/",gsub(".out",paste0("-",chroms,".out"),basename(dataFile)))
}

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

# which par do we plot?
par = 0
if("-par"%in%args) par = as.numeric(args[which(args=="-par")+1])


###################################
# for each chromosome, get set of 'hits' based on recombination distances

# does this data exist already, if we haven't asked to computeLD? If we have asked to computeLD (default) then we run hits and LD functions.
# if there are no hits in this chromosome, and no manhattans are being plotted, then just quit!

hitsDataFile = paste0(basename(dataFile),".chrgenome.maf",minmaf,".miss",maxmiss,".pruned",extraTitle,"-hitsAndLDcalc.RData")
#hitsDataFile = paste0(resultsGENO$Pvalset,"-hitsAndLDcalc.RData")
if("-ldRData"%in%args) hitsDataFile = args[which(args=="-ldRData")+1]
if(hitsDataFile%in%list.files()) {
    load(hitsDataFile,verbose=TRUE)
    HITSchromGENO = HITSchrom
    LDchromGENO = LDchrom
}


if("-ldRDataImp"%in%args) hitsDataFileImputed = args[which(args=="-ldRDataImp")+1] else hitsDataFileImputed=""
if(hitsDataFileImputed%in%list.files()) {
    load(hitsDataFileImputed,verbose=TRUE)
    HITSchromIMP = HITSchrom
    LDchromIMP = LDchrom
}

#don't plot regions unless we have hits information
if( !"HITSchromGENO" %in% ls() ){
    plotRegions=FALSE
    nChromsWithRegions=0
} else {

    nChromsWithRegions = sum( unlist(sapply(chroms,function(x) !dim(HITSchromGENO[[as.character(x)]][[1]])[1]%in%c(0,NULL))) )

    print( paste0(nChromsWithRegions," chromosomes with regions to plot.") )
}

if((!plotManhattans)&( nChromsWithRegions==0 )) {
    print("No hit regions to plot for this chromosome. Stopping.")
    quit()
}


###################################
# read GWAS results!

resultsGENO = read.gwas(dataFile,chrom=chroms,minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF1 = resultsGENO$DF

resultsIMP = read.gwas(bgenData,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF2 = resultsIMP$DF


### read in recombination rates
recombrates = read.recomb(chrom=chroms)

if(par!=0){
    if(par==1) {
        print("only par1")
        resultsIMP$DF = resultsIMP$DF[resultsIMP$DF$BP <= 2799507,]
        resultsGENO$DF = resultsGENO$DF[resultsGENO$DF$BP <= 2799507,]
        resultsIMP$Pvalset = gsub("chr25.","chr25par1.",resultsIMP$Pvalset)
        resultsGENO$Pvalset = gsub("chr25.","chr25par1.",resultsGENO$Pvalset)
        DF1 = resultsGENO$DF
        DF2 = resultsIMP$DF        
    }
    if(par==2) {
        print("only par2")
        resultsIMP$DF = resultsIMP$DF[resultsIMP$DF$BP > 2799507,]
        resultsGENO$DF = resultsGENO$DF[resultsGENO$DF$BP > 2799507,]
        resultsIMP$Pvalset = gsub("chr25.","chr25par2.",resultsIMP$Pvalset)        
        resultsGENO$Pvalset = gsub("chr25.","chr25par2.",resultsGENO$Pvalset)        
        DF1 = resultsGENO$DF
        DF2 = resultsIMP$DF        
    }
}


###################################
# Read in set of snps used in phasing!
if(chroms%in%sexChroms) phaseChroms = names(sexChroms)[match(chroms,sexChroms)] else phaseChroms = chroms
phasedBims = read.phased.snps(phaseChroms)


###################################                 
# plot manhattans with various colours
phasedColour = "red"
rareColour = "blue"
rareThreshold = 0.001

if(plotManhattans){
        
    for( chrom in chroms ){

        print(chrom)
 
        phased = phasedBims[[as.character(chrom)]]
                                        #phased$SNPimp = paste0(phased$V1,":",phased$V4,"_",phased$V5,"_",phased$V6)  # reconstruct ID in imputed data. Assume that alleles are in the same order!
        
        if(chrom==25){
            if(par==1) {
                phased = phased[phased$V4 <= 2699507,]
                
            }
            if(par==2) {
                phased = phased[phased$V4 > 2799507,]
            }
        }


        print("phased SNPs not found in (potentially filtered) imputed data:")
        notInImp = phased[! ( (phased$SNP2%in%resultsIMP$DF$SNP2) | (phased$V2%in%resultsIMP$DF$SNP2) ),]
        print(nrow(notInImp))

        print("phased SNPs not found in unfiltered imputed data:")
        notInImp = phased[! ( (phased$SNP2%in%resultsIMP$DFraw$SNP2) | (phased$V2%in%resultsIMP$DFraw$SNP) ),]
        print(nrow(notInImp))
        
        print("Of them, how many have no matching base pair position on this chromosome:")
        print(sum(!    notInImp$V4%in%resultsIMP$DFraw$BP[resultsIMP$DFraw$CHR==chrom]))

        print("SNPs with bp position but alleles don't match:")
        these = notInImp[notInImp$V4%in%resultsIMP$DFraw$BP[resultsIMP$DFraw$CHR==chrom],]
        matches = resultsIMP$DFraw[resultsIMP$DFraw$CHR==chrom,][match(these$V4,resultsIMP$DFraw$BP[resultsIMP$DF$CHR==chrom]),]
        alleleMismatch = sum( (these$V5!=matches$ALLELE1)|(these$V6!=matches$ALLELE0))
        alleleMismatch == nrow(these)
        print(alleleMismatch)

                                        # writing this set to a file
        if( alleleMismatch > 0 ){
            print("printing this list of snps here:")
            print( paste0(resultsIMP$Pvalset,"-SNPsinPhased-NotInImputedRaw.txt") )
            write.table( cbind(these,matches),file=paste0(resultsIMP$Pvalset,"-SNPsinPhased-NotInImputedRaw.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE )
        }

##### 1. genotypes
        plot.BOLT.pvalues(Data=resultsGENO,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=paste0(".colourPhasedandRareLT",rareThreshold))
                                        # add phased snps (on top)

        inPhased = resultsGENO$DF$SNP%in%phased$V2
        points(resultsGENO$DF$BP[inPhased]/1000000,-log10(resultsGENO$DF$P2[inPhased]),col=phasedColour,pch=16)
                                        # rare snps (use colour of points border)
        rare = (resultsGENO$DF$MAF < rareThreshold) & (resultsGENO$DF$CHR==chrom)
        points(resultsGENO$DF$BP[rare]/1000000,-log10(resultsGENO$DF$P2[rare]),col=rareColour,pch=1,lwd=2)

                                        # add known hits (if specified)
        if(!is.null(catFile) & (chrom!="genome")) {
            print( "Printing catalogue hits..." )

            catFileSub = catFile[catFile$CHR %in% chrom,]
            print( paste0( sum(-log10(catFileSub$Pvalue)>Ymax)," catalogue hits above ",Ymax) )
            
            catFileSub$Pvalue[-log10(catFileSub$Pvalue)>Ymax] = 10^(-Ymax)
                                        # plot non-European hits differently
            colors = rep("darkgreen",dim(catFileSub)[1])
                                        #        colors[!grepl("European",catFileSub$Ancestry)] = "blue"
            print( table(catFileSub$Ancestry[!grepl("European",catFileSub$Ancestry)]) )
            
            points(catFileSub$BP/1000000,-log10(catFileSub$Pvalue),pch=8,col=colors,cex=4,lwd=1)
            points(catFileSub$BP/1000000,-log10(catFileSub$Pvalue),pch=16,col=colors,cex=2)
            
        }
        
        legend("topright",legend=c("Not phased","Phased",paste0("MAF < ",rareThreshold)),col=c("black",phasedColour,rareColour),bty="n",pch=c(16,16,1),cex=3,pt.lwd=c(1,1,2))
        
        dev.off()


##### 2. imputed
        plot.BOLT.pvalues(Data=resultsIMP,chrom,plotOutDir="plots",plotQQ=FALSE,catFile=NULL,Ymax,highlight=NULL,highlightCols=NULL,moreToPlot=TRUE,extraFilename=paste0(".colourPhasedandRareLT",rareThreshold))
                                        # add phased snps (on top). Match on position and chromosome!

        inPhased = ( (resultsIMP$DF$SNP%in%phased$SNPimp)|(resultsIMP$DF$SNP%in%phased$V2) )&( resultsIMP$DF$CHR==chrom )
        points(resultsIMP$DF$BP[inPhased]/1000000,-log10(resultsIMP$DF$P2[inPhased]),col=phasedColour,pch=16)

                                        # rare snps (use colour of points border)
        rare = (resultsIMP$DF$MAF < rareThreshold) & (resultsIMP$DF$CHR==chrom)
        points(resultsIMP$DF$BP[rare]/1000000,-log10(resultsIMP$DF$P2[rare]),col=rareColour,pch=1,lwd=2)

        legend("topright",legend=c("Not phased (imputed)","Phased",paste0("MAF < ",rareThreshold)),col=c("black",phasedColour,rareColour),bty="n",pch=c(16,16,1),cex=3,pt.lwd=c(1,1,2))
        
        dev.off()

    }


###################################                 
                                        # combine plots into one pdf

    for(chrom in chroms){

        print("combining plots into one pdf")
        
        impplot = list.files(path="./plots",pattern=paste0(resultsIMP$Pvalset,"..*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

        genoplot = list.files(path="./plots",pattern=paste0(gsub("genome",chrom,resultsGENO$Pvalset),"..*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

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

                                        # if we're plotting Standing.height then make one with GIANT in it
    if(grepl("Standing.height",resultsGENO$Pvalset)){
                                        #chroms=16:22
        
        for(chrom in chroms){

            if(chrom%in%sexChroms) next 
            
            print("combining plots into one pdf - with GIANT")
            
            impplot = list.files(path="./plots",pattern=paste0(gsub("chr.*.maf",paste0("chr",chrom,".*.maf"),resultsIMP$Pvalset),".*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

            genoplot = list.files(path="./plots",pattern=paste0(gsub("genome",chrom,resultsGENO$Pvalset),".*colourPhasedandRareLT",rareThreshold,".*manhattan.*.png"))

            giantPlot = list.files(path="./plots",pattern=paste0("GIANT.chr",chrom,".withCatalogueHits-manhattan.*.png"))
            
            filename = paste0("plots/", str_split(genoplot,".pruned")[[1]][1],"_COMBINED-withGIANT_",impplot)
            filename=gsub(".png",".pdf",filename)

            print(impplot)
            print(genoplot)
            
            pdf(filename,width=61,height=36)
            
            par(mar=c(0,0,0,0))
            plot.new()
            
            img1 <- readPNG(paste0("plots/",impplot))
            img2 <- readPNG(paste0("plots/",genoplot))
            img3 <- readPNG(paste0("plots/",giantPlot))
                                        # plot top row, then 2 x 1 matrix below
            rasterImage(img1,0,2/3,1,1)
            rasterImage(img2,0,1/3,1,2/3)
            rasterImage(img3,0,0,1,1/3)
            
            dev.off()

        }
        

    }


}


###################################                 
# plotting regions separately
# IMPORTANT NOTE: Only plot regions as determined by default qc'd version of the genotype files!
# The set of hit regions is determined by the "-ldRData" flag.


if(!plotRegions) {
    
    print("No ld-regions found for this set of results, so not plotting region plots!")
    quit()
}

                                        # get gene information
source("/well/donnelly/ukbiobank_project_8874/clare/sex/scripts/commonScripts/plotGenesFunction.R")
genes = load.genes(refGeneFile,fieldNames=fieldNames)


for (chrom in chroms){

    print(chrom)
    if(chrom%in%c(23:25)) chromName = names(sexChroms)[sexChroms==chrom][1] else chromName=chrom
    if(chrom==25) chromName="X" # This is only used for the genes extraction
    
    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    HITS = HITSchromGENO[[as.character(chrom)]]
    ld = LDchromGENO[[as.character(chrom)]]
    recomb = recombrates[[as.character(chrom)]]
    HITSimp = HITSchromIMP[[as.character(chrom)]]
    
    print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )
#### plot the regions
    ldIMP=NULL
    if( "LDchromIMP"%in%ls() ){
        ldIMP = LDchromIMP[[as.character(chrom)]]
        if(!is.null(ldIMP)) ldIMP = as.data.frame(ldIMP)
    }
    
    if(length(HITS$topVariants) > 0){
        
                                        # get exon information for this chromosome
        exons = get.exons(genes,chromosome = chromName)

        sapply(1:length(HITS$topVariants),FUN=function(s) {
                                        # s = 17;
                                        #for(s in h) source(s)
                                        #      s=49;
            # s=51 (chromosome 2, this is the ACP1 region)
            
            regionToPlot = c(HITS$regs[s,1],HITS$regs[s,2])
            
            # WITHOUT IMPUTED DATA
            plot.gwas.region.v2(DF1,topSNP=HITS$topVariants[s],catFile,ld=ld,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=Ymax,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG1"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=NULL,
                                extraPoints=DF1[(DF1$MAF < rareThreshold)&(DF1$CHR==chrom)&(DF1$BP<=regionToPlot[2])&(DF1$BP>=regionToPlot[1]),],
                                extraPointsCol="blue",
                                extraPointsPch=5)
                                        # plot blue circles around rare snps
            
            dev.off()

                                        # WITH IMPUTED DATA. Plot the LD around the top hit based on the highest pvalue.
            Legend = list("topright","legend"=c("Genotype data","Imputed data","Known hit (GWAS catalogue)"),"pch"=c(23,21,17),"col"=c("darkgray","lightgray","blue"),"pt.bg"=c("darkgray","lightgray","blue"),"bty"="n","cex"=2,"pt.lwd"=c(1,1,1))
            
            plot.gwas.region.v2(DF1,HITS$topVariants[s],catFile,ld=ld,ldIMP=ldIMP,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=NULL,width=(10^6)/3,moreToPlot=TRUE,Pvalset= paste0(Pvalset,"-forASHG2"),addGenes=list("genes"=genes,"exons"=exons),
                                showRegion=regionToPlot,
                                basicCol="black",
                                imputedSNPs=DF2,
                                #extraPoints=DF1[(DF1$MAF < rareThreshold)&(DF1$CHR==chrom),],
                                #extraPointsCol="blue",
                                        #extraPointsPch=5,
                                ldIMPvariant = HITSimp$topVariants[s],
                                Legend=Legend,axisScale=1)

            dev.off()            

        })
    }
}



print("DONE!!")

quit()




                                        #quit()
    bgenFile="/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chrX.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen"
    bgenSampleFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.sample"
    samplesFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt.samples.tmp"
    
ld = compute.ld.around.snp.bgen(variants=HITSimp$topVariants2[c(1,6)],bgenFile,bgenSampleFile,chrom,paste0(samplesFile),regions=HITSimp$regs[c(1,6),],DF2,rNumber=37614)

ld2 = compute.ld.around.snp.bgen(variants=HITSimp$topVariants2[29],bgenFile,bgenSampleFile,chrom,paste0(samplesFile),regions=HITSimp$regs[29,],DF2,rNumber=56304)

    # after setting sex into .sample file
ld3 = compute.ld.around.snp.bgen(variants=HITSimp$topVariants2[29],bgenFile,bgenSampleFile,chrom,paste0(samplesFile),regions=HITSimp$regs[29,],DF2,rNumber=60119)

ld4 = compute.ld.around.snp.bgen(variants=HITSimp$topVariants2[21],bgenFile,bgenSampleFile,chrom,paste0(samplesFile),regions=HITSimp$regs[21,],DF2,rNumber=46596)

ld5 = compute.ld.around.snp.bgen(variants=HITSimp$topVariants2[1],bgenFile,bgenSampleFile,chrom,paste0(samplesFile),regions=HITSimp$regs[1,],DF2,rNumber=80467,callthreshold=0.1,sex="F")


    
toRun=paste0(plink," --bfile genotypes.",rNumber,".tmp --keep-allele-order --keep ",samplesFile," --recode A --out genotypes.",rNumber,".tmp")

toRun=paste0(plink," --bfile genotypes.",rNumber,".tmp --keep-allele-order --keep ",samplesFile," --freq --out genotypes.",rNumber,".tmp")

toRun=paste0(plink," --bfile genotypes.",rNumber,".tmp --keep-allele-order --freq --out genotypes.",rNumber,".tmp2")
system(toRun)

raw = read.table(paste0("genotypes.",rNumber,".tmp.raw"),stringsAsFactors=FALSE,header=TRUE)
    
fr = read.table(paste0("genotypes.",rNumber,".tmp.frq"),stringsAsFactors=FALSE,header=TRUE)
fr2 = read.table(paste0("genotypes.",rNumber,".tmp2.frq"),stringsAsFactors=FALSE,header=TRUE)

ss = read.table(paste0("genotypes.",rNumber,".tmp.snp-stats"),stringsAsFactors=FALSE,header=TRUE)

sn="23:45234922_G_A"
    ldThisSNP = read.table(paste0("taggedSNPs.",rNumber,".tmp.ld"),header=TRUE)


forSystem4 = paste0( plink," --bgen ",tempBgen," --hard-call-threshold 0.1 --sample ",bgenSampleFile," --keep-allele-order --filter-females --make-bed --out genotypes.",rNumber,".tmp")
system(forSystem4)




fr2[fr2$SNP==sn,]

maf = DF2$MAF[match(fr2all$SNP,DF2$SNP2)]


######
# Check info scores for x

# compute info scores males + females separately



dataFile9="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15/Standing.height-BOLT-LMM-v15-chr9.out"

resultsIMP9 = read.gwas(dataFile9,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)

dataFile10="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15/Standing.height-BOLT-LMM-v15-chr10.out"

resultsIMP10 = read.gwas(dataFile10,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)


png("plots/info.10.hist.png",res=150,width=1000,height=1000)
hist(resultsIMP10$DFraw$INFO,breaks=100,xlim=c(0,1),main=paste0("chromosome 10: ",length(resultsIMP10$DFraw$INFO)," snps"))
dev.off()

png("plots/info.9.hist.png",res=150,width=1000,height=1000)
hist(resultsIMP9$DFraw$INFO,breaks=100,xlim=c(0,1),main=paste0("chromosome 9: ",length(resultsIMP9$DFraw$INFO)," snps"))
dev.off()

png("plots/info.X.hist.png",res=150,width=1000,height=1000)
hist(resultsIMP$DFraw$INFO,breaks=100,xlim=c(0,1),main=paste0("X chromosome: ",length(resultsIMP$DFraw$INFO)," snps"))
dev.off()

png("plots/info.X.hist2.png",res=150,width=1000,height=1000)
hist(resultsIMP$DFraw$INFO,breaks=100,main=paste0("X chromosome: ",length(resultsIMP$DFraw$INFO)," snps"))
dev.off()
