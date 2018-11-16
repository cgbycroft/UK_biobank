source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

library(dplyr)
library(qqman)
library(stringr)

sexChroms = c(23,24,25,26)
names(sexChroms) = c("X","Y","XY","MT")

plot.BOLT.pvalues <- function(GWASdata,chrom,minmaf,mininfo,maxmiss,plotOutDir,plotQQ=TRUE,extraTitle="",Ymax=FALSE,catFile = NULL,QCexclude=c(),...) {

    DF = dplyr::tbl_df(DFraw)
    print(head(DF))
    if("INFO"%in%colnames(DFraw)){
        Pvalset = paste(basename(GWASdata),".chr",chrom,".maf",minmaf,".info",mininfo,".pruned",extraTitle,sep="")
        DF = dplyr::filter(DF, MAF > minmaf & INFO > mininfo)
        
    } else {
        Pvalset = paste(basename(GWASdata),".chr",chrom,".maf",minmaf,".miss",maxmiss,".pruned",extraTitle,sep="")
        DF = dplyr::filter(DF, MAF > minmaf & F_MISS < maxmiss)
    }

    if(chrom!="genome") DF = dplyr::filter(DF, CHR %in% chrom)  # subset by chromosome

    if( "-qc" %in%args ) DF = dplyr::filter(DF, !SNP %in% QCexclude)  # exclude SNPs in array,imageArtefact, or concordance lists. Only relevant with clare's plink gwas data.

    
    print(Pvalset)
    
    
    if("P_BOLT_LMM_INF" %in% colnames(DF)) {
        print("using BOLT_LMM_INF")
        DF = dplyr::rename(DF, P = P_BOLT_LMM_INF) } else {
        DF = dplyr::rename(DF, P = P_LINREG)
    }
    
    if(chrom%in%c("X","XY","Y","MT")) DF$CHR = sexChroms[chrom] else DF$CHR = as.numeric(DF$CHR)

    maxP = round(max(-log10(DF$P),na.rm=T))
    print(paste('max -log(pval) = ',maxP))

    nFail = length(which(-log10(DF$P) > 8))
    percentFailed = nFail/nrow(DF) * 100

    if(!Ymax){

        Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10
        Ymax = min(Ymax,50)

    }
    
    png(paste(plotOutDir,"/",Pvalset,"-manhattan%02d.png",sep=""),width=41,height=12,units="in",res=150)
    par(las=1,font.main=1,cex.axis=2,cex.lab=2,mar=c(7 ,7, 5 ,2))
    myManhattan(DF,ymax=Ymax,suggestiveline = FALSE,xpd=NA,cex=1,...)
#    myManhattan(DF,ymax=Ymax,suggestiveline = FALSE,xpd=NA,cex=1,col="transparent")

    # add extra catalogue hits?
    if(!is.null(catFile) & (chrom!="genome")) {
        print( "Printing catalogue hits..." )

        catFileSub = catFile[catFile$CHR %in% chrom,]
        print( paste0( sum(-log10(catFileSub$Pvalue)>Ymax)," catalogue hits above ",Ymax) )
        
        catFileSub$Pvalue[-log10(catFileSub$Pvalue)>Ymax] = 10^(-Ymax)
        # plot non-European hits differently
        colors = rep("red",dim(catFileSub)[1])
#        colors[!grepl("European",catFileSub$Ancestry)] = "blue"
        print( table(catFileSub$Ancestry[!grepl("European",catFileSub$Ancestry)]) )

        # do we have these SNPs in UKBiobank? match on chrom and position
        #inHere = (catFileSub$BP %in% DF$BP)&(catFileSub$CHR == DF$CHR)
        #catFileSub$Pvalue[inHere] = DF$Pvalue[]
        
        points(catFileSub$BP,-log10(catFileSub$Pvalue),pch=8,col=colors,cex=4,lwd=1)
        points(catFileSub$BP,-log10(catFileSub$Pvalue),pch=16,col=colors,cex=2)
        
    }
    
    dev.off()

    
########### qq plot p-values
    if(plotQQ){
        png(paste(plotOutDir,"/",Pvalset,"-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
        DF$P2=DF$P
        DF$P2[DF$P<(10^-Ymax)] = 10^-Ymax
        qqman::qq(DF$P2)
        dev.off()
    }

########## plot effect sizes
    DF$index = DF$BP

    if( length(unique(DF$CHR)) > 1 ){
        for(i in unique(DF$CHR)){        
            if(i>1) DF$index[DF$CHR==i] = DF$index[DF$CHR==i] + max(DF$BP[DF$CHR==(i - 1)])
        }
    }        
    snps = which((DF$P < 5e-8)&(!is.na(DF$P)))
    if("BETA"%in%colnames(DF) & ( length(snps) > 0 )){
        
        beta = DF
        beta$BETA[DF$A1FREQ > 0.5] = -beta$BETA[DF$A1FREQ > 0.5]
        beta = beta[snps,]

        png(paste(plotOutDir,"/",Pvalset,"-EffectSizes.png",sep=""),height=1000,width=1000,res=150)
        myManhattan(beta,p="BETA",logtransform=FALSE,genomewideline=0,suggestiveline=FALSE)
        
        dev.off()
        
    }
}




######## Which chromosome and what data are we plotting?

args = commandArgs(TRUE)

#args = c("test.out","1","plots", "-ymax","50","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-Euro-hits")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v3/Standing.height-BOLT-LMM-v3.out","5","plots", "-ymax","50","-qc","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-Euro-hits-for-talk")

print(args)

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""


# Should we highlight some SNPs?
highlightSNPs=NULL
highlightCols=NULL


if(( "-hi"%in%args )|("-qc" %in%args )){
    
    print("Reading QC snps lists...")
    QCSNPList = read.SNPQC.files(justSNPs=TRUE)
    
#    QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$imageSNPs,QCSNPList$concordanceSNPs))    
     QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$concordanceSNPs))    # this is only required if using -hi or clare's versions of plink genotype files.
    
    if( "-hi"%in%args ){
        highlightSNPs = unique(unlist(QCSNPList))
        colors = rep("black",length(highlightSNPs))

        colors[highlightSNPs%in%QCSNPList$batchHweSNPs] = "green"   # HWE (apply first)
        colors[highlightSNPs%in%c(QCSNPList$plateSNPs,QCSNPList$batchSNPs)] = "purple"  # BATCH/PLATE
        colors[highlightSNPs%in%c(QCSNPList$imageSNPs)] = "orange" # IMAGE ARTEFACT
        colors[highlightSNPs%in%c(QCSNPList$arraySNPs)] = "red" # ARRAY
        colors[highlightSNPs%in%c(QCSNPList$concordanceSNPs)] = "blue" # CONCORDANCE
        highlightCols = colors
        
        print(table(highlightCols))
    }
        
}


# get data

print("printing the following chromosomes")
print( chroms )


#DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
if(!grepl("%%",dataFile)){
    print("reading in GWAS output file")
    DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
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

# which chroms?
if(chroms!="genome") {
    if(chroms=="all") chroms = 1:22 else chroms = parse.range.string(chroms)
}

for(chrom in chroms){
    
    if(grepl("%%",dataFile)) {
        DFraw = read.table(gsub("%%",chrom,dataFile),sep="",header=TRUE,stringsAsFactors=FALSE)
        DFraw$MAF = DFraw$A1FREQ
        DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
    }

    print(chrom)
    
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
    
    plot.BOLT.pvalues(GWASdata=gsub("%%",chrom,dataFile),chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss=maxmiss,plotOutDir=plotOutDir,extraTitle=extraTitle,highlight=highlightSNPs,highlightCols=highlightCols,Ymax=Ymax,catFile=catFile)
    
}

############# EXTRAS #############

# snps
#hiEF = catFile[(catFile$V20>5)&(!is.na(catFile$V20)),c("V5","CHR","BP","Pvalue","V20","V21")]
#DFraw$SNPID = paste0(DFraw$CHR,".",DFraw$BP)
#hiEF$SNPID = paste0(hiEF$CHR,".",hiEF$BP)

#snps = DFraw[match(hiEF$SNPID,DFraw$SNPID),]
#write.table(snps$SNP,file="Height.hi.OR.snps.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#system("plink --bfile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc --keep-allele-order --extract Height.hi.OR.snps.txt --recode AD --out Height.hi.OR.snps")

#geno = read.table("Height.hi.OR.snps.raw",header=TRUE)
#map = read.table("Height.hi.OR.snps.map",header=FALSE)
#pheno = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt",header=TRUE)

### get height phenotype
#geno$pheno = pheno$Standing.height[match(geno$IID,pheno$IID)]

#snp = "Affx-20256845"
#snp1 = 5 + 2*which(map$V2==snp)
#png(paste0("Standing.height-effect-",snp,".png"),height=1000,width=1000,res=150)
#boxplot(log(geno$pheno)~geno[,snp1])
#dev.off()

#lm(geno$pheno~geno[,snp1])

############# TESTING & USAGE #############

#Rscript plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-quant-Age-v1.out all plots > Logs/plot-Ychrom-BOLT-LMM-quant-Age-v1.log &

#Rscript plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LREG-quant-Age-v1.out all plots > Logs/plot-Ychrom-BOLT-LREG-quant-Age-v1.log &

#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-quant-Age-v1.out all plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-Age-v1.log &

    
#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-all-snps-quant-Age-v1.out all plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-all-snps-Age-all-v1.log &

#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-all-snps-quant-Age-v1.out genome plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-all-snps-Age-all-v1.log &


# dataFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v1/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-all-snps-v1.out"
