source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


plot.BOLT.pvalues <- function(GWASdata,chrom,minmaf,mininfo,maxmiss,plotOutDir,plotQQ=TRUE,extraTitle="",Ymax=FALSE,useLinReg=FALSE,...) {

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

    if( "-qc" %in%args ) DF = dplyr::filter(DF, !SNP %in% QCexclude)  # exclude SNPs in array,imageArtefact, or concordance lists

    
    print(Pvalset)
    
    
    if( ("P_BOLT_LMM_INF" %in% colnames(DF))&(!useLinReg) ) {
        print("using BOLT_LMM_INF")
        DF = dplyr::rename(DF, P = P_BOLT_LMM_INF) } else {
        print("using LINREG")                    
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
    myManhattan(DF,ymax=Ymax,...)
    dev.off()

    
########### qq plot p-values
    if(plotQQ){
        png(paste(plotOutDir,"/",Pvalset,"-qqplot%02d.png",sep=""),height=1000,width=1000,res=150)
        DF$P2=DF$P
        DF$P2[DF$P<(10^-Ymax)] = 10^-Ymax
        qqman::qq(DF$P)
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
    
    if(("BETA"%in%colnames(DF) ) & (length(snps) > 0) ){
        
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
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v11/Internal.Pico..ng.uL.-BOLT-LMM-v11-%%.out", "all", "plots", "-qc", "-ymax", "50", "-title", "-QCfiltered","-sex")

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
    
    if( "-sex" %in% args ) type = "sexchrom"  else type = "autosome"
    QCSNPList = read.SNPQC.files(justSNPs=TRUE,type=type)
    
    QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$imageSNPs,QCSNPList$concordanceSNPs))    
# NOTE: I'm excluding image SNPs as these weren't excluded by sample/SNP in my GWAS data
    
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


#DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
if(!grepl("%%",dataFile)){
    DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
}

# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

# do we use linear regression results?
useLinReg=FALSE
if("-linReg"%in%args) useLinReg=TRUE

# which chroms?
if(chroms!="genome") {
    if(chroms=="all") {
        
        chroms = 1:22
        if( "-sex" %in% args ) chroms = 23:26
        
    } else {
        
        chroms = parse.range.string(chroms)
        
    }
}

print("printing the following chromosomes")
print( chroms )


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
    
    plot.BOLT.pvalues(GWASdata=gsub("%%",chrom,dataFile),chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss=maxmiss,plotOutDir=plotOutDir,extraTitle=extraTitle,highlight=highlightSNPs,highlightCols=highlightCols,Ymax=Ymax,useLinReg=useLinReg)

 
}



############# TESTING & USAGE #############

#Rscript plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-quant-Age-v1.out all plots > Logs/plot-Ychrom-BOLT-LMM-quant-Age-v1.log &

#Rscript plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LREG-quant-Age-v1.out all plots > Logs/plot-Ychrom-BOLT-LREG-quant-Age-v1.log &

#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-quant-Age-v1.out all plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-Age-v1.log &

    
#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-all-snps-quant-Age-v1.out all plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-all-snps-Age-all-v1.log &

#Rscript ../otherGWAS/plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.v1/Ychrom-BOLT-LMM-all-snps-quant-Age-v1.out genome plots -hi -title -QCcolors > Logs/plot-Ychrom-BOLT-LMM-quant-all-snps-Age-all-v1.log &


# dataFile="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v1/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-all-snps-v1.out"
