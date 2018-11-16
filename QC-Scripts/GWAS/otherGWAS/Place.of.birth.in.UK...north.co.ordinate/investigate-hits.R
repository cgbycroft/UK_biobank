source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
plink="/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink"
for(s in h) source(s)

library(dplyr)
library(qqman)
library(stringr)
library(binom)
library(sp)
library(rgdal)

sexChroms = c(23,24,25,26)
names(sexChroms) = c("X","Y","XY","MT")


dataFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v1/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-all-snps-v1.out"
# read phenodata
phenoData = paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v1.txt")
pheno = read.table(phenoData,header=TRUE,stringsAsFactors=FALSE)


plotOutDir = "plots"
extraTitle=""



# Should we highlight some SNPs?
highlightSNPs=NULL
highlightCols=NULL

if("-hi"%in%args){
    print("Reading QC snps lists...")
    snpList = read.SNPQC.files()

    highlightSNPs = unique(unlist(snpList))
    highlightCols = rep("black",length(highlightSNPs))
    highlightCols[highlightSNPs%in%snpList$hweSNPs] = "green"    # HWE (apply first)
    highlightCols[highlightSNPs%in%snpList$imageSNPs] = "orange"  # IMAGE ARTEFACT
    highlightCols[highlightSNPs%in%snpList$arraySNPs] = "red" # ARRAY  
    highlightCols[highlightSNPs%in%snpList$concordanceSNPs] = "blue"  # CONCORDANCE
    
    print(table(highlightCols))
}


# get data

print("printing the following chromosomes")
print( chroms )


#DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
if(!grepl("%%",dataFile)){
    DFraw = read.table(dataFile,sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
}

# which chroms?
if(chroms!="genome") {
    if(chroms=="all") chroms = 1:22 else chroms = parse.range.string(chroms)
}


if(grepl("%%",dataFile)) {
    DFraw = read.table(gsub("%%",chr,dataFile),sep="",header=TRUE,stringsAsFactors=FALSE)
    DFraw$MAF = DFraw$A1FREQ
    DFraw$MAF[DFraw$A1FREQ > 0.5] = 1-DFraw$A1FREQ[DFraw$A1FREQ > 0.5]
}

#print(chrom)
chrom="genome"
minmaf = 0.001
mininfo = 0.3
maxmiss = 0.05  # maximum 5% missing data    
GWASdata = dataFile

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

Ymax = ceiling(max(-log10(DF$P[DF$P!=0]),na.rm=T)) + 10
Ymax = min(Ymax,50)


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
for(i in unique(DF$CHR)){        
    if(i>1) DF$index[DF$CHR==i] = DF$index[DF$CHR==i] + max(DF$BP[DF$CHR==(i - 1)])
}

if("BETA"%in%colnames(DF)){
    
    snps = which((DF$P < 5e-8)&(!is.na(DF$P)))    
    beta = DF
    beta$BETA2 = beta$BETA
    beta$BETA2[DF$A1FREQ > 0.5] = -beta$BETA2[DF$A1FREQ > 0.5]
    beta = beta[snps,]

    png(paste(plotOutDir,"/",Pvalset,"-EffectSizes.png",sep=""),height=1000,width=1000,res=150)
    myManhattan(beta,p="BETA2",logtransform=FALSE,genomewideline=0,suggestiveline=FALSE)
    
    dev.off()
    
}

########## plot haplotypes


reg = "chrom2"
region = (beta$CHR==2) & (beta$BP >=135755629)
pos = beta$BP[region]
alleles = beta$SNP[region]
write.table(alleles,file="tmp",quote=FALSE,col.names=FALSE,row.names=FALSE)

system( paste0(plink," --bfile ",baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-oxfordqc --extract tmp --keep-allele-order --recode AD --out b1__b11-b001__b095-autosome-oxfordqc-hit-region-",reg ))

raw = read.table( paste0("b1__b11-b001__b095-autosome-oxfordqc-hit-region-",reg,".raw"),header=TRUE,stringsAsFactors=FALSE)

# subset for samples with non-missing phenotype values
samples = pheno$IID[!is.na(pheno$Place.of.birth.in.UK...north.co.ordinate)]
raw = raw[raw$IID%in%samples,]
snpNames = t(read.table( paste0("b1__b11-b001__b095-autosome-oxfordqc-hit-region-",reg,".raw"),header=FALSE,stringsAsFactors=FALSE,nrow=1)[,-c(1:6)])[,1]
snpNames = sapply(snpNames[!grepl("HET",snpNames)],function(x) str_split(x,"_")[[1]][1])

betas = beta$BETA[match(snpNames,beta$SNP)]

getScore = function(hap){
#    sum(hap*betas,na.rm=TRUE)
        sum(hap*betas)
}

HAPS = raw[,-c(1:6,grep("HET",colnames(raw)))]
scores = apply(HAPS,1,getScore)

strings = apply(HAPS,1,paste,collapse="")

s =sort(strings[!grepl("NA",strings)])

# NOTE: Betas are fold-departure from mean of phenotype, for every copy of the derived allele (Allele1)

these = match(raw$IID,pheno$IID)
colors = colour.scale(scores)

x = pheno$Place.of.birth.in.UK...east.co.ordinate[these]
y = pheno$Place.of.birth.in.UK...north.co.ordinate[these]
y[y==-1]=NA
x[x==-1]=NA

#NOTE: grid coordinates are in OSGB1936
# see http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=118122
#"Measurements refer to easting and northing with a reference point near the
# Isles of Sicily"

png(paste(plotOutDir,"/",Pvalset,"-haplotypeScoreMap-",reg,".png",sep=""),width=1500,height=2000,res=150)
plot(x,y,col=colors,pch=16)
dev.off()    

# genotype counts
HOM0 = apply(HAPS,2,function(x) sum(x ==0,na.rm=TRUE))
HOM1 = apply(HAPS,2,function(x) sum(x ==2,na.rm=TRUE))
HET = apply(HAPS,2,function(x) sum(x ==1,na.rm=TRUE))

freqs = (2*HOM1 + HET)/(2*(HOM1 + HET + HOM0))

hw = hwe.test(HOM0,HET,HOM1)

# print some cluster plots

system( paste0("Rscript /well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(snpNames,collapse=",")," /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/clusterPlots_",reg," -b Batch_b034-Batch_b038") )

#system( "wget https://cran.r-project.org/src/contrib/rgdal_1.1-10.tar.gz -P /well/ukbiobank/qcoutput.V2_QCed.sample-QC/src" )
packageurl = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/rgdal_1.1-10.tar.gz"
install.packages(packageurl,repos=NULL,type="source")
install.packages('rgdal', type = "source", configure.args=c('--with-proj-include=/usr/local/include','--with-proj-lib=/usr/local/lib'))

# plot frequency in bins of north coordinate (takes a few minutes)
binwidth=10000
bins = seq(min(y,na.rm=TRUE),max(y,na.rm=TRUE),binwidth)
na = rep(NA,dim(HAPS)[2] )
binFreq = sapply(bins,function(bin,snp){
    print(bin)
    s = (y < (bin + binwidth))&(y >= bin)&(!is.na(y))
    if(sum(s)>1){ # only when more than one individual in bin
        h = HAPS[s,]        
        totals = apply(h,2,function(k) 2*sum(!is.na(k)) )
        f = binom.confint(colSums(h,na.rm=TRUE),n=totals, tol = 1e-8,method="exact")
    } else {
#        f=na
        f = NA
    }
    return(f)
})

png(paste(plotOutDir,"/",Pvalset,"-freqsByNorthCoordinate-",reg,"-%02d.png",sep=""),width=1500,height=2000,res=150)
for(i in 1:nrow(binFreq[[1]])){
    Y = sapply(binFreq,function(x){
        if(!is.na(x)) x[i,"mean"] else NA
    })
    int1 = sapply(binFreq,function(x){
        if(!is.na(x)) x[i,"upper"] else NA
    })
    int2 = sapply(binFreq,function(x){
        if(!is.na(x)) x[i,"lower"] else NA
    })

    r = round(range(Y,na.rm=TRUE),2)
    plot(bins/1000,Y,main=paste0(snpNames[i],"\nFrequency range: ",r[1]," - ",r[2]) ,ylab="Frequency of A1 allele",xlab = paste0("Distance from Isles of Sicily (",binwidth/1000,"km bins)"))
    segments(bins/1000,int2,bins/1000,int1)
}
dev.off()


## Create spatial grid
bb <- rbind(c(0,max(x,na.rm=TRUE)),c(0,max(y,na.rm=TRUE)))
cs <- c(100,100) # cell size (in units of proj4string(MySmallerMapSimpledf))
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd

sp_grd <- SpatialGridDataFrame(grd,data=data.frame(id=1:prod(cd)))
Points <- cbind(x,y)
    
#sp_grd2 =  spTransform(sp_grd,CRS("+proj=merc +units=m"))
#Points2 = spTransform(Points,CRS("+proj=merc +units=m"))

# Get distances between each cell and each point on the map (in metres)
distances <- gDistance(Points2,sp_grd2,byid=T) # this is a long step if you have many cells


# Get lm results for bivariate (2D) outcome variable at each SNP

Y = cbind(x,y)
d1 = lm(Y~HAPS[,1])
d2 = lm(Y~HAPS[,2])

# copy this data to my computer and analyse from there...
pheno2 = pheno[,c("IID","Place.of.birth.in.UK...north.co.ordinate","Place.of.birth.in.UK...east.co.ordinate")]
save(raw,HAPS,pheno2,file="dataForLaptop.RData")
