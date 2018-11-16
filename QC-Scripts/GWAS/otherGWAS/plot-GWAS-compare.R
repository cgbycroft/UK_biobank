h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

# any generic QC exclusions??
QCexclude = c()

# THIS COMPARES only sets where the output has exactly the same format (e.g only genotypes, or only imputed data).

##########
if(!"args"%in%ls()) args = commandArgs(TRUE)

                                        #args = c( "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Array.as.binary/BOLTLMM.v13/Array.as.binary-BOLT-LMM-v13-chr2.out", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Array.as.binary/BOLTLMM.v15/Array.as.binary-BOLT-LMM-v15-chr2.out", "2","plots","-lreg")


print(args)

dataFile1 = args[1]
dataFile2 = args[2]

chroms = args[3]
plotDir = args[4]

if("-noprint"%in%args) printOut=FALSE else printOut=TRUE

if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

if(chroms%in%c("all","genome")) chrom="genome" else chrom = parse.range.string(chroms)
if("-ymax"%in%args) ymax = as.numeric( args[which(args=="-ymax")+1]) else Ymax=50


if("-lreg" %in% args ) useLmmInf1 = useLmmInf2 = FALSE else useLmmInf1 = useLmmInf2 = TRUE
if("-lreg1" %in% args ) useLmmInf1 = FALSE else useLmmInf1=TRUE
if("-lreg2" %in% args ) useLmmInf2 = FALSE else useLmmInf2=TRUE

### read in GWAS results
if("-qc" %in% args){
                                        # Default QC THRESHOLDS  (used in plots made before 6-07-2016)
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


resultsGENO1 = read.gwas(dataFile1,chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf1)
DF1 = resultsGENO1$DF    

resultsGENO2 = read.gwas(dataFile2,chrom=chrom,minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf2)
DF2 = resultsGENO2$DF    

vers1 = gsub("\\.out","",str_split(dataFile1,"-v")[[1]][2])
vers2 = gsub("\\.out","",str_split(dataFile2,"-v")[[1]][2])
    
pvalset1=resultsGENO1$Pvalset

filename=gsub("-v.*.out",paste0("v",vers1,"-with-v",vers2),pvalset1)
if(!useLmmInf2) filename=gsub("-v.*.out",paste0("v",vers1,"-with-v",vers2,"-lreg"),pvalset1)
if(!useLmmInf1) filename=gsub("-v.*.out",paste0("v",vers1,"-lreg-with-v",vers2),pvalset1)



############################################
# Get the p-values where they overlap

x =  DF1$P2
y =  DF2$P2[match(DF1$SNP2,DF2$SNP2)]

print( paste0(dim(resultsGENO1$DFraw)[1]," raw snps in ",vers1) )
print( paste0(dim(resultsGENO2$DFraw)[1]," raw snps in ",vers2) )
print( paste0(dim(resultsGENO1$DF)[1]," filtered snps in ",vers1) )
print( paste0(dim(resultsGENO2$DF)[1]," filtered snps in ",vers2) )

print( paste0(sum(!is.na(y))," snps overlap."))

# don't plot anything below a particular value in both versions
excl = (x>10^-3)&(y>10^-3)

# Annotate with the snps that are significant in one but not even suggestive in the other!
xBetter = which((x<1e-3)&(y>=10^-2))
#yBetter = which((y<5e-8)&(x>=5e-8))
yBetter = which((y<1e-3)&(x>=10^-2))

xBet = DF1[xBetter,]
yBet = DF2[match(DF1$SNP2,DF2$SNP2),][yBetter,]

save(xBet,yBet,file=paste0(filename,"-comparison-outliers.RData"))



############################################
# Get snp classifications (rare etc.)

phasedColour = "red"
rareColour = "blue"
infoColour = "orange"
rareThreshold = 0.001

imputedData=FALSE
if( grepl("chr", dataFile1) ){ # are we comparing results from imputation?
    print("Data are imputed.")
    imputedData=TRUE
    chrom = as.numeric( gsub("\\.out","",str_split(basename(dataFile1),"chr")[[1]][2]) )
    phased = read.phased.snps(chrom)[[1]]
                                        # all based on DF1 values
    inPhased = DF1$alternate_ids%in%phased$V2
    poorInfo = DF1$INFO < 0.3
    
} else {
    
    p = read.phased.snps(chroms)
    phased = rbind_all(p)
    
    inPhased = DF1$SNP%in%phased$V2
    
    poorInfo = DF1$F_MISS >= 0.05
}

rare = DF1$MAF < rareThreshold

colorsGeno = rep("black",dim(DF1)[1]) 
colorsGeno[inPhased]=phasedColour # in phased



############################################
# Plot the p-values where they overlap


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

x =  DF1$P
y =  DF2$P[match(DF1$SNP2,DF2$SNP2)]


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
