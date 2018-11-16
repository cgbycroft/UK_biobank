# script to compare GWAS results from two versions

source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-pheno","Standing.height","-versions","v3,v5")
h = args[-c(which(args%in%c("-versions","-pheno")),which(args%in%c("-versions","-pheno"))+1)]
for(s in h) source(s)

print(args)

genotypeFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc"


#####################
versions = str_split(args[which(args=="-versions")+1],",")[[1]]
pheno = args[which(args=="-pheno")+1]
#####################

setwd(paste0( baseSampleQCDir,"/QC-Scripts/GWAS/otherGWAS/",pheno))

print("Reading QC snps lists...")
QCSNPList = read.SNPQC.files(justSNPs=TRUE)

QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$concordanceSNPs))    
# NOTE: anything created after Friday June 24th uses this QCexclude. Anything created before that also excluded the 'image' SNPs

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")


extraTitle="compare"


# read in results
results = list()

for(vers in versions){
    print(vers)
    dataFile = paste0("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/",pheno,"/BOLTLMM.",vers,"/Standing.height-BOLT-LMM-",vers,".out")
    r = read.gwas(dataFile,chrom="genome",minmaf=0,mininfo=0.3,maxmiss = 1,Ymax=50,QCexclusions=QCexclude,extraTitle=extraTitle)

    results[[vers]] = r
}

DF1 = results[[versions[1]]]$DF
DF2 = results[[versions[2]]]$DF

# how many SNPs in common?
snps1 = DF1$SNP
snps2 = DF2$SNP


################### ad-hoc
# is it the input data that's the problem?
v1inputBim = read.table(paste0( baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-oxfordqc.bim"),header=FALSE)
v2inputBim = read.table("/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas.bim",header=FALSE)

inBimv1 = v1inputBim[(v1inputBim$V1 %in% 1:22) & (!v1inputBim$V2 %in% QCexclude),]
inBimv2 = v2inputBim[(v2inputBim$V1 %in% 1:22),]
dim(inBimv1)
dim(inBimv2)

length(intersect(inBimv1$V2,inBimv2$V2))

###### in clare's gwas files and not in Colin's files.
notBimv2 = inBimv1[!inBimv1$V2%in%inBimv2$V2,]
g = sort(table(notBimv2$V4))
sum(g==2)
length(g)
# 418 SNPs which are in clare's gwas files and not in Colin's files.
# all but one have more than one snp associated with them!! I.e there are two snp ids for the same positions.

notBimv2[notBimv2$V4==(names(g)[g==1]),]
# namely: "Affx-92040534".  Actually this is not in the gwas results file. Must have been filtered out by bolt! > 10% missing ===> this is fine!

###### In Colin's gwas files and not in Clare's files.
notBimv1 = inBimv2[!inBimv2$V2%in%inBimv1$V2,]
# 62 SNPs not in clare's data but in Colin's.
# why did clare filter them out???? e.g Affx-35956907. THis snp fails hwe in all UKBiLEVE batches, and that is all. 
sum(notBimv1$V2%in%BB.ps2snp$AffySNPID)
sum(notBimv1$V2%in%BL.ps2snp$AffySNPID)

# They failed in all batches in which they would ever appear, so Colin's data is complete missing data, and Clare's data the SNPs aren't there at all. 

DF2$F_MISS[DF2$SNP%in%notBimv1$V2]

################### 


# compare common snps
toCompare = intersect(snps1,snps2)
DF1a = DF1[DF1$SNP%in%toCompare,]
DF2a = DF2[DF2$SNP%in%toCompare,]
sum(DF1a$SNP!=DF2a$SNP)


colors = rep(NA,dim(DF1a)[1])
colors[DF1a$MAF < 0.001] = "red"
colors[(DF1a$MAF >= 0.001) & (DF1a$MAF < 0.01)] = "purple"
colors[(DF1a$MAF >= 0.01) & (DF1a$MAF < 0.1)] = "blue"
colors[(DF1a$MAF >= 0.1) & (DF1a$MAF < 0.2)] = "green"
colors[(DF1a$MAF >= 0.2)] = "yellow"
Order = order.by.number.occurrences(colors)
diffs = -log10(DF1a$P) + log10(DF2a$P)
    

png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-log10pvalues.png"),width=1000,height=1000,res=150)
x = -log10(DF1a$P2)
y = -log10(DF2a$P2)
limits = c(0,max(c(x,y)))
plot(x,y,xlab=paste0("-log10(pval) ",versions[1]),ylab=paste0("-log10(pval) ",versions[2]),xlim=limits,ylim=limits)
abline(0,1,col="red")
dev.off()

png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-log10pvaluesDiff.png"),width=1000,height=1000,res=150)
plot(log10(DF1a$F_MISS[Order]),-log10(x[Order]) + log10(y[Order]),xlab=paste0("log10 Fraction missing ",versions[1]),ylab=paste0("-log10(pval): ",versions[1]," - ",versions[2]),col=colors[Order],pch=16)
abline(h=0,col="black",lty=3)
dev.off()

meanTest = t.test(-log10(DF1a$P2),-log10(DF2a$P2))

# any snps above the suggestive line moved in rank?
snps = (DF1a$P < 1e-05) | (DF2a$P < 1e-05)
ranks1 = rank(DF1a$P[snps])
ranks2 = rank(DF2a$P[snps])

sum( (ranks1-ranks2) > 500)


png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-log10pvaluesDiff-suggestiveSNPs.png"),width=1000,height=1000,res=150)
plot(log10(DF1a$F_MISS[snps]),-log10(x[snps]) + log10(y[snps]),xlab=paste0("log10 Fraction missing ",versions[1]),ylab=paste0("-log10(pval): ",versions[1]," - ",versions[2]),col=colors[snps],pch=16)
abline(h=0,col="black",lty=3)
dev.off()


png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-log10pvaluesDiff-MissDiff.png"),width=1000,height=1000,res=150)
plot(-log10(DF1a$F_MISS[snps]) + log10(DF2a$F_MISS[snps]),-log10(x[snps]) + log10(y[snps]),xlab=paste0("-log10(missing): ",versions[1]," - ",versions[2]),ylab=paste0("-log10(pval): ",versions[1]," - ",versions[2]),col=colors[snps],pch=16)
abline(h=0,col="black",lty=3)
dev.off()


png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-maf.png"),width=1000,height=1000,res=150)
limits = c(0,0.5)
plot(DF1a$MAF,DF2a$MAF,xlab=paste0("MAF ",versions[1]),ylab=paste0("MAF ",versions[2]),xlim=limits,ylim=limits)
abline(0,1,col="red")
dev.off()

png(paste0("plots/",pheno,"-compareGWAS-",versions[1],"-and-",versions[2],"-log10pvaluesDiff-MAFDiff.png"),width=1000,height=1000,res=150)
plot(log10(DF1a$MAF[snps]) - log10(DF2a$MAF[snps]),-log10(x[snps]) + log10(y[snps]),xlab=paste0("log10(MAF): ",versions[1]," - ",versions[2]),ylab=paste0("-log10(pval): ",versions[1]," - ",versions[2]),col=colors[snps],pch=16)
abline(h=0,col="black",lty=3)
dev.off()



toCheck = DF1a[(abs(-log10(x) + log10(y)) > 0.01) & snps,]
