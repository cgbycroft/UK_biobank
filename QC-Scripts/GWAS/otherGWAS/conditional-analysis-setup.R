#############
# This script takes a list of SNPs that are 'top' snps in a GWAS, and produces a phenotype file that contains genotypes for these SNP as a covariate.
#############
args = commandArgs(TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

#args = c("/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr%%","-vers","v16","-snps","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/Standing.height-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData")

for(s in h){
    source(s)
}

print(args)

vers = args[which(args=="-vers")+1]
#phenoFile = paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt")
phenoFile = paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt")  # NOTE THIS IS HARD-CODED! 
genoFile = args[1]
if(grepl("%%",genoFile)) genoFileRaw = genoFile

snpFile = args[which(args=="-snps")+1]
setwd('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS')

BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")


############
# read in pheno data to get order and list of samples
pheno = read.table(phenoFile,header=TRUE,stringsAsFactors=FALSE)


# define regions to do conditional GWAS on.
load(file=snpFile,verbose=TRUE)


# run plink to get genotypes in 0,1,2 format, and by chromosome (as need to read in this file and reorder rows)
allTopVariantsFile = gsub("-hitsAndLDcalc.RData","-topVariants.txt",snpFile)
for(i in 1:22){
    print(i)

    if("genoFileRaw"%in%ls()) genoFile = gsub("%%",i,genoFileRaw)

    system(paste0(plink,' --bfile ',genoFile,' --chr ',i,' --extract ',allTopVariantsFile,' --keep-allele-order --recode A --out ',dirname(snpFile),'/genotypes-topVariants-',vers,'-chr',i))   
}

# all variants
system(paste0(plink,' --bfile ',genoFile,' --extract ',allTopVariantsFile,' --keep-allele-order --recode A --out ',dirname(snpFile),'/genotypes-topVariants-',vers))

# sort by phenotype table so we just have to cat the appropriate columns of this file in each GWAS run
for(i in 1:22){
    
    print(i)
    g = read.table(paste0(dirname(snpFile),'/genotypes-topVariants-',vers,'-chr',i,'.raw'),header=TRUE,check.names=FALSE)
    g2 = g[match(pheno$IID,g$IID),-c(1,3:6)]
    #n = str_split(colnames(g2)[-1],"_")
    #n = sapply(n,function(x) x[[1]])
    #snpNames = gsub("\\.","-",n)
    snpNames = colnames(g2)[-1] # this only works when check.names=FALSE
    print( sum(!snpNames%in%BB.ps2snp$AffySNPID) )
    colnames(g2)[-1] = snpNames
    write.table(g2,file=paste0(dirname(snpFile),'/genotypes-topVariants-',vers,'-phenoOrder-chr',i,'.txt'),quote=FALSE,row.names=FALSE,col.names=TRUE)
    print( sum(g2$IID!=pheno$IID,na.rm=TRUE) )
}



# get set of contiguous hit regions
regionsChrom = list()
uniqueRegionsChrom = list()
for(i in 1:22){
    print(i)

    HITS = HITSchrom[[as.character(i)]]
    regs = HITS$regs#[order(HITS$regs[,1]),]
    regsOrder = order(regs[,1])
                                        #    d = regs[regsOrder,][-1,1] - regs[regsOrder,][-nrow(regs),2]
    start = min(regs[,1])
    end = min(regs[,2])
    regsStart = regsEnd = index = rep(NA,nrow(regs))
    r=1
    while(sum(is.na(regsStart)) > 0){
        these = (regs[,1]<=end) & ( is.na(regsStart) )
        if(sum(these) > 0 ){
            these = (regs[,1]<=end) & (regs[,1]>=start)
            end = max(regs[these,2])
            regsStart[these] = start
            regsEnd[these] = end
            index[these] = r
        }
        if(sum(these)==0){
            r = r + 1
            start = min(regs[regs[,1]>end,1])
            end = min(regs[regs[,1]>end,2])       
        }        
    }
    regions = cbind(regs,regsStart,regsEnd,index)
    print( paste0( length(table(index)), " non-overlapping regions"))
    regionsChrom[[as.character(i)]] = regions
    region=list()
    for( r in sort(unique(index))) {
        snps = HITS$topVariants[regions[,5]==r]
        region[[r]] = list(unique(regions[regions[,5]==r,3:4]),snps)
    }
    uniqueRegionsChrom[[as.character(i)]] = region
}

          
# print list of start and end positions
outputFile = gsub("hitsAndLDcalc.RData","snpRegionList.txt", snpFile )
write.table(t(c("chrom region start end snps")),file=outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE)
system(paste0('mkdir ',dirname(outputFile),"/hitRegions"))
sapply(1:length(uniqueRegionsChrom),function(i) {
    R = uniqueRegionsChrom[[as.character(i)]]
    sapply(1:length(R),function(r){
        ranges = uniqueRegionsChrom[[as.character(i)]][[r]][[1]]
        snps = uniqueRegionsChrom[[as.character(i)]][[r]][[2]]
        snpList = paste0(dirname(snpFile),"/hitRegions/chr",i,"-region-",r,"-snpList.txt")
        write.table(snps,file=snpList,quote=FALSE,col.names=FALSE,row.names=FALSE)
        write.table(t(c(i,r,ranges[[1]],ranges[[2]],snpList)),file = outputFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)
    })               
})
    

# check that regions are really non-overlapping
for(chrom in 1:22){
    starts = unlist(sapply(uniqueRegionsChrom[[chrom]],function(y) y[[1]][1]) )
    ends = unlist(sapply(uniqueRegionsChrom[[chrom]],function(y) y[[1]][2]) )
    print( sum((starts[-1] - ends[-length(ends)])<0) )
}


print("DONE!")
