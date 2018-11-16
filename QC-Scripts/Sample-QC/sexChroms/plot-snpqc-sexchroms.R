## Script to plot output from snp-qc, as it relates the sex chromosomes. Also produces some cluster plots.
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

args = commandArgs(trailingOnly=T)

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-out","b1__b11-b001__b095-sexchrom-sampleqc-snpqc")

print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}

outFile = args[which(args=="-out")+1]

# read in Snp-qc files
snpQC = read.SNPQC.files(type="sexchrom",justSNPs=FALSE)

# get snp annotations (position etc)
UKBB = ukbiobank.ps2snp("sexchrom")
UKBL = ukbileve.ps2snp("sexchrom")
UKBL = tbl_df(UKBL)
UKBB = tbl_df(UKBB)

ps2snp = merge(UKBL, UKBB, all = TRUE)
ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]

# set par boundaries (from ucsc browser, hg19)
par1 = 2699520
par2 = 154931044

# plot p-values of m-female test (minimum within each batch) with colours indicating other tests
mfraw = tbl_df(snpQC$mfSNPs)
mfU = unique(mfraw$V2)
lines = sapply(mfU,function(x) {
    m = min(mfraw$V3[mfraw$V2==x])
    line = which((mfraw$V2==x)&(mfraw$V3==m))[1]
})

mf = left_join(mfraw[lines,],ps2snp,by=c("V2"="AffySNPID"))
mf$chrom = 23

png( paste0("plots/",outFile,"-male-female-pvals-X.png") ,width=2000,height=1000,res=150)
myManhattan(mf,p="V3",bp="Position",chr="chrom",snp="V2",ymax=50)
dev.off()

png( paste0("plots/",outFile,"-male-female-pvals-PAR1.png") ,width=2000,height=1000,res=150)
myManhattan(mf,p="V3",bp="Position",chr="chrom",snp="V2",ymax=50,
            xlimits=c(min(ps2snp$Position[ps2snp$IsInPAR1==1]),par1+10^4) )
abline(v=par1,lty=3)
dev.off()

png( paste0("plots/",outFile,"-male-female-pvals-PAR1-boundary.png") ,width=2000,height=1000,res=150)
myManhattan(mf,p="V3",bp="Position",chr="chrom",snp="V2",ymax=50,
            xlimits=c(par1-10^5,par1+10^5) )
abline(v=par1,lty=3)
dev.off()

png( paste0("plots/",outFile,"-male-female-pvals-PAR2.png") ,width=2000,height=1000,res=150)
myManhattan(mf,p="V3",bp="Position",chr="chrom",snp="V2",ymax=50,
            xlimits=c(par2-10^5,max(ps2snp$Position)) )
abline(v=par2,lty=3)
dev.off()

# how many actually fail other tests, in the given batch?

otherFails = apply(mfraw,1,function(x){
    b = x[[1]]
    snp = x[[2]]
    inThis = sapply(names(snpQC)[c(1,2,4,6,7)],function(j){
        print(j)
        i=snpQC[[j]]
        if(j%in%c("batchHweSNPs","plateSNPs","batchSNPs")) t = i[i$V1==b,2] else t = i[,1]
        snp%in%t
    })
})

keep = colSums(otherFails)==0 # snps that don't fail any other test
mf2raw = mfraw[keep,]
mfU = unique(mf2raw$V2)
lines = sapply(mfU,function(x) {
    m = min(mf2raw$V3[mf2raw$V2==x])
    line = which((mf2raw$V2==x)&(mf2raw$V3==m))[1]
})

mf = left_join(mf2raw[lines,],ps2snp,by=c("V2"="AffySNPID"))


#### only plot SNPs that don't fail other tests (77 snps)

# cluster plots of snps in PAR1
# select batches to plot (represent all snps)
getMaxBatches <- function(snps){
    batches = all.batches()
    bs = sapply(snps,function(x) batches%in%mfraw$V1[mfraw$V2==x])
    
    # select a set that covers all snps
    batOrder = batches[order(rowSums(bs),decreasing=TRUE)]
    bat = batOrder[1:2]
    for(b in batOrder[3:length(batOrder)]){
        if( sum(colSums(bs[batches%in%bat,])==0)>1 ) bat = c(bat,b) else break
    }
    return(bat)
}

snps = mf$V2[mf$Position <= par1]
bats = getMaxBatches(snps)[1:4]
write.table(snps,file="snps-temp1.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data snps-temp1.txt ",baseSampleQCDir,"/QC-Scripts/Sample-QC/sexChroms/plots/clusterPlots_PAR1_mf_fails -b ",paste(bats,collapse=",")," -colour Inferred.Gender -orig" ) )

# cluster plots of snps in PAR2
snps = mf$V2[mf$Position >= par2]
bats = getMaxBatches(snps)[1:4]
write.table(snps,file="snps-temp2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data snps-temp2.txt ",baseSampleQCDir,"/QC-Scripts/Sample-QC/sexChroms/plots/clusterPlots_PAR2_mf_fails -b ",paste(bats,collapse=",")," -colour Inferred.Gender -orig") )

# cluster plots of snps nearPAR1
snps = mf$V2[mf$Position <= (par1 + 10^5)]
bats = getMaxBatches(snps)[1:4]
write.table(snps,file="snps-temp3.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data snps-temp3.txt ",baseSampleQCDir,"/QC-Scripts/Sample-QC/sexChroms/plots/clusterPlots_nearPAR1_mf_fails -b ",paste(bats,collapse=",")," -colour Inferred.Gender -orig" ) )

