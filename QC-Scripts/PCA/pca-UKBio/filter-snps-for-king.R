
args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/auxFunctions.R",
#"-in","b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init1","-threshold","0.003")

print(args)
h = args[-c(which(args%in%c("-in","-threshold")),1+which(args%in%c("-in","-threshold")))]
for(helperScript in h){
    source(helperScript)
}

OutputFile = args[which(args=="-in")+1]
threshold = as.numeric(args[which(args=="-threshold")+1])


SnpLoads = paste0(baseSampleQCDir,'/data/PCA/',OutputFile,'.snpload.map')
loads = dplyr::tbl_df(read.table(SnpLoads,header=F,stringsAsFactors=F))

# look at the first 5 PCs
# standard-deviations are very similar, 0.00298768
# what about the distributions?

system(paste0('mkdir plots'))


for(pc in 1:5){
    png(paste0('plots/',OutputFile,'-PC',pc,'-loads-hist-%02d.png'),width=1500,height=1500,res=150)
    hist(loads[[pc + 7]],breaks=1000,xlab=paste0('PC ',pc),main=paste0("SNP-loads PC ", pc,"\nVertical lines shown at ",-threshold," and ",threshold))
#    abline(v=mean(loads[[pc + 7]]),col="red")
    abline(v=c(-threshold,threshold),col="red")
    hist(abs(loads[[pc + 7]]),breaks=1000,xlab=paste0('PC ',pc),main=paste0("Absolute value SNP-loads PC ", pc))
#    abline(v=mean(loads[[pc + 7]]),col="red")
    abline(v=c(threshold),col="red")
    dev.off()
}


# select a set of SNPs
for(threshold in c(seq(0.003,0.004,by=0.0005),seq(0.005,0.01,by=0.001)) ){
    print(threshold)

    pc1 = abs(loads[[8]]) > threshold
    pc2 = abs(loads[[9]]) > threshold
    pc3 = abs(loads[[10]]) > threshold
    pc4 = abs(loads[[11]]) > threshold

    hi = pc1 + pc2 + pc3
    exc = sum(hi>0) # extreme in at least one
    inc = sum(hi==0)

    print(paste0(exc," SNPs will be in the exclusion list. That leaves ",inc," remaining, out of ",nrow(loads)))


                                        # 1-4: 34,364 SNPs left.
                                        # 1-3: 50,959 SNPs left. ~ 40 %  
                                        # 1-2: 63,835 SNPs left.
                                        # 1: 81,370 SNPs left. 
    ## NOTE: all are left with around 80,000 SNPs using 0.003 threshold
                                        # exclude remaining SNPs, plus SNPs in LD

    excludeSNPs = loads$V2[hi>0]
    keepSNPs = loads$V2[hi==0]

                                        # save list of SNPs
    write.table(excludeSNPs,file=paste0(OutputFile,"-snpsToExcludePCs-",threshold,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(keepSNPs,file=paste0(OutputFile,"-snpsToKeepPCs-",threshold,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

    
}

write.table(seq(0.003,0.01,by=0.001),paste0(OutputFile,"-snpsToKeepPCs-thresholds.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
