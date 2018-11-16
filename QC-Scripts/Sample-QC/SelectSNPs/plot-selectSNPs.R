args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")


for(helperScript in args){
    source(helperScript)
}

type = "sexchrom"
type="autosome"
batches = all.batches()

sexchrom = attach(paste0('snpIncludeList-sexchrom.Rdata'))
autosome = attach(paste0('snpIncludeList-autosome.Rdata'))


for(type in c("sexchrom","autosome")){
    UKBL = paste0(baseDataDir,'/UKBiLEVEAX_b1/sorted_bed/UKBiLEVEAX_b1-',type,'.bim') # NOTE: this assumes that all batches in UKBiLEVE have the same set
    UKBB = paste0(baseDataDir,'/Batch_b001/sorted_bed/Batch_b001-',type,'.bim') # NOTE: this assumes that all batches in UKBiobank have the same set of SNPs in their plink files!!! Which is true.
    ps2snpBL = ukbileve.ps2snp(type)  # all snps in V2_QCed, see snpqc-tests.R
    ps2snpBB = ukbiobank.ps2snp(type)  # all snps in V2_QCed

    print(UKBL)
    print(UKBB)

    UKBLsnps = read.table(UKBL,header=F,stringsAsFactors=F)
    UKBBsnps = read.table(UKBB,header=F,stringsAsFactors=F)
    
    allsnps1 = unique(c(ps2snpBL$AffySNPID,ps2snpBB$AffySNPID))
    allsnps = unique(c(UKBLsnps$V2,UKBBsnps$V2))
    
    if(type=="sexchrom") {
        sexchrom$allsnps = allsnps
        sexchrom$snpsOnly = unique(c(ps2snpBL$AffySNPID[ps2snpBL$Variant=="SNP"],ps2snpBB$AffySNPID[ps2snpBB$Variant=="SNP"]))
        sexchrom$overlap2 = intersect(UKBLsnps$V2,UKBBsnps$V2) # including INdels
    }
    if(type=="autosome") {
        autosome$allsnps = allsnps
        autosome$snpsOnly = unique(c(ps2snpBL$AffySNPID[ps2snpBL$Variant=="SNP"],ps2snpBB$AffySNPID[ps2snpBB$Variant=="SNP"]))
        autosome$overlap2 = intersect(UKBLsnps$V2,UKBBsnps$V2)

    }
}

#### TOTAL SNPS IN OUTPUT GENOTYPES FILE (31/10/2016) = 805518
# i.e /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.bim

# merge lists for genome total
total = list()
keys <- batches
total[["affyFail"]] <- c(sexchrom$affyFail, autosome$affyFail)
total[["oxFail"]] <- c(sexchrom$oxFail, autosome$oxFail)
total[["lowFreq"]] <- c(sexchrom$lowFreq,autosome$lowFreq)
total[["overlap"]] <- c(sexchrom$overlap, autosome$overlap)
total[["overlap2"]] <- c(sexchrom$overlap2, autosome$overlap2)
total[["imageSNPs"]] <- c(sexchrom$imageSNPs, autosome$imageSNPs)
total[["snpsToInclude"]] <- c(sexchrom$snpsToInclude, autosome$snpsToInclude)
total[["snpsOnly"]] <-  c(sexchrom$snpsOnly, autosome$snpsOnly) # snps only
total[["allsnps"]] <-  c(sexchrom$allsnps, autosome$allsnps) # includes INDELS!

# get xy snps only for sexchroms
xyChroms = list()
ps2snpBLsex = ukbileve.ps2snp("sexchrom")  # all snps in V2_QCed, see snpqc-tests.R
ps2snpBBsex = ukbiobank.ps2snp("sexchrom")  # all snps in V2_QCed

xySnps =unique(c(ps2snpBLsex$AffySNPID[ps2snpBLsex$Chromosome%in%c(23:24)],ps2snpBBsex$AffySNPID[ps2snpBBsex$Chromosome%in%c(23:24)]))


keys <- batches
xyChroms[["affyFail"]] <- intersect(sexchrom$affyFail,xySnps)
xyChroms[["oxFail"]] <- intersect(sexchrom$oxFail, xySnps)
xyChroms[["lowFreq"]] <- intersect(sexchrom$lowFreq, xySnps)
xyChroms[["overlap"]] <- intersect(sexchrom$overlap, xySnps)
xyChroms[["overlap2"]] <- intersect(sexchrom$overlap2, xySnps)
xyChroms[["imageSNPs"]] <- intersect(sexchrom$imageSNPs, xySnps)
xyChroms[["snpsToInclude"]] <- intersect(sexchrom$snpsToInclude, xySnps)
xyChroms[["snpsOnly"]] <-  intersect(sexchrom$snpsOnly, xySnps) # snps only
xyChroms[["allsnps"]] <-  intersect(sexchrom$allsnps, xySnps) # includes INDELS!

    
#1. PRINT OUT SOME NUMBERS:

sink(file="SelectSNPs-numbers-summmary.txt")

for(type in c("sexchrom","autosome","total","xyChroms","sexchrom2","autosome2","xyChroms2") ){

    l = get(type)

    
    uniqueFailures = unique(c(unlist(l$affyFail),unlist(l$oxFail)))
    uniqueFailuresAffy = unique(unlist(l$affyFail))
    uniqueFailuresOxford = unique(unlist(l$oxFail))
    allsnps = unique(l$allsnps)

    excludeSNPs = unique( c(allsnps[!allsnps%in%l$snpsOnly],allsnps[!allsnps%in%l$overlap2],uniqueFailures,l$lowFreq,l$imageSNPs) )
    excludeSNPs = intersect(excludeSNPs,allsnps)

###### CHECK AGAINST THE ACTUAL FILES USED IN SAMPLE QC

    if( !type %in%c( "total","xyChroms","autosome2","sexchrom2","xyChroms2" ) ) {
        print("Checking against actual data files used in sample QC. All below should be TRUE")
        
        combinedPlink = paste0( baseSampleQCDir,'/data/Combined/b1__b11-b001__b095-',type,'-sampleqc.bim')

        combined = read.table(combinedPlink,header=FALSE,stringsAsFactors=FALSE)
        print( dim(combined)[1] == length(allsnps) - length(excludeSNPs) )
        print( dim(combined)[1] == length(l$snpsToInclude) )

        print( sum( combined$V2%in%l$snpsToInclude ) == length(combined$V2) )
        print( sum( l$snpsToInclude%in%combined$V2 ) == length(l$snpsToInclude) )

    }

    
    
    print(paste0(length(allsnps)," snps to start off with."))
    
    print(paste0(type,":"))
    
    print(paste0("# Markers on both arrays (includes Indels):  ",length(unique(intersect(allsnps,l$overlap2) ) )))
    print(paste0("# Markers on both arrays (excludes Indels):  ",length(unique(intersect(allsnps,l$overlap) ) )))

    print(paste0("# Markers NOT on both arrays (includes Indels):  ",length(unique(allsnps[!allsnps%in%l$overlap2]) )))

    print(paste0("# SNPs as opposed to indels:  ",length(unique(intersect(allsnps,l$snpsOnly) ) )))
    print(paste0("# Markers that are NOT snps:  ",length(unique(allsnps[!allsnps%in%l$snpsOnly]))))

    print(paste0("# Markers in image artefact:  ",length(unique(intersect(allsnps,l$imageSNPs) ) )))

    print(paste0("# Markers that failed in at least one batch (on either array): ",length(unique(intersect(allsnps,uniqueFailures)) ) ))
    print(paste0("# Markers that failed in at least one batch (on either array) - AFFY: ",length(unique(intersect(allsnps,uniqueFailuresAffy)) ) ))
    print(paste0("# Markers that failed in at least one batch (on either array) - OXFORD: ",length(unique(intersect(allsnps,uniqueFailuresOxford)) ) ))
    print(paste0("# Markers with MAF less than ",lowMAF," : ",length(unique(intersect(l$lowFreq,allsnps)) ) ) )

    print(paste0("# Markers included in sample QC: ",length(l$snpsToInclude) ) )
    print(paste0("# Markers excluded in sample QC: ",length(excludeSNPs) ) )

    print( length(allsnps) - length(excludeSNPs) )
    
    success = length(l$snpsToInclude)/length(l$overlap2)
    print(paste0("Success rate  = (# Markers included in sample QC) / (# Markers on both arrays) = ",100*success,"%") )



print("\n")

}

sink()

quit()




#################
#AAGGGGH. What was the actual allele frequency threshold that I effectively used!?

freqs = read.genotyped.maf( snpFrequencyData= "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr" )

# all snps excluded based on this criteria
lowf = total[["lowFreq"]]
freqs$AffyID = ps2snpUnion$AffySNPID[match(freqs$SNP,ps2snpUnion$dbSNPRSID)]

freqsLOWF = freqs[freqs$AffyID%in%lowf,]




    freq = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + Counts$AB_all_female)
    # SHOULD HAVE BEEN - with the 2x at the bottom:
#    freq2 = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + 2*Counts$AB_all_female)

freqmaf = freq; freqmaf[freq>0.5] = 1-freq[freq>0.5]
freq2maf = freq2; freq2maf[freq2>0.5] = 1-freq2[freq2>0.5]



print( paste0(length(unique(lowFreq))," SNPs with MAF lower than or equal to ",lowMAF ) )

png("plots/test.png")
plot(freqmaf,freq2maf,xlim=c(0,0.5),ylim=c(0,0.5))
abline(0,1,col="red")
dev.off()

excl = rownames(Counts)%in%lowf

    
png("plots/test2.png")
plot(freqmaf,freq2maf,xlim=c(0,0.01),ylim=c(0,0.1))
points(freqmaf[excl],freq2maf[excl],col="purple")
points(freqmaf[exclany],freq2maf[exclany],col="red")

abline(v=0.0001,col="blue")
dev.off()

####################################################################
# What I should have done:

for(type in c("autosome","sexchrom")){
    
    print(type)
    # use genotype counts file created in SNP-QC
    rm("Counts")

    for (batch in batches){
        print(batch)
            
        countsFile = paste0(baseDataDir,"/",batch,"/genotype_counts/",batch,"-",type,"-genotype_counts.txt")
        counts = read.table(countsFile,header=T,colClasses=c("character",rep("numeric",8),rep("NULL",16)))
        totals = counts[,c("AA_all_female","AB_all_female","BB_all_female","NC_all_female")] + counts[,c("AA_all_male","AB_all_male","BB_all_male","NC_all_male")] 

        totals[is.na(totals)] = 0  # set NAs to zero
        
        if("Counts"%in%ls()) Counts = totals + Counts else Counts = totals
    }
    
    rownames(Counts) = counts[,1]
    
#### THIS IS ACTUALLY WRONG (but I used it...)!!!    
   # freq = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + Counts$AB_all_female)
    # SHOULD HAVE BEEN - with the 2x at the bottom:
    freq = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + 2*Counts$AB_all_female)

    names(freq) = rownames(Counts)
    assign(paste0(type,".freq"),freq)
}


for(type in c("autosome","sexchrom")){

#5. VERY RARE SNPs (near monomorphic - 1/10000)
    ##############
    lowMAF = 0.0001 # <==== this is the threshold I used with the bunged-up version
    ##############
#    lowMAF = 0.001
    
    freq = get(paste0(type,".freq"))
#    if(paste0(type,".freq")%in%ls()) freq = get(paste0(type,".freq"))
    # This is not acutally only females. The names of the columns are just inherited
    lowFreq = names(freq)[(freq <= lowMAF)|((1-freq) <= lowMAF)]
    
    print( paste0(length(unique(lowFreq))," SNPs with MAF lower than or equal to ",lowMAF ) )

    assign(paste0(type,".lowFreq"),lowFreq)

    l=get(type)
    UKBL = paste0(baseDataDir,'/UKBiLEVEAX_b1/sorted_bed/UKBiLEVEAX_b1-',type,'.bim') # NOTE: this assumes that all batches in UKBiLEVE have the same set
    UKBB = paste0(baseDataDir,'/Batch_b001/sorted_bed/Batch_b001-',type,'.bim') # NOTE: this assumes that all batches in UKBiobank have the same set of SNPs in their plink files!!! Which is true.
    
    print(UKBL)
    print(UKBB)

    UKBLsnps = read.table(UKBL,header=F,stringsAsFactors=F)
    UKBBsnps = read.table(UKBB,header=F,stringsAsFactors=F)
    
    snpsToInclude = Reduce(intersect, list(l$overlap,UKBLsnps$V2,UKBBsnps$V2))
    snpsToInclude = snpsToInclude[!snpsToInclude%in%c(l$affyFail,l$oxFail,lowFreq,l$imageSNPs)] # remove SNPs that failed at least once and have low frequency
    
    assign(paste0(type,".snpsToInclude"),snpsToInclude)

}


sexchrom2 = attach(sexchrom)
autosome2 = attach(autosome)

sexchrom2$lowFreq = sexchrom.lowFreq
autosome2$lowFreq = autosome.lowFreq
sexchrom2$snpsToInclude = sexchrom.snpsToInclude
autosome2$snpsToInclude = autosome.snpsToInclude

xyChroms2 = xyChroms
xyChroms2$lowFreq = intersect(sexchrom2$lowFreq, xySnps)
xyChroms2$snpsToInclude = intersect(sexchrom2$snpsToInclude, xySnps)


# threshold at 0.001 (doesn't work!)
sexchrom3 = attach(sexchrom)
autosome3 = attach(autosome)
sexchrom3$lowFreq = sexchrom.lowFreq
autosome3$lowFreq = autosome.lowFreq
sexchrom3$snpsToInclude = sexchrom.snpsToInclude
autosome3$snpsToInclude = autosome.snpsToInclude


length(intersect(autosome2$snpsToInclude,autosome$snpsToInclude))
length(intersect(autosome3$snpsToInclude,autosome$snpsToInclude))

extra = autosome2$snpsToInclude[!autosome2$snpsToInclude%in%autosome$snpsToInclude]
left = autosome$snpsToInclude[!autosome$snpsToInclude%in%autosome2$snpsToInclude]

extra = autosome3$snpsToInclude[!autosome3$snpsToInclude%in%autosome$snpsToInclude]
left = autosome$snpsToInclude[!autosome$snpsToInclude%in%autosome3$snpsToInclude]


hist(freqs$MAF[freqs$AffyID%in%extra],plot=F)



#################


#2. PLOT STATS BY SNP TYPES (not implemented)

# just need to read in one of the arrays, as we are only looking at overlapping SNPs anyway.
ps2snpBBsex = ukbiobank.ps2snp("sexchrom")
ps2snpBBaut = ukbiobank.ps2snp("autosome")

ps2snp = rbind(ps2snpBBsex,ps2snpBBaut)


                                        #3. PLOT STATS BY FAILURE TYPE AND BATCH
type="autosome"
snpqc = read.SNPQC.files(type=type,justSNPs=FALSE)
nSNPs = rep(dim(ps2snpBB)[1],length(batches))
ps2snpBoth = ps2snpBB[ps2snpBB$AffySNPID%in%intersect(ps2snpBB$AffySNPID,ps2snpBL$AffySNPID),]
names(nSNPs) = batches
nSNPs[ukbileve.batches()] = dim(ps2snpBL)[1]

Fracs = sapply(names(snpqc)[!sapply(snpqc,function(x) is.null(x))],function(y){
    print(y)
    test = snpqc[[y]]
    if( !y %in% c("hweSNPs","arraySNPs","concordanceSNPs","imageSNPs") ){
        
        nPerBatch = table(test[,1])   
        fracs = nPerBatch/nSNPs[names(nPerBatch)]
        
    }
    if(y=="imageSNPs") fracs = NULL#length(unique(test))/dim(ps2snpBoth)[1]
    if(y=="hweSNPs") fracs=NULL
    if(y %in% c("arraySNPs","concordanceSNPs") ) fracs=length(unique(test[,1]))/dim(ps2snpBoth)[1]

                                        #    out = cbind.data.frame(rep(y,length(fracs)),fracs)
    out = fracs
#    if(ncol(out)==3) out = out[,c(1,3)]

    return(out)
},simplify=FALSE)

fracPerBatch = abind(Fracs,along=1,hier.names=TRUE)
                                        # nSNPs per batch
toPlot = cbind.data.frame(fracPerBatch,gsub("\\.Batch.*$|\\.UKB.*$","",names(fracPerBatch)))

png('plots/boxplot-fractionSNPsPerTest.png',width=1200,height=1000,res=150)
par(mar=c(10,10,2,2))
e=factor(toPlot[,2],levels=unique(toPlot[,2])[order(sapply(unique(toPlot[,2]),median))])
boxplot(toPlot[,1]~e,las=2,ylab="Fraction of markers failing")

dev.off()


#4. PLOT STATS SNP BY GENOME REGION (windows of 100 MB)
l = total
uniqueFailuresAffy = unique(unlist(l$affyFail))
uniqueFailuresOxford = unique(unlist(l$oxFail))
uniqueRare = unique(unlist(l$lowFreq))

# for chunks 100MB compute the fraction of failed snps in each category (out of all overlapping SNPs)


windowSize=100000

DENS = list()
for(chrom in 1:26){
    print(chrom)
    
    snps = ps2snp$AffySNPID[ps2snp$Chromosome==chrom]
    overlap = l$overlap[l$overlap%in%snps]
    positions = ps2snp$Position[match(overlap,ps2snp$AffySNPID)]

    nWindows = round(( max(positions)- min(positions) )/windowSize,0)
    
    densities = sapply(1:nWindows,function(x){

        thesePos = c(min(positions) + (x-1)*windowSize,min(positions) + x*windowSize)
        midPoint = thesePos[1] + (thesePos[2] - thesePos[1])/2
        # include left-handside but not right-hand side
        s = ( positions >= thesePos[1] ) & ( positions < thesePos[2] )
        nsnps = sum(s)
        dens = c( sum(overlap[s]%in%uniqueFailuresAffy),sum(overlap[s]%in%uniqueFailuresOxford) )/nsnps
        
        return(c(dens,midPoint))
    })
    DENS[[chrom]] = densities
}


# genome-wide failure rate
failRateAffy = length(intersect(l$overlap,uniqueFailuresAffy))/length(l$overlap)
failRateOxford = length(intersect(l$overlap,uniqueFailuresOxford))/length(l$overlap)

for(chrom in 1:26){
    print(chrom)

    dens = DENS[[chrom]]
    
    x = dens[3,] # midpoint positions
    y1 = dens[2,] # Oxford failures
    y2 = dens[1,] # Affy failures
    
    colours = c("green","blue")
    labels = c("Proportion failed Oxford QC","Proportion failed Affy QC")

    png(paste0('plots/genome-local-failure-rates-',windowSize/1000000,'MB-window-chrom-',chrom,'.png'),width=1600,height=1000,res=150)
    plot(NULL,xlim=c(min(x),max(x)), ylim=c( min(c(y1,y2),na.rm=T), max(c(y1,y2),na.rm=T) ), ylab="Proportion of SNPs on both arrays",xlab=paste0("Position on chromosome ",chrom))
    lines(x,y1, col=colours[1])
    lines(x,y2, col=colours[2])
    legend("top", horiz=T,legend=labels,col=colours,inset=c(0,-0.1),lty=1,bty="n",xpd=NA)
    dev.off()

    
    png(paste0('plots/genome-local-failure-rates-',windowSize/1000000,'MB-window-chrom-',chrom,'-AffyOxford-scatter.png'),width=1600,height=1000,res=150)
    plot(dens[1,],dens[2,],xlab="Affy failure rate",ylab="Oxford failure rate")
    dev.off()

}



# plot affymetrix failures by batch
# failedAffyOnce



# plot some properties of included SNPs
    



