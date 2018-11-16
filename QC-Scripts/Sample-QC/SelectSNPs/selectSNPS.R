################################
# For UKBiobank sample QC pipeline step 1
################################
# Clare Bycroft, May 2016
#
# This script selects the SNPs for base plink files in sample QC. Lists of SNPs failing certain filters are created in SNP-QC pipeline (Colin Freeman), so we use these lists to exclude SNPs.
################################

args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

for(helperScript in args){
    source(helperScript)
}

batches = all.batches()

##### set filters #####

# 1. Passed Affy QC in all batches
# 2. Not Indels or multi-allelic
# 3. Overlap UKBiLEVE and UKBiobank 
# 4. Passed Oxford QC in all batches (p-value thresholds are set by SNP-QC pipeline, i.e Colin Freeman)
# 5. Not a very rare SNP (MAF < 0.0001)
# 6. SNPs in Affymetrix 37,000 SNPs with image problems


# select the rest by chromosome type

for( type in c("autosome","sexchrom") ){

    # preliminaries
    print(type)

    ps2snpBL = ukbileve.ps2snp(type)  # all snps in V2_QCed, see snpqc-tests.R
    ps2snpBB = ukbiobank.ps2snp(type)  # all snps in V2_QCed

    outputfile=paste0('snpIncludeList-',type,'.txt')

    
#1. FAILING AFFY
    
    # read in Colin's lists.
    affyFail = read.table(paste0(baseDataDir,'/SNP-QC_summary/V2_QCed.',type,'.Affy-SNPolisher-QC.fail.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt'),stringsAsFactors=FALSE,header=FALSE)[,1]

    print( paste0(length(unique(affyFail))," SNPs fail Affy QC in at least one batch (Colin's list)."))
    

#2. NON INDELS
    print('Getting list of SNPs (i.e not indels)')

    BLsnps = ps2snpBL$AffySNPID[ps2snpBL$Variant=="SNP"]
    BBsnps = ps2snpBB$AffySNPID[ps2snpBB$Variant=="SNP"]


#3. OVERLAP
    print('Getting list of SNPs on both arrays')
    overlap = intersect(BLsnps,BBsnps)
    print(paste0(length(unique(overlap)), " SNPs on both arrays") )
    
#4. PASSED OXFORD QC IN ALL BATCHES (for version v1 -- probably include some changes later)

    # read in Colin's list
    oxFail = read.table(paste0(baseDataDir,'/SNP-QC_summary/V2_QCed.',type,'.WTCHG.fail.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt'),stringsAsFactors=FALSE,header=FALSE)[,1]

    
#5. VERY RARE SNPs (near monomorphic - 1/10000)
    ##############
    lowMAF = 0.0001
    ##############
    
    # use genotype counts file created in SNP-QC
    rm("Counts")

    for (batch in batches){
#        print(batch)
            
        countsFile = paste0(baseDataDir,"/",batch,"/genotype_counts/",batch,"-",type,"-genotype_counts.txt")
        counts = read.table(countsFile,header=T,colClasses=c("character",rep("numeric",8),rep("NULL",16)))
        totals = counts[,c("AA_all_female","AB_all_female","BB_all_female","NC_all_female")] + counts[,c("AA_all_male","AB_all_male","BB_all_male","NC_all_male")] 

        totals[is.na(totals)] = 0  # set NAs to zero
        
        if("Counts"%in%ls()) Counts = totals + Counts else Counts = totals
    }
    
    rownames(Counts) = counts[,1]
    
#### THIS IS ACTUALLY WRONG (but I used it...)!!!    
    freq = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + Counts$AB_all_female)
    # SHOULD HAVE BEEN - with the 2x at the bottom:
#    freq = (2*Counts$AA_all_female + Counts$AB_all_female)/(2*Counts$AA_all_female + 2*Counts$BB_all_female + 2*Counts$AB_all_female)

    # This is not acutally only females. The names of the columns are just inherited
    lowFreq = rownames(Counts)[(freq <= lowMAF)|((1-freq) <= lowMAF)]
    print( paste0(length(unique(lowFreq))," SNPs with MAF lower than or equal to ",lowMAF ) )
    

#6. IMAGE PROBLEM SNPs FROM AFFY
    imageSNPs = read.table("probesets_in_subimage_37555.txt",header=TRUE,stringsAsFactors=FALSE)[,"affy_snp_id"]

    print( paste0(length(imageSNPs)," snps with image problem." ) )
    
#7. INTERSECT WITH SNPs IN PLINK FILES
    print('Getting intersection with SNPs in original plink files...')

    UKBL = paste0(baseDataDir,'/UKBiLEVEAX_b1/sorted_bed/UKBiLEVEAX_b1-',type,'.bim') # NOTE: this assumes that all batches in UKBiLEVE have the same set
    UKBB = paste0(baseDataDir,'/Batch_b001/sorted_bed/Batch_b001-',type,'.bim') # NOTE: this assumes that all batches in UKBiobank have the same set of SNPs in their plink files!!! Which is true.

    print(UKBL)
    print(UKBB)

    UKBLsnps = read.table(UKBL,header=F,stringsAsFactors=F)
    UKBBsnps = read.table(UKBL,header=F,stringsAsFactors=F)# aggh! actually only read in UKBiLEVE SNPS! but doesn't matter, because anything on UKBiobank only won't be kept (i.e only the intersect of the two...)

# check that these snps are indeed a subset of the snps on the whole array (i.e not any rsIDs involved)
    allsnps = c(ps2snpBL$AffySNPID,ps2snpBB$AffySNPID)

    if( sum(!UKBLsnps$V2%in%allsnps) != 0 ) print("WARNING: not all snps in .bim files are in ps2snp UKBiLEVE! Something wrong with SNP ID convention in plink files?")
    if( sum(!UKBBsnps$V2%in%allsnps) != 0 ) print("WARNING: not all snps in .bim files are in ps2snp UKBiobank! Something wrong with SNP ID convention in plink files?")

    snpsToInclude = Reduce(intersect, list(overlap,UKBLsnps$V2,UKBBsnps$V2))

    if(sum(!snpsToInclude%in%affyFail)==0) print("WARNING: no SNPs in plink files failed Affy once... suspiciously wrong!")

#    snpsToInclude = snpsToInclude[!snpsToInclude%in%c(unlist(failedAffyOnce),unlist(failedOxfordOnce),lowFreq)] # remove SNPs that failed at least once

    snpsToInclude = snpsToInclude[!snpsToInclude%in%c(affyFail,oxFail,lowFreq,imageSNPs)] # remove SNPs that failed at least once and have low frequency

    print(paste0(length(snpsToInclude),' snps included in list'))


#7. OUTPUT TO FILE

    print(paste0('writing list of SNPs to output file ',outputfile))
    write.table(snpsToInclude,file=outputfile,quote=F,row.names=F,col.names=F)


#8. SAVE RDATA FOR PLOTTING
    save(snpsToInclude,overlap,oxFail,affyFail,lowFreq,lowMAF,imageSNPs,file = paste0('snpIncludeList-',type,'.Rdata'))

}




