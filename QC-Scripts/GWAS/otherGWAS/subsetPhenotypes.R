#############
# This script takes a list of phenotype fields, an optional file with sample IDs, and returns a nice text file for running SNP-test or BOLT, or any other test for that matter.
# It uses the extraction script created by W. Rayner (wrayner@well.ox.ac.uk), A. Mahajan (anubha@well.ox.ac.uk) 2015
#############

args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

h=args
for(s in h){
    source(s)
}


setwd('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS')

#############
vers="v0"  # All samples!
#############
#############
vers="v1"  # Using data from preliminary final run. PCs init.
#############
#############
vers="v2"  # Using data from preliminary final run. White British subset (from testing PCA). Create quantile normalised version of each trait as well. PCs round2
#############
vers="v3"  # As v2 but using PCs from White British subset run
#############
vers="v11"  # As v3 but special for X,Y & MT
#############
vers="v17"  # Same samples as v3 but random vectors.
#############


#### %%%%%%%%% ####
#### If >= v3 then skip all between %%%%%%%%%

genoInfo = read.multiple.batch.info(fields = batchInfoFields)

# get Array
genoInfo$Array = "UKBB"
genoInfo$Array[genoInfo$Batch%in%ukbileve.batches()] = "UKBL"

genoInfo$Array.as.binary = 0
genoInfo$Array.as.binary[genoInfo$Batch%in%ukbileve.batches()] = 1

# get Y:X ratio
genoInfo$Y.X.intensity.ratio = genoInfo$Y.intensity/genoInfo$X.intensity

# get averate IOP
genoInfo$Intra.ocular.pressure..Goldmann.correlated..mean = (genoInfo$Intra.ocular.pressure..Goldmann.correlated..left. + genoInfo$Intra.ocular.pressure..Goldmann.correlated..right.)/2
genoInfo$Intra.ocular.pressure..corneal.compensated..mean = (genoInfo$Intra.ocular.pressure..corneal.compensated..left. + genoInfo$Intra.ocular.pressure..corneal.compensated..right.)/2

# read in data for het/missing
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"),verbose=T)
hetmisslist = read.table(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt"),header=FALSE,stringsAsFactors=FALSE)[,1]

# read in data for excessive relatedness
relexcess = read.table(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-excess-relatives.txt"),header=FALSE,stringsAsFactors=FALSE)[,1]

# read in data for PCs (%%%%%% CHANGE THIS TO FULL SET LATER %%%%%%%%)
#load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData"),verbose=T)
if(vers %in% c("v0","v1","v2") ) load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-round2.RData"),verbose=T)
if(vers %in% c("v3") ) load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-White_British.RData"),verbose=T)

# fold in the PCs
genoInfo = tbl_df(genoInfo)
genoInfo = left_join(genoInfo,PCs,by=c("PIID"="PIID"))
# only include PCs 1-20
#genoInfo = genoInfo[,!colnames(genoInfo)%in%paste0("PC",21:60)]

# get 'Caucasian' list
CaucasianFile ="/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt"
ancestrySubset = read.table(CaucasianFile,header=FALSE)[,1]


###### get list of samples to use in GWAS
# Keep just males, (start with White British set, too) and exclude any related samples

# get list of QC'd samples
## NOTE: this is preliminary as of May 14th. I.e haven't run aberrant. Version 2 uses aberrant outliers list.
hetMiss = Table
highHet = hetMiss$IID[hetMiss$het.corrected > 0.2]
highMiss = hetMiss$IID[hetMiss$logit.miss > -2.75]
sexMismatch = genoInfo$PIID[genoInfo$Submitted.Gender!=genoInfo$Inferred.Gender]

if(vers=="v0") {
    genoInfoGood = genoInfo
} else {

                                        # do some exclusions
    
    if(vers=="v1") excludeSamples = unique(c(highHet,highMiss,sexMismatch,other))

    if(vers%in%c("v2")) {
        noncaucasian = genoInfo$PIID[!genoInfo$PIID%in%ancestrySubset] # non white british
        excludeSamples = unique(c(noncaucasian,relexcess,hetmisslist))
                                        # 344,397 after excluding samples for ethnic background (78,790), QC (1,167), and relatedness (65,128)
    }

    print( paste0(length(excludeSamples)," samples excluded for QC reasons, before finding unrelated set." ) )

    write.table(excludeSamples,file="excludeSamples.tmp",quote=FALSE,row.names=FALSE,col.names=FALSE)

                                        # get unrelated samples within the QC'd set (this script also finds individuals with many relatives)
    args <- c(h,"-in",paste0(baseSampleQCDir,"/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0"), "-ibs0","0.0012", "-exclude", "excludeSamples.tmp","-outdir",".","-outfile",paste0("b1__b11-b001__b095-pair_batches.filtered-",vers) )

    source( paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness/postking/get-unrelated.R")) 


                                        # subset genoInfo data
    goodSamples = read.table(paste0("b1__b11-b001__b095-pair_batches.filtered-",vers,"-unrelated.txt"),header=FALSE,stringsAsFactors=FALSE)[,1]

    genoInfoGood = filter(genoInfo,genoInfo$PIID%in%goodSamples)

## V1: 412,781 samples (all ethnicities)
## V2: 344,397 samples (Caucasian...etc.)

}


# Quantile normalise the quantitative traits (where samples are non-missing); and plot the histograms

quantTraits = colnames(genoInfoGood)[sapply(colnames(genoInfoGood),FUN=function(x) class(genoInfoGood[[x]])%in%c("numeric","integer") )]
quantTraits = quantTraits[!quantTraits%in%c("Reproducibility","Concordance","Encoded.anonymised.participant.ID","EID","Chars.x","Chars.y","Array.as.binary",grep("PC",quantTraits,value=TRUE) )]

system('mkdir plots')
for(x in quantTraits){
    print(x)
    # for X.Y.intensity, height, and BMI, then quantile normalise within sexes
    
    X = genoInfoGood[[x]]
    if(sum(!is.na(X))==0) next

    if(x %in% c(grep("intensity",quantTraits,value=TRUE),"Body.mass.index..BMI.","Standing.height") ){
        genoInfoGood$this = NA
        f = genoInfoGood$Inferred.Gender=="F"
        m = genoInfoGood$Inferred.Gender=="M"
        X2 = quantile.norm(X[f])   
        genoInfoGood$this[f] = X2
        X2 = quantile.norm(X[m])   
        genoInfoGood$this[m] = X2

    } else {
        X2 = quantile.norm(X)   
        genoInfoGood$this = X2
    }
    colnames(genoInfoGood)[ncol(genoInfoGood)] = paste0(x,".qnorm")

    nonmissing = sum(!is.na(X2))
    png(paste0("plots/",x,"-",vers,"-histogram.png"),height=500,width=500,res=100)
    hist(X,breaks=100,xlab=x,main=paste0(nonmissing," non-missing samples"))
    dev.off()
    png(paste0("plots/",x,"-",vers,"-qnorm-histogram.png"),height=500,width=500,res=100)
    hist(X2,breaks=100,xlab=x,main=paste0(nonmissing," non-missing samples"))
    dev.off()

}


# Now write phenotype file for plink (or BOLT-LMM). E.g continuous vs categorical; missing value code etc.
BOLTCols = c("PIID","PIID",grep("PC",colnames(genoInfo),value=TRUE),batchInfoFields[c(1:8,10,11,17,19,20,23,24,28,29,32:52,54,56)],"Y.X.intensity.ratio","Array","Array.as.binary","Intra.ocular.pressure..Goldmann.correlated..mean","Intra.ocular.pressure..corneal.compensated..mean",grep("qnorm",colnames(genoInfoGood),value=TRUE))

toPrint = genoInfoGood[,BOLTCols]
colnames(toPrint)[1:2] = c("FID","IID")

# Write any blank spaces as NA!!!
toPrint[toPrint==""] = NA

p1 = sapply(colnames(toPrint),function(x) sum(toPrint[[x]] < 0,na.rm=TRUE))
p2 = sapply(colnames(toPrint),function(x) sum(is.na(toPrint[[x]])))
p3 = sapply(colnames(toPrint),function(x) sum(grepl(" ",toPrint[[x]]),na.rm=TRUE))
p4 = sapply(colnames(toPrint),function(x) sum(toPrint[[x]]=="",na.rm=TRUE) )

write.table(toPrint,file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

# can we read it? i.e are the black spaces sorted out?
test = read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt"),header=TRUE)


#### %%%%%%%%% ####

if(vers %in% c("v3")){
    
    genoInfoGood = read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v2.txt"),header=TRUE)
    genoInfoGood = tbl_df(genoInfoGood)
    genoInfoGood = genoInfoGood[,!grepl("PC",colnames(genoInfoGood))]
    
    # replace PCs
    load(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-pca-UKbio-White_British.RData"),verbose=T)
    PCs = PCs[,c("PIID",grep("PC",colnames(PCs),value=TRUE))]
    genoInfoGood2 = left_join(genoInfoGood,PCs,by=c("IID"="PIID"))
    
    toPrint = genoInfoGood2
    
    write.table(toPrint,file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

}

if(vers %in% c("v11")){

    # For each of X (23), PAR (25), Y (24) and MT (26) create different phenotype files based on v3
    genoInfoGood = read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt"),header=TRUE)
    genoInfoGood = tbl_df(genoInfoGood)

    # get list of odd sex chrom samples (same list as used in phasing - made using create-exclusion-flats.R)
    sampleFilterFile = paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Flags/b1__b11-b001__b095-sampleQC-flags.txt")
    exclSex = read.table(sampleFilterFile,header=TRUE,stringsAsFactors=FALSE)
    exclSexSamples = exclSex$PIID[exclSex$possible.sex.aneuploidy==1]
    # check that none of the het-miss outliers etc. are actually in the file! As we used slightly different files to generate version v3. (yep, this is correct)
    
    length(exclSexSamples) # 652 samples; 
    sum(exclSexSamples%in%genoInfoGood$IID)  # 458 are in genoInfoGood (so will be removed on X and Y chroms)
    
    for( chr in 23:26){
        
        print(chr)
        
        if( chr == 26 ) toPrint = genoInfoGood  # MT just keep the same set of samples
            
        if( chr %in% c(23,25) ) {
            # just exclude odd sex chrom samples on X
            toPrint = genoInfoGood[!genoInfoGood$IID%in%exclSexSamples,]
        }
        if( chr == 24) {
            # Exclude odd sex chrom samples on X + all inferred females!
            toPrint = genoInfoGood[(!genoInfoGood$IID%in%exclSexSamples)&(genoInfoGood$Inferred.Gender=="M"),]
                                        # 146 males excluded here; all females (obviously)
        }
        
        print( nrow(toPrint) ) 
        
        write.table(toPrint,file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,"-",chr,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

    }    

}

if(vers %in% c("v17")){
                                        # construct a set of random phenotypes ~N(0,1)
    genoInfoGood = read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt"),header=TRUE)
    genoInfoGood = tbl_df(genoInfoGood)

    set.seed(123456)
    Rnorm.0.1 = rnorm(n=dim(toPrint)[1],mean=0,sd=1)
    Rlnorm.0.1 = rlnorm(n=dim(toPrint)[1],meanlog=0,sdlog=1)
    
    toPrint = cbind(genoInfoGood[,c("FID","IID")],Rnorm.0.1,Rlnorm.0.1)
    write.table(toPrint,file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)


}


#toPrint = test

##############################################################################
########## SUBSETS FOR CREDIBLE SET ANALYSIS (RANDOM ELEMENT HERE!! USING RANDOM SEED)
vers = "v19" # <== check runBOLT_v19-LMM.sh for which phenotype file was used in this version.
pheno = "Standing.height"
toPrint=read.table(paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt"),header=TRUE)
imputedSample = read.table("/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample",header=FALSE,skip=2)

set.seed(123456)
subVersions = paste0(vers,letters[1:4])
N = c(1000,10000,50000,150000)
names(N) = subVersions

noPheno = is.na(toPrint[,pheno])| !toPrint$IID %in% imputedSample[,1]

for(subVersion in subVersions){
    print(subVersion)
    n = min(N[subVersion],nrow(toPrint[!noPheno,]))

    subset = sample(toPrint$IID[!noPheno],n,replace=FALSE)
                                        # actually want a list that excludes samples from the whole dataset
    toPrintExclusions = imputedSample[,1][!imputedSample[,1] %in% subset]
    length(toPrintExclusions)

    write.table(cbind(toPrintExclusions,toPrintExclusions),file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/SampleExclusionsForBOLT-",subVersion,".txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)

}


########## SUBSETS FOR ILLUSTRATION (RANDOM ELEMENT HERE!! USING RANDOM SEED)
set.seed(123456)
subVersions = paste0(vers,letters[1:4])
N = c(1000,10000,50000,150000)
names(N) = subVersions
for(subVersion in subVersions){
    print(subVersion)
    UKB = toPrint$Batch%in%ukbiobank.batches() # just sample from UKBiobank where possible
    n = min(N[subVersion],nrow(toPrint))
    subset = sample(toPrint$IID[UKB],n,replace=FALSE)
    # actually want a list that excludes samples from the whole dataset
    toPrintExclusions = genoInfo$PIID[!genoInfo$PIID %in% subset]
    length(toPrintExclusions)
    # are they all in the UKBiobank?
    print( table(genoInfo$Batch[!genoInfo$PIID%in%toPrintExclusions]) )
    write.table(cbind(toPrintExclusions,toPrintExclusions),file=paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/SampleExclusionsForBOLT-",subVersion,".txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
}



####### FOR SNPTEST: BELOW NOT DONE YET #######


# Create files for snp.test. They must have the same set of samples as in the plink files (just put them in the same order for each chromosome, for safety!)
toPrint = genoInfo[,BOLTCols] # must include everyone. Add extra lines with NA for extra samples (they'll be excluded anyway)
fam = read.table("/well/donnelly/ukbiobank_project_8874/interim_genotype_data/binaryPed/ukb8874.chr22.calls.fam",header=FALSE)
extraSamples =fam[,1][!fam[,1]%in%toPrint[,1]]
extra = toPrint[1:3,]; extra[,1] = extraSamples
toPrint = rbind(toPrint,extra)

toPrint$missing = 0
colnames(toPrint)[1:2] = c("ID_1","ID_2")
toPrint = toPrint[,c(1,2,ncol(toPrint),3:(ncol(toPrint)-1))]
#extraHeader = c(0,0,0,"B",rep("C",15),"P","P","B","C","D") version 2 only
extraHeader = c(0,0,0,"B",rep("C",15),"P","P","B","C","D","P")
cbind(extraHeader,colnames(toPrint))

rownames(toPrint) = toPrint$ID_1
# make separate phenotype file for each chromosome, in order of .fam file
Order22 = as.character(fam[,1])

for (chr in 1:22){
    print(chr)
    f = read.table(paste0("/well/donnelly/ukbiobank_project_8874/interim_genotype_data/binaryPed/ukb8874.chr",chr,".calls.fam"),header=FALSE)
    Order = as.character(f[,1])
    toPrintThisChrom = toPrint[Order,]
    print(sum(toPrintThisChrom$ID_1!=Order))
    print(sum(toPrintThisChrom$ID_1!=Order22))
    
    File = paste0("data/YchromPhenotypesForSNPTest-",vers,"-chrom",chr,".sample")
                                        #File = paste0("data/YchromPhenotypesForSNPTest-",vers,".sample")
    write.table(t(colnames(toPrint)),file=File,col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(t(extraHeader),file=File,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(toPrintThisChrom,file=File,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

exclusions = fam[,1][!fam[,1]%in%genoInfoGood[,1]]

write.table(exclusions,file=paste0("data/YchromPhenotypesForSNPTest-",vers,".exclusions"),col.names=FALSE,row.names=FALSE,quote=FALSE)


# ====> now run BOLT, or snptest
