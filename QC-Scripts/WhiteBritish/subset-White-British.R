rm(list=ls( ))
args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/auxFunctions.R","-in","b1__b11-b001__b095-pca-UKbio-proj40")


print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}

setwd(paste0(baseSampleQCDir,"/QC-Scripts/WhiteBritish"))

library(aberrant)
## 'aberrant' lambda: the same value for lambda can lead to different 'separation' or 'strictness'
## criterion for identifying outliers. This seems to depend on how 'extreme' the outliers are. 
## The easiest solution to this problem is to run 'aberrant' on all points together.
lambda = 40


## I will assume that
## 1. PCA has been performed to compute SNP loadings
##    Some samples could have been excluded from this computation,
##    e.g. hetmiss outliers and related individuals.
## 2. Then all samples have been projected into the already computed
##    principal component space, including hetmiss outliers and 
##    related individuals.

BatchesFile = paste0(baseSampleQCDir,"/QC-Scripts/batchList.txt")
InputPCs =  args[which(args%in%c("-in"))+1]
OutputFile = InputPCs
batches = all.batches()

BatchList = read.table(BatchesFile,header = FALSE, stringsAsFactors = FALSE)
nBatches = nrow(BatchList)

Info = NULL

# Load in the PCs (computed using $basedir/QC-Scripts/PCA/submit-pca-UKBio.sh)
load(paste0(baseSampleQCDir,"/data/PCA/",InputPCs,".RData"),verbose=T)

fields = c("Country.of.birth..UK.elsewhere.",
    "Country.of.Birth..non.UK.origin.",
    "Ethnic.background",
    "Place.of.birth.in.UK...north.co.ordinate",
    "Place.of.birth.in.UK...east.co.ordinate")
Info = read.multiple.batch.info(fields)
Info = dplyr::tbl_df(Info)


# Add the Info to PCs
PCs = dplyr::left_join(PCs,Info,by=c("PIID"="PIID"))


### Function that runs aberrant and saves results. Takes as input a set of PCs subset by some ethnicity, as well as a name for the output files.

plotSubsets <- function(PCs.Subset,Indices.PCs.outliers,Name){
    
    Samples.PCs.outliers = PCs.Subset$PIID[Indices.PCs.outliers]
    Samples.PCs.inliers = PCs.Subset$PIID[-Indices.PCs.outliers]

    PCs.Subset$Colors2 = PCs.Subset$Colors
    PCs.Subset$Colors2[Indices.PCs.outliers] = "gray"

    
    if(length( Samples.PCs.outliers > 0 )){

        Order = order.by.number.occurrences(PCs.Subset$Colors2)
        PCs.Subset = PCs.Subset[Order,]
        PCs.Subset = PCs.Subset[c(which(PCs.Subset$Colors2=="gray"),which(PCs.Subset$Colors2!="gray")),]
        
        png(file=paste0("plots/",OutputFile,"-",Name,"-with_outliers-pc%02d.png"),
            width=8,height=7.5,units="in",res=150)
        plot(PCs.Subset$PC1,PCs.Subset$PC2,
             xlab="PC1",ylab="PC2",
             col=PCs.Subset$Colors2,pch=PCs.Subset$Chars)
        plot(PCs.Subset$PC3,PCs.Subset$PC4,
             xlab="PC3",ylab="PC4",
             col=PCs.Subset$Colors2,pch=PCs.Subset$Chars)
        plot(PCs.Subset$PC5,PCs.Subset$PC6,
             xlab="PC5",ylab="PC6",
             col=PCs.Subset$Colors2,pch=PCs.Subset$Chars)
        dev.off( )
    } else {
        print( paste0("WARNING: There are no outliers in ",Name,". Do you expect this?") )
    }
    
    if(length( Samples.PCs.inliers > 0 )){
        
        PCs.White = dplyr::filter(PCs.Subset, PIID %in% Samples.PCs.inliers)
        Order = order.by.number.occurrences(PCs.White$Colors2)
        PCs.White = PCs.White[Order,]
        
        png(file=paste0("plots/",OutputFile,"-",Name,"-wout_outliers-pc%02d.png"),
            width=8,height=7.5,units="in",res=150)
        plot(PCs.White$PC1,PCs.White$PC2,
             xlab="PC1",ylab="PC2",
             col=PCs.White$Colors,pch=PCs.White$Chars)
        plot(PCs.White$PC3,PCs.White$PC4,
             xlab="PC3",ylab="PC4",
             col=PCs.White$Colors,pch=PCs.White$Chars)
        plot(PCs.White$PC5,PCs.White$PC6,
             xlab="PC5",ylab="PC6",
             col=PCs.White$Colors,pch=PCs.White$Chars)
        dev.off( )        
        
    } else {
        print( paste0("WARNING: There are no inliers in ",Name,". Do you expect this?") )
    }

    # print list of all samples to include on this basis
    write.table(Samples.PCs.inliers,
                file=paste0(OutputFile,"-",Name,".txt"),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    # print list of all samples to exclude on this basis
    write.table(PCs$PIID[!PCs$PIID%in%Samples.PCs.inliers],
                file=paste0(OutputFile,"-",Name,"-exclusions.txt"),
                quote=FALSE,row.names=FALSE,col.names=FALSE)

}


runAberrantOnSubset <- function(PCs.Subset,Name){

## Unfortunately, 'aberrant' takes only two 'features'. So instead,
## run 'aberrant' three times -- this will take some time to finish
## Iterations = Default = 10,000. Takes about half a minute per ~500 iterations.

    print(paste0("Running Aberrant on ",Name ,"..." ))
    
    set.seed(123456) ## Set the seed for reproducibility. 
    aberrant.PCs12 = aberrant(data.frame(PC1=PCs.Subset$PC1,PC2=PCs.Subset$PC2),
        lambda,uncorr=FALSE,standardize=FALSE)
    aberrant.PCs34 = aberrant(data.frame(PC3=PCs.Subset$PC3,PC4=PCs.Subset$PC4),
        lambda,uncorr=FALSE,standardize=FALSE)
    aberrant.PCs56 = aberrant(data.frame(PC5=PCs.Subset$PC5,PC6=PCs.Subset$PC6),
        lambda,uncorr=FALSE,standardize=FALSE)
    ## Save the aberrant results, just in case
    save(PCs.Subset,aberrant.PCs12,aberrant.PCs34,aberrant.PCs56,
         file=paste0(OutputFile,"-",Name,".RData"))


    ## Use the leading 6 PCs to remove outliers in the subset of 
    ## self-declared British samples
    Indices.PCs.outliers = c(aberrant.PCs12$outlier)
    Indices.PCs.outliers = c(aberrant.PCs34$outlier,Indices.PCs.outliers)
    Indices.PCs.outliers = c(aberrant.PCs56$outlier,Indices.PCs.outliers)
    Indices.PCs.outliers = unique(Indices.PCs.outliers)


    ## This amounts to a fairly strict criterion for declaring a sample to be
    ## "White British"
    writeLines( paste0("The number of self-declared ",Name) )
    print(nrow(PCs.Subset))
    writeLines( paste0("The number of PCA outliers, among the self-declared ",Name) )
    print(length(Indices.PCs.outliers))


    print(paste0("Plotting results on ",Name ,"..." ))
    
    plotSubsets(PCs.Subset,Indices.PCs.outliers,Name)
    # just print the outliers for each set of PCs
    plotSubsets(PCs.Subset,aberrant.PCs12$outlier,paste0(Name,"-PC12"))
    plotSubsets(PCs.Subset,aberrant.PCs34$outlier,paste0(Name,"-PC34"))
    plotSubsets(PCs.Subset,aberrant.PCs56$outlier,paste0(Name,"-PC56"))

}


## Actually run this function.

## Choose the self-declared ethnicity criterion
##### TESTING ONLY
#PCs.Subset = dplyr::filter(PCs, Pops == "British")
#runAberrantOnSubset(PCs.Subset[1:10000,],"White_British_TEST") # test function on 10000 samples
#load(paste0(OutputFile,"-",Name,".RData"),verbose=T)
#####
#Name = "White_British"
#load(paste0(OutputFile,"-",Name,".RData"),verbose=T)
#Indices.PCs.outliers = c(aberrant.PCs12$outlier)
#Indices.PCs.outliers = c(aberrant.PCs34$outlier,Indices.PCs.outliers)
#Indices.PCs.outliers = c(aberrant.PCs56$outlier,Indices.PCs.outliers)
#Indices.PCs.outliers = unique(Indices.PCs.outliers)

#plotSubsets(PCs.Subset,Indices.PCs.outliers,Name)
                                        # just print the outliers for each set of PCs
#plotSubsets(PCs.Subset,aberrant.PCs12$outlier,paste0(Name,"-PC12"))
#plotSubsets(PCs.Subset,aberrant.PCs34$outlier,paste0(Name,"-PC34"))
#plotSubsets(PCs.Subset,aberrant.PCs56$outlier,paste0(Name,"-PC56"))

#quit()


## For the interim release, we chose "British" only
PCs.Subset = dplyr::filter(PCs, Pops == "British")
print( table(PCs.Subset$Pops) )
plotSubsets(PCs,which(PCs$Pops != "British"),"White_British_all") # just plot all the self-reported subset
runAberrantOnSubset(PCs.Subset,"White_British")


## Alternatively, choose samples with "Caucasian" background
thesePops = c("British","Irish","Any other white background")
PCs.Subset = dplyr::filter(PCs, Pops %in% thesePops)
print( table(PCs.Subset$Pops) )
plotSubsets(PCs,which(!PCs$Pops %in% thesePops),"Caucasian_all") # just plot all the self-reported subset
runAberrantOnSubset(PCs.Subset,"Caucasian")

## Alternatively, choose samples with "British/Irish/Scottish/Wales" background. Probably similar to "Caucasian", but only includes 'any other white' if they're born in UK + Ireland.
PCs.Subset = dplyr::filter(PCs, Pops == "British" | Pops == "Irish" | (Pops %in% c("Any other white background","Other/Unknown") & !PCs$Country.of.birth..UK.elsewhere. %in% c("Do_not_know","Elsewhere","Prefer_not_to_answer","") )  )
print( table(PCs.Subset$Pops) )
print( table(PCs.Subset$Country.of.birth..UK.elsewhere.) )
plotSubsets(PCs,which(!PCs$PIID%in%PCs.Subset$PIID),"Ireland_UK_all") # just plot all the self-reported subset
runAberrantOnSubset(PCs.Subset,"Ireland_UK")


## Use the same set of aberrant outliers to print a set of 'Irish/Scottish/Wales/England' samples. Will run a separate PCA on these samples, and even project on POBI PCs


print(paste0("DONE!"))
print(paste0("Files saved to: ",OutputFile," .txt"))
