#####################
# Script to plot posterior distributions and calculate some summaries.
#####################

#####################
# Preliminaries

args = commandArgs(trailingOnly=TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)
library(plyr)

# args = c("-posteriors","Posteriors-Standing.height-BOLT-LMM-v20-chr%%.out.chr%%.maf0.info0.pruned.giantRegs.0.2.RData,Posteriors-Standing.height-BOLT-LMM-v16.out.chr%%.maf0.miss1.pruned.giantRegs.0.2.RData,Posteriors-GIANT.chr%%.giantRegs.0.2.RData","-chr","1-22","-plotdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height","-outdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts")

sourceNames = c("IMP"="UK Biobank \n(imputed)","IMPSUB"="UK Biobank \n(imputed subset)","GENO"="UK Biobank \n(genotyped)","GIANT"="GIANT (2014)")

print(args)

title = args[which(args=="-title")+1]
dataFiles = str_split(args[which(args=="-posteriors")+1],",")[[1]]

if(length(dataFiles) != 3) {
    print("NEED 3 input files! Stopping.")
    quit()
}

outDir = args[which(args=="-outdir")+1]
print(outDir)

chroms = parse.range.string(args[which(args=="-chr")+1])


############################
##### Load in the data.

notGood=FALSE
first=TRUE
types=c()
for(file in dataFiles){

    if(grepl("v19|v20|IMP",file)) type="IMP"
    if(grepl("v19b|v20b|IMPSUB",file)) type="IMPSUB"
    if(grepl("v16|GENO",file)) type="GENO"
    if(grepl("GIANT",file)) type="GIANT"
    print(file)
    print(type)
    types = c(types,type)
    
    if(first) {
        load(paste0(outDir,"/",gsub("%%","22",file)))
        regionFileBase = regionFile
        first=FALSE
    }

     myPosteriors = lapply(chroms,function(chr){
         
         print(chr)
         
         f = paste0(outDir,"/",gsub("%%",chr,file))
         
         if(file.exists(f)){
             
             load(f,verbose=FALSE)
             
             if(regionFile != regionFileBase) {
                 print("ERROR: You have not specified posteriors using the same region files for each data source!")
                 quit()
             }
             
             return(list("DFsubset"=DFsubset,"regions"=regions))
             
         } else {
             print(paste0("Missing file for this chromosome: ",f))
             return(NULL)
         }
         
     })
         assign(paste0(type,".data"),myPosteriors)
 }

print(paste0("Processing regions for these data sources: ",paste0(types,collapse=",")))

outName = paste0("Posteriors-merged-",paste0(types,collapse="_"),".",title,".RData")
print(outName)

#quit()

############################
##### Find overlapping SNPS in the regions 

# First try:
# 1. Exclude any region that has no hits in the Imputed data.
# 2. Pick the top X non-overlapping regions in GIANT.

AllRegions = lapply(get(paste0("GIANT.data")),function(x) x$regions)

nChroms = length(AllRegions)
print(paste0("Processing regions in ",nChroms," chromosomes."))
print(paste0("Regions are from this file: ",regionFileBase))

matchFields = c('region.index','CHR','BP','SNP2')


MergedData = lapply(1:nChroms,function(chrRegions){

    print(paste0("Merging regions for the ",chrRegions,"th chromosome"))
    
    regs = AllRegions[[chrRegions]]

    myData = lapply(types,function(H){
        
        d = get(paste0(H,".data"))
        
        if(length(d)!=nChroms) {
            print("ERROR: your posterior lists are not aligned to the main set for source: ",H,".")
            quit()
        }
        dat = d[[chrRegions]]$DFsubset
        
                                        # Check that there are no region indices outside the set of regions.
        check1 = sum(!unique(dat$DFsubset$region.index)%in%1:nrow(regs))
        if(check1 > 0) print("ERROR: your posterior region lists have regions not included in the main set for source: ",H,".")

        #if(H!="GIANT") dat$SNP2.swap = dat$SNP2 # The swapping refers to GIANT, but this makes the merging easier later

        colnames(dat)[!colnames(dat)%in%matchFields] = paste0(H,"..",colnames(dat)[!colnames(dat)%in%matchFields])
                
        return(dat)
        
        
    })

    # Merge SNPs where possible (using chrom:BP_REF_ALT to match up the SNPs)
    # First merge non-GIANT data.
    mergedCounts1 = plyr::join_all(myData[types!="GIANT"], by=matchFields, type='full')
    giantData = myData[[which(types=="GIANT")]]

    toMatch1 = paste0(mergedCounts1$region.index,":",mergedCounts1$SNP2)
    toMatch2 = paste0(giantData$region.index,":",giantData$SNP2)
    toMatch3 = paste0(giantData$region.index,":",giantData$GIANT..SNP2.swap)

    genoIndexA = match(toMatch1,toMatch2)
    genoIndexB = match(toMatch1,toMatch3)

    genoIndex = genoIndexA; genoIndex[is.na(genoIndexA)] = genoIndexB[is.na(genoIndexA)]

    print(paste0(sum(!is.na(genoIndex))," snps (might not be unique) in IMP/GENO found in GIANT."))
    print(paste0(sum(!is.na(genoIndexB))," of these needed to be flipped in GIANT."))

    colsToAdd = !colnames(giantData)%in%c(matchFields)
    mergedCounts2 = cbind(mergedCounts1,giantData[genoIndex,colsToAdd])

    # Add any SNPs not found in GENO or IMP (shouln't be many)
    extras  = (!toMatch2%in%toMatch1)&(!toMatch3%in%toMatch1)
    print(paste0(sum(extras)," extra snps (might not be unique) in GIANT not found in IMP or GENO."))

    # Combine it all together
    mergedCounts3 = rbind.fill(mergedCounts2,giantData[extras,])

    # Which ones do we need to swap GIANT?
    mergedCounts3$GIANT.swapped = FALSE
    mergedCounts3$GIANT.swapped[which(mergedCounts3$SNP2==mergedCounts3$GIANT..SNP2.swap)] = TRUE
    if(sum(mergedCounts3$GIANT.swapped) != sum(!is.na(genoIndexB))) print("WARNING: Something wrong in merge...")
    

    # unique snp_region combinations:
    nUniq = unique(paste0(mergedCounts3$region.index,mergedCounts3$SNP2))
    nUniqBP = unique(paste0(mergedCounts3$region.index,mergedCounts3$BP))

    print(paste0("Number of unique region_SNP combinations: ",length(nUniq)))
    print(paste0("Rows in merged data frame: ",nrow(mergedCounts3)))

    if(nrow(mergedCounts3)!=length(nUniq)) print("WARNING: merging of SNPs in regions didn't work as expected.")

    postCols = paste0(types,"..Posterior")
    overlappingSNPs = !is.na(mergedCounts3[,postCols])
    colnames(overlappingSNPs) = paste0(types,"..hasit")

    print("Number of SNPs 1,2,3, datasets:")
    print(table(rowSums(overlappingSNPs)))

    return(mergedCounts3)
})



save(MergedData,AllRegions,chroms,regionFileBase,dataFiles,types,file=paste0(outDir,"/",outName))


print("DONE!")
print(outName)
