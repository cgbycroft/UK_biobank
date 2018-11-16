#######
# script to get hit regions and compute LD. Normally this is done in the region.plots.new.R script, but here we have manually selected regions for given phenotype
#######

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

BB.ps2snp = rbind(ukbiobank.ps2snp("autosome"),ukbiobank.ps2snp("sexchrom"))
BL.ps2snp = rbind(ukbileve.ps2snp("autosome"),ukbileve.ps2snp("sexchrom"))
QCexclude=c()

args = commandArgs(trailingOnly=TRUE)

genotypeFile = "/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095"


#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v13/Standing.height-BOLT-LMM-v13-chr2.out","2","plots","-ymax","50","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgen","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr2.hrc+uk10k.I4.v1.1.bgen","-sample","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr2.hrc+uk10k.I4.v1.1.sample","-regions","Standing.height-BOLT-LMM-v14.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15/Standing.height-BOLT-LMM-v15-chr2.out","2","plots","-ymax","50","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-lreg","-title","-raw","-bgen","/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr2.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen","-sample","/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample","-regions","Standing.height-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v16/Internal.Pico..ng.uL.-BOLT-LMM-v16.out","1-22", "plots" ,"-ymax", "50", "-qc", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-QCFiltered" ,"-regions", "../Internal.Pico..ng.uL..qnorm/Internal.Pico..ng.uL..qnorm-BOLT-LMM-v16.out.chrgenome.maf0.miss1.pruned-raw-lreg-hitsAndLDcalc.RData")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v15s/Standing.height-BOLT-LMM-v15s-chr23.out","23", "plots" ,"-ymax", "50", "-qc", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-raw","-bgen","/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chrX.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen","-sample","/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample","-regions", "Standing.height-BOLT-LMM-v16s.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-callthreshold","0.3")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr23.out","23", "plots" ,"-ymax", "50", "-qc", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-raw","-bgen","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrX.v4.v1_1.bgen","-sample","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.v2.sample","-regions", "Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-callthreshold","0.1")

# PAR v2 bgen
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr25.out","25", "plots" ,"-ymax", "50", "-qc", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-QCFiltered","-bgen","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR1.v4.v1_1.bgen","-sample","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrPAR1.v2.sample","-regions", "Standing.height-BOLT-LMM-v16s-25.out.chr25.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-callthreshold","0.1")

# X bgen
#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr23.out","23", "plots" ,"-ymax", "50", "-qc", "-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-lreg", "-title", "-QCFiltered","-bgen","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrX.v4.v1_1.bgen","-sample","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.v2.sample","-regions", "Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData","-callthreshold","0.1")


print(args)
print(getwd())

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

chroms = parse.range.string(chroms)
if((length(chroms)==1)&(!grepl("chr",dataFile))){
    if(chroms%in%c(23:26,"X","XY","Y","MT")){
        dataFile = paste0(dirname(dataFile),"/",gsub(".out",paste0("-",chroms,".out"),basename(dataFile)))
        args = c(args,"-lreg")
    }    
}

# any specified bgen or genotypes file for computing LD?
#bgenData=args[which(args=="-bgenFile")+1]

# look at lreg or lmm results?
if("-lreg" %in% args ) {
    useLmmInf = FALSE
    extraTitle = paste0(extraTitle,"-lreg")
} else {
    useLmmInf=TRUE
}

if("-bgen"%in%args){
    bgenFile=args[which(args=="-bgen")+1]
    bgenSampleFile=args[which(args=="-sample")+1]
#    if(grepl("chrX",bgenFile)){
                                        #        bgenSampleFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.sample"
 #       bgenSampleFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chrX.v2.sample"
 
        print("Using the corresponding sample file for the sex chromosome.")
   # }
    print(bgenSampleFile)
    print(bgenFile)
#    if(grepl("PAR1",bgenFile)) names(sexChroms)[names(sexChroms)=="PAR"]="PAR1"
#    if(grepl("PAR2",bgenFile)) names(sexChroms)[names(sexChroms)=="PAR"]="PAR2"
    print(sexChroms)
}

                                        # What callthreshold are we using for bgen to plink conversion for ld stats?
if("-callthreshold"%in%args) {
    callT = as.numeric(args[which(args=="-callthreshold")+1])
    extraTitle = paste0(extraTitle,"-ct",callT)
} else {
    callT=0.1
}
# 0.1 is the plink default: https://www.cog-genomics.org/plink/1.9/input#oxford

# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

###  Set QC thresholds!
if("-qc" %in% args){
                                        # QC THRESHOLDS  (used in plots made before 6-07-2016)
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


#### Define regions to look for hits in.
regionsToLookFile = args[which(args=="-regions")+1]

if( grepl("RData",regionsToLookFile) ) {
    
    # load in results from region.plost.new.R (i.e LDchrom; HITSchrom)
    if(regionsToLookFile%in%list.files()) {
        load(regionsToLookFile,verbose=TRUE)
    } else {
        regionsToLookFile = gsub(".out.chr.*.maf",paste0("-",chroms[1],".out.","chr",chroms[1],".maf"),regionsToLookFile)
        load(regionsToLookFile,verbose=TRUE)
    }
    
    HITSchromOrig = HITSchrom[names(HITSchrom)!=""]
    regionsToLook = sapply(names(HITSchromOrig),function(i) {
#        print(i)
        r = HITSchromOrig[[i]]$regs
        nRegions = dim(r)[1]; if(is.null(nRegions)) return(NULL)
        cbind(as.numeric(rep(i,nRegions)),r)
    },simplify=FALSE)
    regionsToLook = abind(regionsToLook,along=1)

    onlyFindTopSNP=TRUE  # i.e keep the regions as they are, but just find the top SNP within the region based on imputed data.
    
} else {
    regionsToLook = read.table(regionsToLookFile,header=FALSE,stringsAsFactors=FALSE)
    
    if("-fixRegions"%in%args) {
        # Fix the regions to look for the highest snp (i.e don't combine regions that overlap.)
        onlyFindTopSNP=TRUE
        HITSchromOrig = vector("list", 22); names(HITSchromOrig)=as.character(1:22)

        for( i in unique(regionsToLook[,1]) ){
            re = regionsToLook[regionsToLook[,1]==i,-1]
            colnames(re)=c("starts","ends")
            HITSchromOrig[[as.character(i)]]$regs = re
        }
        
    } else {
        onlyFindTopSNP=FALSE
    }
}


chroms = intersect(chroms,unique(regionsToLook[,1])) # only look in chromosomes that are both in the regions to look, and the chromosomes asked for!

if(length(chroms)==0) {
    print("No regions in this chromosome. Stopping.")
    quit()
}

###################################a
# read GWAS results! (GENO doesn't necessarily mean genotypes in this context.)

resultsGENO = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
DF1 = resultsGENO$DF

#resultsIMP = read.gwas(bgenData,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle,useLmmInf=useLmmInf)
#DF2 = resultsIMP$DF


### read in recombination rates
print('Reading recombination rate files from /well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr*')

recombrates = read.recomb(chroms)


###################################a
# Find hit regions, compute LD and save data.

HITSchrom = vector("list", 22); names(HITSchrom)=as.character(1:22)
LDchrom = vector("list", 22); names(LDchrom)=as.character(1:22)
computeLD=TRUE
samplesFile = args[which(args=="-phenoFile")+1]
system(paste0("cut -f1,2 -d' ' ",samplesFile," | tail -n +2 > ",samplesFile,".samples.tmp"))
saveData=TRUE


for (chrom in chroms){
#for (chrom in 5:22){

    print(chrom)
    if(chrom%in%c(23:25)) chromName = names(sexChroms)[sexChroms==chrom] else chromName=chrom
    genotypeFile = paste0("/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr",chromName)
    print(genotypeFile)
    
    # any regions to look in this chromosome?
    if(! chrom%in%regionsToLook[,1] ) {
        LDchrom[[as.character(chrom)]]=NULL
        HITSchrom[[as.character(chrom)]]=NULL
        next
    }
    
    Pvalset = paste0(gsub("genome",chrom,resultsGENO$Pvalset))
    
#### Find hit regions automatically within these regions

    snpsToCheck = unique(unlist(apply( regionsToLook[regionsToLook[,1]==chrom,,drop=FALSE],1,function(x){
        s1 = as.numeric(x[2]);s2 = as.numeric(x[3]);
        #print(s1)
        #print(s2)
        keep = which((resultsGENO$DF$CHR==chrom) & (resultsGENO$DF$BP >= s1) & (resultsGENO$DF$BP <= s2))
    })))

    print(length(snpsToCheck))
    
    recomb = recombrates[[as.character(chrom)]]
    
    # NOTE: If hit regions are already provided, this finds the top SNP within that region.
    if( onlyFindTopSNP ){
        
        HITS = HITSchromOrig[[as.character(chrom)]]
            
        DFchrom = resultsGENO$DF[snpsToCheck,]
        hits = (DFchrom$P < 5e-8) & (DFchrom$CHR == chrom)
        if(sum(hits)==0) hits = (DFchrom$CHR == chrom)
        DFthese = DFchrom[hits,]

        HITS$DFthese = DFthese
        
        HITS$topVariants = apply(HITS$regs,1,function(x){
            
            s1 = as.numeric(x[1]);s2 = as.numeric(x[2]);
            inReg = (DFthese$BP >= s1)&(DFthese$BP <= s2)
            topSNP = DFthese$SNP[inReg][order(DFthese$P[inReg])][1]

        })
        HITS$topVariants2 = apply(HITS$regs,1,function(x){
            
            s1 = as.numeric(x[1]);s2 = as.numeric(x[2]);
            inReg = (DFthese$BP >= s1)&(DFthese$BP <= s2)
                                        # if this is the imputed data, then use alternative IDs.
            topSNP = DFthese$SNP2[inReg][order(DFthese$P[inReg])][1]
        })
        
            
    } else {

        print("looking for hit regions first...")
        HITS = find.hits.recomb(resultsGENO$DF[snpsToCheck,],chrom,minGap=0.125,buffer=25,recombRates=recomb)
        
    }
    
    HITSchrom[[as.character(chrom)]] = HITS

#### compute LD with top SNPs (if they exist, and were asked for)
    if( length(HITS$topVariants) > 0 ){
        
        if(computeLD){
            if( "-bgen" %in% args){
               #if(chrom==23) sex="F" else sex =NULL # only compute LD on females if chrom 23
                sex = NULL # <=== males and females together regardless.
                print("Computing LD based on imputed data...")
                print(bgenFile)
                print(bgenSampleFile)

                # use alternative SNP IDs
                ld = compute.ld.around.snp.bgen(variants=HITS$topVariants2,bgenFile,bgenSampleFile,chrom,paste0(samplesFile,".samples.tmp"),regions=HITS$regs,DF1,callthreshold=callT,sex=sex)
#                ld = compute.ld.around.snp.bgen(variants=HITS$topVariants2[c(7,44)],bgenFile,bgenSampleFile,chrom,paste0(samplesFile,".samples.tmp"),regions=HITS$regs[c(7,44),],DF1,callthreshold=callT,sex=sex,rNumber=34542)
                

            } else {

                ld = compute.ld.around.snp(HITS$topVariants,genotypeFile,chrom,paste0(samplesFile,".samples.tmp"),bgenSampleFile,DF1)
                
            }
            LDchrom[[as.character(chrom)]] = ld
        }
    } else {
        LDchrom[[as.character(chrom)]]=NULL
        print( paste0("No hits found for chromosome ",chrom) )
    }   
}                


save( LDchrom,HITSchrom,file=paste0(resultsGENO$Pvalset,"-hitsAndLDcalc.RData") )
# ====> now you can plot these regions as you like!


quit()

## check for duplicates (chrom 22)
DFraw = resultsGENO22$DFraw
isdup = duplicated(DFraw$SNP)
isdupos = DFraw$BP[isdup]
dupsDFraw = DFraw[DFraw$BP%in%isdupos,]


write.table(unique(dupsDFraw$SNP),file="chr22.v13.gwas.duplicate.SNPs.txt",quote=F,col.names=F,row.names=F)

system1 = paste0(bgenix," -g ../../../../data/imputed/chr22.hrc+uk10k.I4.v1.1.bgen -incl-rsids chr22.v13.gwas.duplicate.SNPs.txt > chr22.v13.gwas.duplicate.SNPs.bgen")
system(system1)

system1 = paste0(qctool," -g chr22.v13.gwas.duplicate.SNPs.bgen -snp-stats chr22.v13.gwas.duplicate.SNPs.snp-stats")
system(system1)


# read snp-stats
dupStats = read.delim("chr22.v13.gwas.duplicate.SNPs.snp-stats",header=TRUE,stringsAsFactors=FALSE,sep="\t")

tab=table(dupStats$position)
table(tab)

#   3    4 
#1159    4

# 1159 positions are there three times; 4 positions are there 4 times.

dupStats[dupStats$position%in%as.numeric(names(which(tab==4))),]

                                        # count alleles
posdup=unique(dupStats$position)
alleleCounts = sapply(posdup,function(p){
    as = dupStats[dupStats$position==p,c("A_allele","B_allele")]
    length(unique(c(as[,1],as[,2])))
})

twoAs = posdup[alleleCounts==2]
check = dupStats[dupStats$position%in%twoAs,]
 # are these exceptions all indels?
indelsCheck = apply(check[,c("A_allele","B_allele")],1,function(i){
    max(c(nchar(i[1]),nchar(i[2])))
})

min(indelsCheck) # = 2  ===> all are indels. 


write.table(check[,-1],file="chr22.v13.gwas.duplicate.SNPs.snp-stats.INDELS.txt",quote=F,col.names=T,row.names=F,sep="\t")
