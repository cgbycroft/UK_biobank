source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


genotypeFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc"

#print("Reading QC snps lists...")
#QCSNPList = read.SNPQC.files(justSNPs=TRUE)

#QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$concordanceSNPs))                                            # NOTE: anything created after Friday June 24th uses this QCexclude. Anything created before that also excluded the 'image' SNPs    
QCexclude = c()

##########

args = commandArgs(TRUE)

#args = c("test.out","1","plots", "-ymax","50","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-Euro-hits")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v3/Standing.height-BOLT-LMM-v3.out","22","plots", "-ymax","50","-qc","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-recombRegions")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Intra.ocular.pressure..Goldmann.correlated..mean/BOLTLMM.v3/Intra.ocular.pressure..Goldmann.correlated..mean-BOLT-LMM-v3.out" ,"genome", "plots" ,"-ymax", "50" ,"-qc", "-hits" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Intra.ocular.pressure..Goldmann.correlated..mean.RData" ,"-phenoFile", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt", "-title", "-recombRegions")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v6/Standing.height-BOLT-LMM-v6-chr2.out","all","plots","-ymax","50","-qc","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt","-title","-recombRegions","-bgen","/well/ukbiobank/imputation/final/chr2/out/chr2.test2.hrc.I4.bgen","-sample","/well/ukbiobank/phasing/final/phased_chunks/chr2.test2.sample","-dontComputeLD")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v11/Standing.height-BOLT-LMM-v11-23.out","23", "plots","-ymax","50","-qc","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-23.txt","-title","-recombRegions","-sex")
    
print(args)

type = "autosome"
if( "-sex" %in% args ) {
    type="sexchrom"
    # only function for the X chromosome at the moment!
    genotypeFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/byChrom/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095-sexchroms-%%"
}

BB.ps2snp = ukbiobank.ps2snp(type)
BL.ps2snp = ukbileve.ps2snp(type)

dataFile = args[1]
chroms = args[2]
# which chroms?
if( chroms %in% c("all","genome") ) {
    
    chroms = 1:22
    if( "-sex" %in% args ) chroms = 23:26
    
} else {
    
    chroms = parse.range.string(chroms)
    
}


plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

plotRegions=TRUE
computeLD=TRUE
allhits=FALSE # print region plots for all hits above genome-wide significance.
if("-dontPlotRegions"%in%args) plotRegions=FALSE
if("-dontComputeLD"%in%args) computeLD=FALSE
if("-allhits"%in%args) allhits=TRUE


# any specified bgen or genotypes file for computing LD?
bgenSampleFile=NULL
if("-bgen"%in%args) {
    genotypeFile=args[which(args=="-bgen")+1]
    if("-sample"%in%args) bgenSampleFile = args[which(args=="-sample")+1]
}

#########
#vers = "v3"
#########
#########
#vers = "v5"
#########
system( paste0("mkdir ",plotOutDir,"/Region.plots") )

if(computeLD){
    samplesFile = args[which(args=="-phenoFile")+1]
    system(paste0("cut -f1,2 -d' ' ",samplesFile," | tail -n +2 > samples.tmp"))
}

######## Get GWAS catalogue information
# NOTE: field descriptons are here: http://genome.ucsc.edu/cgi-bin/hgTables
catFile = NULL
if("-hits"%in%args) {
    # NOTE: this overrides the -hi for QC colours
    catInfo = args[which(args=="-hits")+1]
    load(catInfo,verbose=TRUE)
    
    catFile = catPheno  # just europeans
    colnames(catFile)[ncol(catFile)] = "Ancestry"
        
    catFile$Pvalue = catFile$V18
    catFile$SNP = catFile$V5
    catFile$BP = catFile$V4 # this is the chromEnd field
    catFile$CHR = gsub("chr","",catFile$V2)
    catFile$CHR[catFile$CHR%in%names(sexChroms)] = sexChroms[catFile$CHR[catFile$CHR%in%names(sexChroms)]]
    catFile$CHR = as.numeric(catFile$CHR)

    # adjust any Xchrom hits for PAR regions
    parHits = ( catFile$CHR == 23 ) & ( ( catFile$BP < par1 ) | ( catFile$BP > par2 ) )
    print( paste0( sum(parHits)," hits in PAR X") )
    catFile$CHR[parHits] = 25
    
    print( head(catFile) )
}


# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])


### read in GWAS results
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


#results = read.gwas(dataFile,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=50,QCexclusions=QCexclude,extraTitle=extraTitle) used pre: 6-07-2016
results = read.gwas(dataFile,chrom="genome",minmaf=minmaf,mininfo=mininfo,maxmiss = maxmiss,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle)

DF = results$DF


### read in recombination rates
sex = FALSE
if( "-sex" %in% args ) sex =TRUE
recombrates = read.recomb(sex)


# for each chromosome, get set of 'hits' based on recombination distances

if(!computeLD){
    hitsDataFile = paste0(results$Pvalset,"-hitsAndLDcalc.RData")
    if("-ldRData"%in%args) hitsDataFile = args[which(args=="-ldRData")+1]
    if(hitsDataFile%in%list.files()) load(hitsDataFile,verbose=TRUE)
}

if( !( "HITSchrom"%in%ls() )&(!"LDchrom"%in%ls() ) ){

    print("finding hits...")
    
    HITSchrom = vector("list", length(chroms) )
    LDchrom = vector("list", length(chroms) )

    for (chrom in chroms){

        print(chrom)
        
        Pvalset = paste0(gsub("genome",chrom,results$Pvalset))

#### Find hit regions
        recomb = NULL
        if( as.character(chrom) %in% names(recombrates) ) recomb = recombrates[[as.character(chrom)]] else print("no recombination map supplied. Using base-pair distance and genome-averate recombination rate instead.")
        
        HITS = find.hits.recomb(DF,chrom,minGap=0.125,buffer=25,recombRates=recomb)
        HITSchrom[[as.character(chrom)]] = HITS

#### compute LD with top SNPs (if they exist, and were asked for)
        if( length(HITS$topVariants) > 0 ){
            
            if(computeLD){
                if(grepl(".bgen",genotypeFile)){
                }
                ld = compute.ld.around.snp(HITS$topVariants,gsub("%%",chrom,genotypeFile),chrom,"samples.tmp",bgenSampleFile,DF)
                LDchrom[[as.character(chrom)]] = ld
            }
        } else {
            LDchrom[[as.character(chrom)]]=NULL
            print( paste0("No hits found for chromosome ",chrom) )
        }   
    }
}

################################### 
# plot hit regions over whole chromosomes (always do this - only if hits found!)

    
for (chrom in chroms){

    print(chrom)
    
    Pvalset = paste0(gsub("genome",chrom,results$Pvalset))
    HITS = HITSchrom[[as.character(chrom)]]
    
    if( length(HITS$topVariants) > 0 ){
        
        colors = rep("black",dim(HITS$DFthese)[1])
        for(i in 1:nrow(HITS$regs)){
            p = (HITS$DFthese$BP <= HITS$regs[i,2]) & (HITS$DFthese$BP >= HITS$regs[i,1])
            colors[p]=c("red","blue")[i%%2 + 1]
        }

#### plot the whole chromosome
        png(paste0("plots/Region.plots/",Pvalset,"-hitRegionsToPlot.png"),width=41,height=12,units="in",res=150)
        par(las=1,font.main=1,cex.axis=2,cex.lab=2,mar=c(7 ,7, 5 ,2))
        myManhattan(HITS$DFthese,ymax=Ymax,suggestiveline = FALSE,xpd=NA,cex=1,col=colors)
        points(HITS$DFthese$BP[HITS$DFthese$SNP%in%HITS$topVariants],-log10(HITS$DFthese$P2[HITS$DFthese$SNP%in%HITS$topVariants]),col="green",pch=17)
        dev.off()

    } else {
        print( paste0("No hits found for chromosome ",chrom) )
    }
}

    
################################### 

# plotting regions separately
#print(LDchrom)

if(plotRegions){
    
    for (chrom in chroms){

        print(chrom)
        
        Pvalset = paste0(gsub("genome",chrom,results$Pvalset))
        HITS = HITSchrom[[as.character(chrom)]]
        ld = LDchrom[[as.character(chrom)]]
        recomb = recombrates[[as.character(chrom)]]
        
        print( paste0("Processing ",length(HITS$topVariants)," top snps in chromosome ",chrom) )

#### plot the regions
        if(length(HITS$topVariants) > 0){
            sapply(1:length(HITS$topVariants),FUN=function(s) {
                regionToPlot = c(HITS$regs[s,1],HITS$regs[s,2])
                plot.gwas.region(DF,HITS$topVariants[s],catFile,ld,plotDir=paste0(plotOutDir,"/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=50,width=10^6,moreToPlot=TRUE,Pvalset=Pvalset)

                outWidth = diff(regionToPlot)
                abline(v=regionToPlot,col="darkgray",lty=3,lwd=3)
                text(regionToPlot[1],y=par("usr")[4],labels=paste0(round(outWidth/1000000,1),"MB region"),xpd=NA,cex=2,pos=3)
                dev.off()
            })
        }
    }
}


################################### REQUIRED FOR CONDITIONAL ANALYSIS
# output save
save(HITSchrom,LDchrom,file=paste0(results$Pvalset,"-hitsAndLDcalc.RData"))

# print set of SNPs for conditional analysis
allHitSNPs = unlist(sapply(HITSchrom,function(x) x$topVariants))

# convert RSIDs back to AffyIDs
allHitSNPs2 = sapply(allHitSNPs,function(s){
    if(s %in% BL.ps2snp$dbSNPRSID)  s = BL.ps2snp$AffySNPID[BL.ps2snp$dbSNPRSID==s]
    if(s %in% BB.ps2snp$dbSNPRSID)  s = BB.ps2snp$AffySNPID[BB.ps2snp$dbSNPRSID==s]
    return(s)
})
allHitSNPs2 = unlist(allHitSNPs2)

write.table(allHitSNPs2,file=paste0(results$Pvalset,"-topVariants.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)

quit()  # finish here if running from the command line



# plot the cluster plots for rare SNPs
outdir = paste0(getwd(),"/clusterplots_rare",rare,"_hits")
snps = paste(rareHitSnps,collapse=",")
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir)
system(forSystem)

# combine with region plots
sapply(rareHitSnps,combine.cluster.regions.plots,prefix="recombRegions")






# make cluster plot for some chosen snps.
########################
# CLUSTER PLOTS FOR ALL TOP HITS WITH NO KNOWN HITS NEARBY

novelSnps = c()
for (chrom in 1:22){
    print(chrom)
    Pvalset = paste0(gsub("genome",chrom,results$Pvalset))
    HITS = HITSchrom[[chrom]]
    isKnown = apply(HITS$regs,1,function(region){
        known = catFile[(catFile$BP >= region[1]) & (catFile$BP <= region[2]) & (catFile$CHR == chrom),]
        return(nrow(known))
    })
    snps = HITS$topVariants[isKnown == 0]  # no known hits
    print(length(snps))
    novelSnps = c(novelSnps,snps)
}

snps = paste(novelSnps, collapse=",")
outdir = paste0(getwd(),"/clusterplots_novel_hits_allbatches")
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir)
system(forSystem)

## combine region and cluster plots
chrom=2
toPlot = novelSnps[novelSnps%in%DF$SNP[DF$CHR!=chrom]]
sapply(toPlot,combine.cluster.regions.plots,prefix="recombRegions")

sapply("Affx-18444914",combine.cluster.regions.plots,prefix="recombRegions")




########################
# CLUSTER PLOTS FOR ODD HITS (identified by visually inspecting region plots)

snpToPlot = read.delim("RegionsToInvestigate.txt",header=TRUE,sep="\t")

sum(snpToPlot$affyid%in%QCexclusions)
snps = paste(unique(snpToPlot$affyid[snpToPlot$type=="lw"]),collapse=",")

outdir = paste0(getwd(),"/clusterplots_singleton_hits")    
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir,' -b Batch_b001-Batch_b095,UKBiLEVEAX_b1-UKBiLEVEAX_b11 -orig')

system(forSystem)

# all batches
outdir = paste0(getwd(),"/clusterplots_singleton_hits_allbatches")    
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir)

system(forSystem)

# non lone wolf
snps = paste(unique(snpToPlot$affyid[snpToPlot$type!="lw"]),collapse=",")
outdir = paste0(getwd(),"/clusterplots_other_odd_hits_allbatches")    
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir)

system(forSystem)




########################
# Aggregate measures (do this on all SNPs)

hits = DF$P < 5e-8

png(paste0("plots/",Pvalset,"-histCompare-missing.png"),width=1000,height=1000,res=150)
hist(DF$F_MISS[!hits],breaks=seq(0,0.05,0.001), col=add.alpha("blue",0.8),xlim=c(0,0.05))
hist(DF$F_MISS[hits],breaks=seq(0,0.05,0.001), col=add.alpha("red",0.8),add=TRUE)
dev.off()

Hits = rep("> 5e-8",nrow(DF)); Hits[hits] = "< 5e-8"
png(paste0("plots/",Pvalset,"-boxplotCompare-missing.png"),width=1000,height=1000,res=150)
boxplot(log10(DF$F_MISS)~Hits,col=c("red","blue"),ylab="log10(missing)")
dev.off()
png(paste0("plots/",Pvalset,"-boxplotCompare-maf.png"),width=1000,height=1000,res=150)
boxplot(log10(DF$MAF)~Hits,col=c("red","blue"),ylab="log10(MAF)")
dev.off()

lmMiss = lm(DF$F_MISS~Hits)
lmPmaf = glm(factor(Hits)~log10(DF$F_MISS) + DF$MAF,family="binomial")


# LREG vs LMM
x = DF$P2
y = DF$P_LINREG
y[-log10(y)>Ymax] = 10^-Ymax
png(paste0("plots/",Pvalset,"-lregLMMCompare-pvals.png"),width=1000,height=1000,res=150)
plot(-log10(x),-log10(y),ylab="Linear regression",xlab="LMM",xlim=c(0,50),ylim=c(0,50))
abline(v=-log10(5e-8),h=-log10(5e-8),col="red",lty=3)
abline(0,1,col="red")
dev.off()






##############################
# Specific regions a bit wider -- comparing Height to POB

######
# read in POB GWAS data
DFheight = DF
dataFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v3/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v3.out"
DFpob = DF
DFpob$Pvalue=DFpob$P
DF = DFheight

### set wider width for plotting
width = 2*10^6
##########
chrom=4 # around "Affx-24457234" on Chromosome 4
topSNP ="Affx-24457234"
region = c(DF$BP[DF$SNP==topSNP] - width/2,DF$BP[DF$SNP==topSNP] + width/2)
these = (DF$BP >= region[1]) & (DF$BP <= region[2]) & (DF$CHR == chrom) & (DF$P < 10^-8)
spVariants = DF$SNP[these]
##########
chrom=2
topSNP= "Affx-36377258"
topSNP= "Affx-18554975"
region = c(DF$BP[DF$SNP==topSNP] - width/2,DF$BP[DF$SNP==topSNP] + width/2)
these = (DF$BP >= region[1]) & (DF$BP <= region[2]) & (DF$CHR == chrom) & (DF$P < 10^-8)
spVariants = DF$SNP[these]
##########


Pvalset = paste(basename(dataFile),".chr",chrom,".maf",minmaf,".miss",maxmiss,".pruned",extraTitle,sep="")
ld = compute.ld.around.snp(spVariants,genotypeFile,chrom,"samples.tmp")
sapply(spVariants,FUN=function(s) plot.gwas.region(DF,s,catFile,ld,plotDir=paste0(plotOutDir,"/Region.plots"),region=region) )


# plot the pob SNPs with height SNPs in this region. Are they mirror image??
t = (DF$BP >= region[1]) & (DF$BP <= region[2]) & (DF$CHR == chrom)
t2 = (DFpob$BP >= region[1]) & (DFpob$BP <= region[2]) & (DFpob$CHR == chrom)

ylimit = round(max(c(-log10(DF$P[t]),-log10(DFpob$P[t2]))))

plot.gwas.region(DF,topSNP,catFile,ld,plotDir=paste0(plotOutDir,"/Region.plots"),region=region,extra="withPOBresults",moreToPlot=TRUE,ylim=c(-ylimit,ylimit))
points(DFpob$BP[t2],log10(DFpob$P[t2]),col="blue",pch=16,cex=2)
abline(h=log10((5e-8)),col="red")
abline(h=0,col="gray")
legend("bottomleft",legend="Place of birth",cex=2,bty="n")
legend("topleft",legend="Standing height",cex=2,bty="n")
dev.off()



########################
# Compare GWAS versions

v1="v3"
v2="v4"

dataFile2 = gsub(vers,v2,dataFile)

resultsV = read.gwas(dataFile2,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=50,QCexclusions=QCexclude)

DFV = resultsV$DF
x = DF$P
y = DFV$P
sum(DFV$SNP!=DF$SNP)
l = gsub(v2,paste0(v1,"-and-",v2),resultsV$Pvalset)
colors = rep(NA,dim(DF)[1])
colors[DF$MAF < 0.001] = "red"
colors[(DF$MAF >= 0.001) & (DF$MAF < 0.01)] = "purple"
colors[(DF$MAF >= 0.01) & (DF$MAF < 0.1)] = "blue"
colors[(DF$MAF >= 0.1) & (DF$MAF < 0.2)] = "green"
colors[(DF$MAF >= 0.2)] = "yellow"
Order = order.by.number.occurrences(colors)
diffs = -log10(x) + log10(y)
    
png(paste0("plots/",l,"-pvaluesCompare.png"),width=1000,height=1000,res=150)
plot(x[Order],y[Order],xlab="Original covariates",ylab="Plus dQC and CR as covariates",col=colors[Order],pch=16)
dev.off()

png(paste0("plots/",l,"-pvaluesLogCompare.png"),width=1000,height=1000,res=150)
plot(-log10(DF$P2[Order]),-log10(DFV$P2[Order]),xlab="Original covariates",ylab="Plus dQC and CR as covariates",col=colors[Order],pch=16)
dev.off()

png(paste0("plots/",l,"-pvaluesLogDiffCompare.png"),width=1000,height=1000,res=150)
plot(-log10(DF$P[Order]),-log10(x[Order]) + log10(y[Order]),xlab="Pvalue Original",ylab="( Original covariates ) - ( Plus dQC and CR as covariates )",col=colors[Order],pch=16)
dev.off()

png(paste0("plots/",l,"-pvaluesLogDiffMissCompare.png"),width=1000,height=1000,res=150)
plot(log10(DF$F_MISS[Order]),-log10(x[Order]) + log10(y[Order]),xlab="log10 Fraction missing",ylab="( Original covariates ) - ( Plus dQC and CR as covariates )",col=colors[Order],pch=16)
dev.off()

png(paste0("plots/",l,"-pvaluesHistCompare.png"),width=1000,height=1000,res=150)
hist(-log10(x[Order]) + log10(y[Order]),xlab="( Original covariates ) - ( Plus dQC and CR as covariates )",breaks=100,ylim=c(0,10^4))
dev.off()

#  how many SNPs went up or below genome-wide significance?
sum((DF$P < 5e-8) & (DFV$P >= 5e-8))  # 10 moved above
sum((DF$P >= 5e-8) & (DFV$P < 5e-8))  # 6 moved below

(x/y)[(DF$P < 5e-8) & (DFV$P >= 5e-8)]
(x/y)[(DF$P >= 5e-8) & (DFV$P < 5e-8)]


#  what are the properties of the SNPs that changed?
bigDiff = abs(-log10(DF$P) + log10(DFV$P)) > 0.1; bigDiff[is.na(bigDiff)] = FALSE
snps = paste(DF$SNP[bigDiff],collapse=",")
outdir = paste0(getwd(),"/clusterplots_v3_v4_differences")    
forSystem = paste0( 'Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ',snps,' ',outdir,' -b Batch_b001,Batch_b064,Batch_b070,Batch_b095,UKBiLEVEAX_b1,UKBiLEVEAX_b5')

system(forSystem)

