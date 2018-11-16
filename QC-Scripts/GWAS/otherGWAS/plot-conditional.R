source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

args = commandArgs(TRUE)

args = c("Standing.height-BOLT-LMM-v3.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-snpRegionList.txt","-readRaw","-hits","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData","-chr","all","-qc","-ymax","50","-resultsPrefix","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v3/Standing.height-BOLT-LMM-v3-cond","-phenoFile","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt")

regionList = args[1]
ldStats = gsub("snpRegionList.txt","hitsAndLDcalc.RData",regionList)
chroms = args[which(args=="-chr")+1]
resultsPrefix = args[which(args=="-resultsPrefix")+1]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""

# do we fix the y-axis?
Ymax = FALSE
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])


### read in recombination rates
recombrates = list()
for(i in c(1:22,"X_nonPAR","X_PAR1","X_PAR2")){
    print(i)
    r = read.table(paste0("/well/donnelly/spain_structure/phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr",i,"_combined_b37.txt"),header=TRUE)
    recombrates[[i]] = r
}

# read QC and snp annotations
print("Reading QC snps lists...")
QCSNPList = read.SNPQC.files(justSNPs=TRUE)
QCexclude = unique(c(QCSNPList$arraySNPs,QCSNPList$concordanceSNPs))                             # NOTE: anything created after Friday June 24th uses this QCexclude. Anything created before that also excluded the 'image' SNPs    . QCexclude only required if using GWAS versions < v5
BB.ps2snp = ukbiobank.ps2snp("autosome")
BL.ps2snp = ukbileve.ps2snp("autosome")

genotypeFile = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-oxfordqc"

samplesFile = args[which(args=="-phenoFile")+1]
system(paste0("cut -f1,2 -d' ' ",samplesFile," | tail -n +2 > samples.tmp"))


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
    print( head(catFile) )
}



# load output from region.plots.R (this is based on the non-conditional results)
load(file=ldStats,verbose=TRUE)


# load list of regions
regions = read.table(regionList,header=TRUE)


# which chroms?
if(chroms!="genome") {
    if(chroms=="all") chroms = 1:22 else chroms = parse.range.string(chroms)
}

# read in unconditional gwas
uncondGWAS = gsub("-cond",".out",resultsPrefix)
DF = read.gwas(uncondGWAS,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle)


# do we read raw data from conditional gwas??
if("-readRaw" %in% args) readRawFiles=TRUE else readRawFiles=FALSE


# function to read in region files
get.region.data <- function(reg,regs,readRawFiles=FALSE){
    
    print(reg)
    chrom = unique(regs$chrom)
    index = regs[reg,"region"]
    start = regs[reg,"start"]
    end = regs[reg,"end"]
    dataFile = paste0(resultsPrefix,"-chr",chrom,"-reg",index,"-",start,"-",end,".out")
    print(dataFile)
                                        # do we need to read in the raw data?
    filteredDataFile = gsub(paste0(end,".out"),paste0(end,"-sub.out"),dataFile)
    print(filteredDataFile)
    if( ( !basename(filteredDataFile) %in% list.files(dirname(filteredDataFile)) ) | readRawFiles ){
        if( !basename(dataFile) %in% list.files(dirname(dataFile)) ) {
            return(NULL)
        }
        print('reading raw data')
        df = read.gwas(dataFile,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle)        
                                        # first just filter the output files and save with list
        if(is.null(df)) return(NULL)
       
        filteredDF =  filter(df$DFraw,(BP >= start)&(BP <= end)&(CHR == chrom))        
        write.table(filteredDF,file=filteredDataFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    }
    
    print('reading small data')
    filteredDF = read.gwas(filteredDataFile,chrom="genome",minmaf=0.001,mininfo=0.3,maxmiss = 0.05,Ymax=Ymax,QCexclusions=QCexclude,extraTitle=extraTitle)        
    
    return(filteredDF)
}


# combine all filtered data 
fChrom = sapply(chroms,function(chrom){
    regs = regions[regions$chrom==chrom,]    
    f = sapply(1:nrow(regs), function(reg){
        filteredDF = get.region.data(reg,regs,readRawFiles=readRawFiles)
    },simplify=FALSE)
})


# get ld stats and hit regions
HITSCondChrom = list()
LDCondChrom = list()
for(chrom in chroms){
    print(chrom)    
    recomb = recombrates[[as.character(chrom)]]
    regs = regions[regions$chrom==chrom,]

    if(nrow(regs)==0) {
        print("No regions in this chromosome!")
        HITSCondChrom[[chrom]]=NULL
        LDCondChrom[[chrom]]=NULL
        next
    }

    print(paste0(nrow(regs)," regions in this chromosome."))
                                        # combine results for all regions at once
    thisDF = Reduce(rbind,sapply(fChrom[[chrom]],function(x) x$DF,simplify=FALSE))
    
                                        # get ld stats for any hits in conditional regions
    HITS = find.hits.recomb(thisDF,chrom,minGap=0.125,buffer=25,recombRates=recomb)
    HITSCondChrom[[chrom]] = HITS
    ld = compute.ld.around.snp(HITS$topVariants,genotypeFile,chrom,"samples.tmp")
    LDCondChrom[[chrom]] = ld
}

######## RE-RUN LD and HITS bit for region chr1:66 and chr2:1 ==> re-run plots for these only

# save original ldstats along with conditional ld stats
save(HITSchrom,LDchrom,HITSCondChrom,LDCondChrom,LDrare,rareHitSnps,file=gsub("-hitsAndLDcalc.RData","-withConditional-hitsAndLDcalc.RData",ldStats))


plot.conditional <- function(chrom,theseRegs=NULL){
    outFile = gsub("-hitsAndLDcalc.RData","",ldStats)
    outFile = gsub("genome",chrom,outFile)
    
    print(chrom)
    
    HITS = HITSchrom[[chrom]]
    HITSCond = HITSCondChrom[[chrom]]    
    regs = regions[regions$chrom==chrom,]    
    recomb = recombrates[[as.character(chrom)]]
    ld = LDchrom[[chrom]]
    ldCond = LDCondChrom[[chrom]]
    
    if(is.null(HITSCond)) return(NULL)

#### plot the whole chromosome
    colors = rep("black",dim(HITSCond$DFthese)[1])
    for(i in 1:nrow(HITSCond$regs)){
        p = (HITSCond$DFthese$BP <= HITSCond$regs[i,2]) & (HITSCond$DFthese$BP >= HITSCond$regs[i,1])
        colors[p]=c("red","blue")[i%%2 + 1]
    }

    png(paste0("plots/Region.plots/",outFile,"-cond-hitRegionsToPlot.png"),width=41,height=12,units="in",res=150)
    par(las=1,font.main=1,cex.axis=2,cex.lab=2,mar=c(7 ,7, 5 ,2))
    # match xlimits for the original hit plots
    xlimits = c(floor(max(HITS$DFthese$BP) * -0.03),ceiling(max(HITS$DFthese$BP) * 1.03))
    myManhattan(HITSCond$DFthese,ymax=Ymax,xlimits=xlimits,suggestiveline = FALSE,xpd=NA,cex=1,col=colors)
    points(HITSCond$DFthese$BP[HITSCond$DFthese$SNP%in%HITSCond$topVariants],-log10(HITSCond$DFthese$P2[HITSCond$DFthese$SNP%in%HITSCond$topVariants]),col="green",pch=17)
    dev.off()
    
    
#### plot each of the original regions separately - before and after conditioning.
    if(is.null(theseRegs)) theseRegs = 1:nrow(regs)
    
    sapply(theseRegs,FUN=function(s) {
        start = regs[s,"start"]
        end = regs[s,"end"]
        regionToPlot = c(start,end) + c(-10^4,10^4)
        print(regionToPlot)

        hitSNPsOrig = HITS$topVariants[(HITS$regs$start <= end) & (HITS$regs$start >= start )]
        hitSNPsCond = HITSCond$topVariants[(HITSCond$regs$start <= end) & (HITSCond$regs$start >= start )]
        PvalsetOrig = paste0(outFile,"-reg",start,"-",end,"-orig")
        PvalsetCond = paste0(outFile,"-reg",start,"-",end,"-cond")
                                        # only plot if there are still conditional hits
        if( length(hitSNPsCond) == 0 ){
            
            print("not conditional in this region")
            return(NULL)
            
        } else {
            print("plotting this region")
                                        # original hits
            sapply(1:length(hitSNPsOrig),function(t){
                
                plot.gwas.region(DF$DF,hitSNPsOrig[t],catFile,ld=ld,plotDir=paste0("plots/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=Ymax,width=10^3,moreToPlot=TRUE,Pvalset=PvalsetOrig)
                w = c(HITS$regs[HITS$topVariants==hitSNPsOrig[t],1],HITS$regs[HITS$topVariants==hitSNPsOrig[t],2])
                outWidth = diff(w)
                abline(v=w,col="darkgray",lty=3,lwd=3)
                text(w[1],y=par("usr")[4],labels=paste0(round(outWidth/1000000,1),"MB region"),xpd=NA,cex=2,pos=3)
                dev.off()
            })
                                        # original hits
            sapply(1:length(hitSNPsCond),function(t){
                
                thisDF = Reduce(rbind,sapply(fChrom[[chrom]],function(x) x$DF,simplify=FALSE))
                
                plot.gwas.region(thisDF,hitSNPsCond[t],catFile,ld=ldCond,plotDir=paste0("plots/Region.plots"),recombRates=recomb,region=regionToPlot,Ymax=Ymax,width=10^3,moreToPlot=TRUE,Pvalset=PvalsetCond)
                w = c(HITSCond$regs[HITSCond$topVariants==hitSNPsCond[t],1],HITSCond$regs[HITSCond$topVariants==hitSNPsCond[t],2])
                outWidth = diff(w)
                abline(v=w,col="darkgray",lty=3,lwd=3)
                text(w[1],y=par("usr")[4],labels=paste0(round(outWidth/1000000,1),"MB region"),xpd=NA,cex=2,pos=3)
                dev.off()
            })
            
            # combine the plots
            print('printing plots combined in a pdf file...')
            plotNames = paste0(outFile,"-reg",start,"-",end)
            plots = grep(plotNames,list.files(path="plots/Region.plots/"),value=TRUE)
            pdf(paste0("plots/Region.plots/",plotNames,"-combinedPlots.pdf"),width=41,height=12)
            par(mar=c(0,0,0,0))
        
            for( i in plots[grep("orig",plots)] ){
                plot.new()
                img <- readPNG(paste0("plots/Region.plots/",i))
                rasterImage(img,0,0,1,1)
            }
            for( i in plots[grep("cond",plots)] ){
                plot.new()                                
                img <- readPNG(paste0("plots/Region.plots/",i))
                rasterImage(img,0,0,1,1)
            }            
            dev.off()
            
        }
    })                  
}


sapply(chroms,function(chrom) plot.conditional(chrom) )
