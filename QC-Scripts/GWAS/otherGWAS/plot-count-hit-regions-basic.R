#####################
# Script to plot the output of count-hit-regions-basic.R
#####################

#####################
# Preliminaries
library(VennDiagram)
library(venneuler) # from http://www.rforge.net/venneuler/files
library(plyr)

args = commandArgs(TRUE)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

sourceNames = c("IMP"="UK Biobank \n(imputed)","GENO"="UK Biobank \n(genotyped)","GIANT"="GIANT (2014)")

# args = c("-counts","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/Standing.height-BOLT-LMM-v19-chr%%.out.chr%%.maf0.info0.pruned.1MbWin.IMP-5e-08.txt,/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/Standing.height-BOLT-LMM-v16.out.chr%%.maf0.miss1.pruned.1MbWin.GENO-5e-08.txt,/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/GIANT.chr%%.1MbWin.GIANT-5e-08.txt","-title","hit-counts-Imp.v16.Geno.v19.GIANT-1MbWin","-chr","1-2","-plotdir","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height","-minmaf","0","-mininfo","0",)

print(args)

title = args[which(args=="-title")+1]
dataFiles = str_split(args[which(args=="-counts")+1],",")[[1]]

if(length(dataFiles) != 3) {
    print("NEED 3 input files! Stopping.")
    quit()
}

plotdir=args[which(args=="-plotdir")+1]
outDir = args[which(args=="-outdir")+1]
print(outDir)

chroms = parse.range.string(args[which(args=="-chr")+1])

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

title = paste0(title,".maf",minmaf,".info",mininfo,".miss",maxmiss)


######
# Read in the counts

myData = lapply(dataFiles,function(file) {

    if(grepl("v19|v20|IMP",file)) type="IMP"
    if(grepl("v16|GENO",file)) type="GENO"
    if(grepl("GIANT",file)) type="GIANT"
    print(file)
    print(type)
    
    theseCounts = lapply(chroms,function(chr){
        
        f = gsub("%%",chr,file)
        
        if(file.exists(f)){
            counts = read.table(f,header=TRUE)
        } else {
            print(paste0("Missing file for this chromosome: ",f))
            return(NULL)
        }
        return(counts)
    })
    
    Counts = abind(theseCounts,along=1)
    # Apply any filters here
    toExclude = (Counts[,"max_INFO"] <= mininfo)|(Counts[,"max_MAF"] <= minmaf)

    print("Not counting hits in this many regions because of filters.")
    print(sum(toExclude,na.rm=TRUE))
    
    Counts[which(toExclude),"counts"] = 0

    colnames(Counts)[!colnames(Counts)%in%c('chr','start','end')] = paste0(colnames(Counts)[!colnames(Counts)%in%c('chr','start','end')],".",type)
    
    Counts = dplyr::tbl_df(data.frame(Counts))
    
    return(Counts)
})


######
# Align the counts by region positions

mergedCounts = plyr::join_all(myData, by=c('chr','start','end'), type='inner')

nRegions = nrow(mergedCounts)

chroms = unique(mergedCounts[["chr"]])
print(paste0("Plotting ",nRegions," regions across chromosomes ",paste0(chroms,collapse=",")))


############
# Save this info

save(mergedCounts,chroms,file=paste0(outDir,"/",title,".RData"))


############
# Count the overlaps

inIMP = mergedCounts[,"counts.IMP"]>0
inGENO = mergedCounts[,"counts.GENO"]>0
inGIANT = mergedCounts[,"counts.GIANT"]>0

imp = sum(inIMP)
geno = sum(inGENO)
giant = sum(inGIANT)

n12 = sum(inIMP&inGENO)
n13 = sum(inIMP&inGIANT)
n23 = sum(inGENO&inGIANT)
n12only = sum(inIMP&inGENO&!inGIANT)
n13only = sum(inIMP&!inGENO&inGIANT)
n23only = sum(!inIMP&inGENO&inGIANT)
n123 = sum(inIMP&inGENO&inGIANT)
n1 = sum(inIMP&!inGENO&!inGIANT)
n2 = sum(!inIMP&inGENO&!inGIANT)
n3 = sum(!inIMP&!inGENO&inGIANT)
nUnion = sum(inIMP|inGENO|inGIANT)
n1all = sum(inIMP)
n2all = sum(inGENO)
n3all = sum(inGIANT)
    
imp = which(inIMP)
geno = which(inGENO)
giant = which(inGIANT)


############
#  Draw a venn diagram of the overlaps

studyCols = c(impColour,genoColour,giantColour)
names(studyCols) = c("IMP","GENO","GIANT")

pdf(paste0(plotdir,"/plots/",title,".pdf"),bg="transparent")

#  draw.triple.venn(n1all, n2all, n3all, n12, n23, n13, n123,category=sourceNames[c("IMP","GENO","GIANT")])

m = data.frame(elements=c(imp,geno,giant), sets=c(rep("IMP",length(imp)),rep("GENO",length(geno)),rep("GIANT",length(giant))))

vd = venneuler(m)


#plot.VennDiagram(vd,main=title)

plot(vd,main=paste0(title,"; ",nRegions," regions"),col.txt=NA,col=studyCols[rownames(vd$centers)])
                                        #text(vd$centers, vd$labels, col=1)
#labelAdj = cbind(c(-0.06,-0.17,0),c(0.25,-0.06,0))
labelAdj = cbind(c(0,0,0),c(0,0,0))
text(vd$centers+labelAdj, sourceNames[vd$labels], col=studyCols[rownames(vd$centers)],font=2)

#text(vd$centers+labelAdj, sourceNames[vd$labels], col=studyCols[rownames(vd$centers)],font=2,pos=2)

plot(0,1,axes=FALSE,col=NA,xlab=NA,ylab=NA)
getPerc <- function(x,total=nUnion) paste(x," (",round(100*x/total,3),"%)")

legend("top",bty="n",legend=c(paste0("All IMP:   ",getPerc(n1all)),
                 paste0("only IMP:   ",getPerc(n1)),
                 paste0("All GENO:   ",getPerc(n2all)),
                 paste0("only GENO:   ",getPerc(n2)),
                 paste0("All GIANT:   ",getPerc(n3all)),
                 paste0("only GIANT:   ",getPerc(n3)),                 
                 paste0("IMP & GENO:   ",getPerc(n12)),
                 paste0("only IMP & GENO:   ",getPerc(n12only)),
                 paste0("IMP & GIANT:   ",getPerc(n13)),
                 paste0("only IMP & GIANT:   ",getPerc(n13only)),
                 paste0("GENO & GIANT:   ",getPerc(n23)),
                 paste0("only GENO & GIANT:   ",getPerc(n23only)),
                 paste0("IMP & GENO & GIANT:   ",getPerc(n123)),
                 paste0("Union:   ",getPerc(nUnion)) ) )

dev.off()


print("DONE!")
print(paste0(plotdir,"/plots/",title,".pdf"))
