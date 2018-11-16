## Script to filter output from King;
## It provides a list of samples excluded from kinship table on grounds of QC
## It also outputs a filtered kinship table, with kinship classes attached
## It also outputs an Rdata file with family configurations (after filtering)

args = commandArgs(trailingOnly=TRUE)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0","-ibs0","0.0012","-outliers","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt")

print(args)
h = args[-c(which(args%in%c("-in","-out","-ibs0","-outliers")),1+which(args%in%c("-in","-out","-ibs0","-outliers")))]
for(helperScript in h){
    source(helperScript)
}

library(igraph)
setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness"))

KingOutFile = args[which(args=="-in")+1]
OutputFile = gsub(".kin0","",basename(KingOutFile))
OutDir = args[which(args=="-out")+1]
OutlierFile = args[which(args=="-outliers")+1]
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])

## FUNCTIONS
otherInfo = read.multiple.batch.info(c("Ethnic.background","Place.of.birth.in.UK...east.co.ordinate","Place.of.birth.in.UK...north.co.ordinate","Pops","Chars","Colors","Batch"))

## get kinship data
kin = read.table(KingOutFile,header=T,stringsAsFactors=F)

## get classes
kin$class = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)
table(kin$class)
kin$col = colDef[kin$class]

## get het-missing data
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"),verbose=TRUE)

## get het.miss outliers
exclusions = read.table(OutlierFile,stringsAsFactors=FALSE)[,1]

# Exclude het-missing outliers
hetmissout = (kin$ID1%in%exclusions)|(kin$ID2%in%exclusions )
kinHetmissout=kin[!hetmissout,]

write.table(kinHetmissout[,c("FID1","ID1","FID2","ID2","N_SNP","HetHet","IBS0","Kinship")],file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-hetmiss-out.kin0"),quote=FALSE,row.names=FALSE,col.names=TRUE )

quit()
# First exclude those with more than 200 relatives, and other poor quality samples
##########
maxNrelatives = 200
##########

# first exclude het-missing outliers, then calculate number of 3rd-degree relatives.
kin = kin[!hetmissout,]

totalKin = table(c(kin$ID1[kin$class=="3rd degree"],kin$ID2[kin$class=="3rd degree"]))

excIndex = (totalKin >= maxNrelatives) | ( names(totalKin)%in%exclusions)
#excIndex = ( names(totalKin)%in%exclusions)
tooManyRelatives = names(totalKin)[totalKin >= maxNrelatives]

print( paste0(sum(totalKin >= maxNrelatives)," samples with >= ",maxNrelatives," relatives.") )
print( paste0(sum(tooManyRelatives%in%exclusions)," of these are already in het/missing exclusion list.") )


exc = (kin$ID1%in%(names(totalKin)[excIndex]) )|(kin$ID2%in%(names(totalKin)[excIndex]) )
kin2=kin[!exc,]   # <===== This is where the filtering for the pruned-200.kin0 file happens.

# list of exclusions, including hetmiss and over 200 relatives after pruning for het miss
allExcludedNames = unique( c(exclusions,tooManyRelatives) )
print( paste0( length(allExcludedNames) ," samples in het-miss outliers and/or >=",maxNrelatives," relatives after pruning het-miss outliers") )

# check properties of other excess related samples
third = kin2$class=="3rd degree"
totalKin2 = table(c(kin2$ID1[third],kin2$ID2[third]))
totalKin10 = totalKin2[totalKin2 > 10 ]

print( paste0( length(totalKin10)," samples left with more than 10 3rd-degree relatives." ) )

# create graph for kinship prior to filtering others out
kin2 = kin2[,!colnames(kin2)%in%c("FID1","FID2")]

print("Checking family formations for dubious pairs.")

                                        # create graph - including 3rd degree
network <- graph.data.frame(kin2, directed=F)
V(network)$size <- log(degree(network))+1
E(network)$color <- E(network)$col
fams = clusters(network)

                                        # create graph - excluding 3rd degree
kin2No3 = kin2[kin2$class!="3rd degree",]

networkNo3 <- graph.data.frame(kin2No3, directed=F)
V(networkNo3)$size <- log(degree(networkNo3))+1
E(networkNo3)$color <- E(networkNo3)$col
famsNo3 = clusters(networkNo3)

                                        # create graph - only 3rd degree
kin2Only3 = kin2[kin2$class=="3rd degree",]

networkOnly3 <- graph.data.frame(kin2Only3, directed=F)
V(networkOnly3)$size <- log(degree(networkOnly3))+1
E(networkOnly3)$color <- E(networkOnly3)$col
famsOnly3 = clusters(networkOnly3)


## plot graphs of large n relatives
theseNames = names(totalKin2)[totalKin2 > 10]

plotting=TRUE
if( plotting ) {
    

    inhere = (kin2$ID1%in%theseNames)|(kin2$ID2%in%theseNames)

    rels = unique(c(kin2$ID1[inhere],kin2$ID2[inhere]))
    theseKin = (kin2$ID1%in%rels)&(kin2$ID2%in%rels)
    bsk.network <- graph.data.frame(kin2[theseKin,], directed=F)
    V(bsk.network)$size <- log(degree(bsk.network))+1
    V(bsk.network)$color <- ifelse(V(bsk.network)$name%in%theseNames,"black","gray")
    E(bsk.network)$color <- E(bsk.network)$col
    h =clusters(bsk.network)
    subgraph =  igraph::induced_subgraph(bsk.network,theseNames,impl="copy_and_delete")


                                        # only consider 3rd degree relatives
    rels3 = V(bsk.network)$name[V(bsk.network)$name%in%V(networkOnly3)$name]
    subgraph2 = igraph::induced_subgraph(networkOnly3,rels3,impl="copy_and_delete")
    V(subgraph2)$color <- ifelse(V(subgraph2)$name%in%theseNames,"black","gray")
    subFam = clusters(subgraph2)


    png(paste0("plots/",OutputFile,"-network-graph-10To200Relatives-beforeExclusions-graph.png"),width=1000,height=1000,res=150)
    plot.igraph(bsk.network,main=paste0("Samples with 10-200 3rd-degree 'relatives',\nand their relatives, before pruning"),vertex.frame.color=NA,vertex.label=NA,vertex.size=1.1)
    dev.off()

    png(paste0("plots/",OutputFile,"-network-graph-10To200Relatives-beforeExclusions-graph2.png"),width=1000,height=1000,res=150)
    plot.igraph(subgraph,main=paste0("Samples with 10-200 3rd-degree 'relatives', before pruning"),vertex.frame.color=NA,vertex.label=NA,vertex.size=1.1)
    dev.off()

    png(paste0("plots/",OutputFile,"-network-graph-10To200Relatives-beforeExclusions-graph3.png"),width=1000,height=1000,res=150)
    plot.igraph(subgraph2,main=paste0("Samples with 10-200 3rd-degree 'relatives', and ther relatives\n3rd-degree only, before pruning"),vertex.frame.color=NA,vertex.label=NA,vertex.size=1.1)
    dev.off()

                                        # write a list of samples that have dubious relatedness metrics
    highrels = c(names(totalKin)[totalKin >= maxNrelatives],theseNames)  # 229 (105 are in the het-miss outliers anyway)
    print( paste0(length(highrels), " samples with > 200 relatives, or > 10 after pruning") )
    highrels = highrels[!highrels%in%exclusions] # 124
    print( paste0(length(highrels), " samples with the above property which aren't also in the het-missing outlier list.") )

    write.table(highrels,file=paste0(OutputFile,"-excess-relatives.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE )


                                        # plot these family units (third degree only) on a map

    these = otherInfo$PIID%in%V(subgraph2)$name
    x = otherInfo$Place.of.birth.in.UK...east.co.ordinate[these]
    y = otherInfo$Place.of.birth.in.UK...north.co.ordinate[these]
    rand = sample(1:nrow(otherInfo),10000)
    x[is.na(x)] = -1
    y[is.na(y)] = -1

    colors = subFam$membership[otherInfo$PIID[these]]
    shapes = rep(1,length(x)); shapes[names(colors)%in%theseNames] = 16

    png(paste0("plots/",OutputFile,"-network-graph-10To200Relatives-beforeExclusions-graph3-POB-UK.png"),width=1500,height=2000,res=200)
    plot(otherInfo$Place.of.birth.in.UK...east.co.ordinate[rand],otherInfo$Place.of.birth.in.UK...north.co.ordinate[rand],col="gray")
    points(x,y,col=colors,pch=shapes,cex=0.5)
    dev.off()

                                        # plot these on het-missing
    y = Table$het.corrected
    x = Table$logit.miss
    colors = rep("gray",length(y))
    colors[Table$IID%in%exclusions] = "red"
    colors[Table$IID%in%names(totalKin)[totalKin >= maxNrelatives]] = "blue"
    colors[Table$IID%in%theseNames] = "green"
    Order = order.by.number.occurrences(colors)

    png(paste0("plots/",OutputFile,"het-missing-beforeExclusions.png"),width=1500,height=2000,res=200)
    plot(x[Order],y[Order],col=colors[Order],xlab="logit(missing)",ylab="PC-corrected heterozygosity")
    legend("topright",legend=c("het-missing outliers","more than 200 3rd-degree relatives","more than 10 3rd-degree relatives after excluding the above samples."),col=c("red","blue","green"),pch=1)
    dev.off()

                                        # are they in the same batches??
    bat = factor(otherInfo$Batch)
    png(paste0("plots/",OutputFile,"-gt10-relatives-by-batch-beforeExclusions.png"),width=2000,height=1000,res=150)
    barplot(table(bat[otherInfo$PIID%in%theseNames]),las=2)
    dev.off()

    chisq.test(table(bat[otherInfo$PIID%in%theseNames]),table(bat[!otherInfo$PIID%in%theseNames]))

}


# Now write a new kinship table after excluding the dubious 3rd-degree pairs, as well as the het/missing exclusions.
kin2$FID1 = kin2$ID1
kin2$FID2 = kin2$ID2

kinFiltered = kin2
kinFiltered = kinFiltered[!((kinFiltered$ID1%in%theseNames)&(kinFiltered$class=="3rd degree")),]
kinFiltered = kinFiltered[!((kinFiltered$ID2%in%theseNames)&(kinFiltered$class=="3rd degree")),]


# excluding those with >10 3rd-degree relatives after pruning 200
write.table(kinFiltered[,c("FID1","ID1","FID2","ID2","N_SNP","HetHet","IBS0","Kinship")],file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-samples-pruned.kin0"),quote=FALSE,row.names=FALSE,col.names=TRUE )

# exclude het.miss outliers, plus anyone with over 200 3rd-degree relatives ==> used for output.
write.table(kin2[,c("FID1","ID1","FID2","ID2","N_SNP","HetHet","IBS0","Kinship")],file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-samples-pruned-",maxNrelatives,".kin0"),quote=FALSE,row.names=FALSE,col.names=TRUE )


# print out list of samples that were effectively excluded from the kinship calculation - i.e hetmiss outliers + those that had more than 200 3rd-degree relatives after removing outliers.
    
write.table(allExcludedNames,file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-samples-pruned-",maxNrelatives,"-sample-effective_exclusions.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE ) # Note that this file will actually include some samples that weren't actually in the raw kinship file, but don't want to trust that they're NOT related to anyone...

write.table(tooManyRelatives,file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-samples-pruned-",maxNrelatives,"-sample-gt",maxNrelatives,"_relatives.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE )  # this is the list after excluding the het-miss outliers

write.table(unique(names(totalKin10)),file=paste0(baseSampleQCDir,"/data/Relatedness/",OutputFile,"-samples-pruned-",maxNrelatives,"-sample-gt10_3rdDeg_relatives.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE )  # this is the list after excluding the het-miss outliers and anyone with more than 200 relatives.

