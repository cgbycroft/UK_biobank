# This script takes the output from filter-king.R and looks for trios. It outputs a table with these groups for use in imputation. It also checks the ages of samples in trios etc. as confirmation.

args = commandArgs(TRUE)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/ForRelease/b1__b11-b001__b095-sampleTable_v4_Kinship.txt","-ibs0","0.0012")

print(args)
h = args[-c(which(args%in%c("-in","-out","-ibs0","-exclude-dupes","-3degree")),1+which(args%in%c("-in","-out","-ibs0","-exclude-dupes","-3degree")))]
for(helperScript in h){
    source(helperScript)
}

library(igraph)
library(plotrix)

KingOutFile = args[which(args=="-in")+1]
OutputFile = gsub("\\.kin0|\\.txt","",basename(KingOutFile))
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness"))

## get kinship data
kin = read.table(KingOutFile,header=T,stringsAsFactors=F)

## get other data
otherInfo = read.multiple.batch.info(batchInfoFields[c(1:59)])
otherInfo = tbl_df(otherInfo)

## get list of duplicates
dupFile = "/well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt"
dup = read.table(dupFile,header=FALSE)[,1]

## get kinship classes
kin$class = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)
table(kin$class)
kin$col = colDef[kin$class]

kin = kin[,!colnames(kin)%in%c("FID1","FID2")]

if("-3degree"%in%args) {
    # exclude any pairs with kinship below given threshold
    
    deg3Thresh = as.numeric(args[which(args=="-3degree")+1])
    OutputFile = paste0(OutputFile,"-3deg-",deg3Thresh)
    print( paste0(sum(!kin$Kinship>=deg3Thresh)," pairs excluded because have kinship below ",deg3Thresh))
    nomore=kin[!kin$Kinship>=deg3Thresh,]
    print( paste0("This comprises ",length(unique(c(nomore$ID1,nomore$ID2)))," individuals.") )
    
    kin=kin[kin$Kinship>=deg3Thresh,]
        
}

if("-exclude-dupes" %in% args){
    print( "This run excludes duplicates" )
    otherInfo = filter(otherInfo,!PIID%in%dup)
    kin = kin[(!kin$ID1 %in% dup )&(!kin$ID2 %in% dup ),]
    OutputFile = paste0(OutputFile,"-nodupes")
}


print("Dimensions of kinship table:")
print( dim(kin) )

## get pairwise network
# including 3rd degree
network <- graph.data.frame(kin, directed=F)
V(network)$size <- log(degree(network))+1
V(network)$color <- otherInfo$Colors[match(V(network)$name,otherInfo$PIID)]
E(network)$color <- E(network)$col
fams = clusters(network)


# excluding 3rd degree
kinNo3 = kin[kin$class!="3rd degree",]

networkNo3 <- graph.data.frame(kinNo3, directed=F)
V(networkNo3)$size <- log(degree(networkNo3))+1
V(networkNo3)$color <- otherInfo$Colors[match(V(networkNo3)$name,otherInfo$PIID)]
E(networkNo3)$color <- E(networkNo3)$col
famsNo3 = clusters(networkNo3)



## Tabulate family configurations (if needed) for families of size > 2

kin$class2 = factor(kin$class)
kinNo3$class2 = factor(kinNo3$class)

fams3= which(fams$csize > 2)
famsNo33= which(famsNo3$csize > 2)


#famSaveFile = gsub(".kin0","-families.RData",basename(KingOutFile))
famSaveFile = paste0(OutputFile,"-families.RData")
if(basename(famSaveFile)%in%list.files(dirname(famSaveFile))){
    load(file=famSaveFile,verbose=TRUE)
} else {

    # takes about 1 second per 50 families
    # including 3rd degree relatives
    fams3Configs = sapply(1:length(fams3),function(i){
        if(i%%100==0) print(i)
        x=fams3[i]
        inds=names(fams$membership)[fams$membership==x]
        config = table(kin$class2[(kin$ID1%in%inds)|(kin$ID2%in%inds)]) # count of pairwise connections.
    })
    colnames(fams3Configs) = fams3

    # excluding 3rd degree relatives (some families drop out altogether, and some will split up.)
    famsNo33Configs = sapply(1:length(famsNo33),function(i){
        if(i%%100==0) print(i)
        x=famsNo33[i]
        inds=names(famsNo3$membership)[famsNo3$membership==x]
        config = table(kinNo3$class2[(kinNo3$ID1%in%inds)|(kinNo3$ID2%in%inds)]) # count of pairwise connections.
    })
    colnames(famsNo33Configs) = famsNo33
    
# save this filtered family configurations file
    save(networkNo3,network,famsNo3,fams,fams3Configs,famsNo33Configs,file=famSaveFile)
}


## PLOTTING 

## plot family sizes
png(paste0("plots/",OutputFile,"-FamilySizes.png"),width=2000,height=1500,res=150)
par(mfrow=c(2,1),mar=c(5,4,1,1))
barplot(table(fams$csize),main="Including 3rd degree",xlab="Family size",ylab="Count of families")
barplot(table(famsNo3$csize),main="Excluding 3rd degree",xlab="Family size",ylab="Count of families")
dev.off()

png(paste0("plots/",OutputFile,"-FamilySizes-Log.png"),width=2000,height=1500,res=150)
par(mfrow=c(2,1),mar=c(5,4,1,1))
barplot(table(fams$csize),main="Including 3rd degree",log="y",xlab="Family size",ylab="Count of families",col="blue")
barplot(table(famsNo3$csize),main="Excluding 3rd degree",log="y",xlab="Family size",ylab="Count of families",col="blue")
dev.off()


# write out table of duplicates/twins
write.table(kin[kin$class=="dupe/twins",],file=paste0(OutputFile,"-duplicates-twins.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)
print( paste0(sum(kin$class=="dupe/twins") ," twins or duplicates written to:"))
print( paste0(OutputFile,"-duplicates-twins.txt") )    
    

## Plot twins - more than two twins in a family ----> duplicates
theseFams = fams3[fams3Configs["dupe/twins",]>1]
inds = names(fams$membership)[fams$membership%in%theseFams]
if(length(inds)>0){
    subGraph = igraph::induced_subgraph(network,inds,impl="copy_and_delete")

    png(paste0("plots/",OutputFile,"-families-morethan2twins-graph.png"),width=1500,height=1500,res=150)
    plot.igraph(subGraph,main=paste0("Families with more than 2 twins"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle")
    dev.off()
}


## Plot examples of families
print("Plotting examples of families graphs...")

for(famSize in c(3:9)){

    print(famSize)
    
    theseFams = which(fams$csize==famSize)
    inds = names(fams$membership)[fams$membership%in%theseFams]
    subGraph = igraph::induced_subgraph(network,inds,impl="copy_and_delete")
    
    configurations = fams3Configs[c("dupe/twins","parent/child","sibs","2nd degree","3rd degree"),as.character(theseFams)]
    strings = apply(configurations,2,function(x) paste(x,collapse=""))
    stringsU = unique(strings)
    
    png(paste0("plots/",OutputFile,"-families-size-",famSize,"-barplot-configs.png"),width=1500,height=1500,res=150)
    confTab = table(strings)
    barplot(sort(confTab,decreasing=TRUE),las=2,main=paste0("Configurations of families with ",famSize," members"))
    legend("topright",legend=c(rev(rownames(configurations))),title="configuration order",bty="n")
    dev.off()
    
    # select random sample of each configuration
    props = confTab/sum(confTab)
    f = c()
    for(config in names(confTab)){
        n = round(50*props[config],0)
        if(n==0) n=1
        if(sum(strings==config)==1) f = c(f,colnames(configurations)[strings==config]) else f = c(f,sample(colnames(configurations[,strings==config]),n,replace=FALSE))        
    }

    inds = names(fams$membership)[fams$membership%in%as.numeric(f)]
    subGraph = igraph::induced_subgraph(network,inds,impl="copy_and_delete")
    
    png(paste0("plots/",OutputFile,"-families-size-",famSize,"-sample-graph.png"),width=1500,height=1500,res=150)
    plot.igraph(subGraph,main=paste0("Families with ",famSize," members\n",length(confTab)," unique configurations"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle")
    dev.off()
}


## Plot large families
theseFams = which(fams$csize>10)
inds = names(fams$membership)[fams$membership%in%theseFams[1:(round(length(theseFams)/2,0))]]
subGraph1 = igraph::induced_subgraph(network,inds,impl="copy_and_delete")
famIDs1 = getFamLabels(inds,fams,subGraph1)

inds = names(fams$membership)[fams$membership%in%theseFams[(round(length(theseFams)/2,0) +1 ):length(theseFams)]]
subGraph2 = igraph::induced_subgraph(network,inds,impl="copy_and_delete")
famIDs2 = getFamLabels(inds,fams,subGraph2)

png(paste0("plots/",OutputFile,"-families-morethan10members-graph-%02d.png"),width=1500,height=1500,res=150)
# with labels
plot.igraph(subGraph1,main=paste0("Families with more than 10 members"),vertex.frame.color=NA,vertex.label=famIDs1,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight,vertex.label.color="black",vertex.label.family="Arial",vertex.label.dist=0.5)
plot.igraph(subGraph2,main=paste0("Families with more than 10 members"),vertex.frame.color=NA,vertex.label=famIDs2,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight,vertex.label.color="black",vertex.label.family="Arial",vertex.label.dist=0.5)
# without labels
plot.igraph(subGraph1,main=paste0("Families with more than 10 members"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight)
plot.igraph(subGraph2,main=paste0("Families with more than 10 members"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight)
dev.off()

                                
## Plot large families > size 5 with 1,2nd degree only; but plot the 3rd degree relatives too.
theseFams = which(famsNo3$csize>6)
inds = names(famsNo3$membership)[famsNo3$membership%in%theseFams]
indsRels=unlist(sapply(inds,function(x) V(network)$name[neighbors(network, x)]))
subGraph1 = igraph::induced_subgraph(network,c(inds,indsRels),impl="copy_and_delete")
famIDs1 = getFamLabels(inds,famsNo3,subGraph1)

theseFams = which(famsNo3$csize==6)
inds = names(famsNo3$membership)[famsNo3$membership%in%theseFams]
indsRels=unlist(sapply(inds,function(x) V(network)$name[neighbors(network, x)]))
subGraph2 = igraph::induced_subgraph(network,c(inds,indsRels),impl="copy_and_delete")
famIDs2 = getFamLabels(inds,famsNo3,subGraph2)

for(class in names(kinWeights)){
    # apply edge-weights
    E(subGraph1)$weight[E(subGraph1)$class%in%class] <- kinWeights[class]/3
    E(subGraph2)$weight[E(subGraph2)$class%in%class] <- kinWeights[class]/3
}

png(paste0("plots/",OutputFile,"-families-no3degree-morethan6members-graph-%02d.png"),width=1500,height=1500,res=150)
# with labels
plot.igraph(subGraph1,main=paste0("Families with 7 or more members \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=famIDs1,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight,vertex.label.color="black",vertex.label.family="Arial",vertex.label.dist=0.5)
plot.igraph(subGraph2,main=paste0("Families with 6 members \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=famIDs2,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph2)$weight,vertex.label.color="black",vertex.label.family="Arial",vertex.label.dist=0.5)
# without labels
plot.igraph(subGraph1,main=paste0("Families with 7 or more members \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight)
plot.igraph(subGraph2,main=paste0("Families with 6 members \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph2)$weight)

dev.off()


#### Plot mixtures of different family configurations
# pairs
fams3Inds = names(fams$membership)[fams$membership%in%which(fams$csize > 2)]
take = sample(which((!kin$class %in% "3rd degree")&(!kin$ID1 %in% fams3Inds)&(!kin$ID2 %in% fams3Inds)),5,replace=FALSE)
pairInds = unique(c(kin$ID1[take],kin$ID2[take]))

# trios
trioFams = colnames(fams3Configs)[(fams3Configs["parent/child",]==2) & (colSums(fams3Configs)==2)]
trioInds = names(fams$membership)[fams$membership%in%as.numeric(trioFams)]

# large sib families
sibFams = colnames(fams3Configs)[(fams3Configs["sibs",] > 4)]
sibInds = names(fams$membership)[fams$membership%in%as.numeric(sibFams)]

# sibs and parents
sibFams2 = colnames(fams3Configs)[(fams3Configs["sibs",] > 1)&(fams3Configs["parent/child",] > 1 )]
sibInds2 = names(fams$membership)[fams$membership%in%as.numeric(sibFams2)]

# networks (no 3rd degree)
netFams = which(famsNo3$csize > 5)
netInds = names(famsNo3$membership)[famsNo3$membership%in%netFams]

# twins
twinFams = colnames(fams3Configs)[(fams3Configs["dupe/twins",]==1) & (colSums(fams3Configs)>1)]
twinInds = names(fams$membership)[fams$membership%in%as.numeric(twinFams)]


set.seed(12345)
f = c(sample(trioFams,4,replace=FALSE),sample(sibFams2,4,replace=FALSE),sample(sibFams,4,replace=FALSE),sample(netFams,3,replace=FALSE),twinFams[1:3])
inds = c(pairInds, names(fams$membership)[fams$membership%in%as.numeric(f)] )

subGraph = igraph::induced_subgraph(network,inds,impl="copy_and_delete")

E(subGraph)$weight <- 1
for(class in names(kinWeights)){
    E(subGraph)$weight[E(subGraph)$class%in%class] <- kinWeights[class]/3.5
}
#E(subGraph)$color[E(subGraph)$class%in%c("dupe/twins")] = "yellow"
#E(subGraph)$color[E(subGraph)$class%in%c("2nd degree")] = "darkgreen"

#l <- layout_with_kk(subGraph)

png(paste0("plots/",OutputFile,"-example-families-graph-%02d.png"),width=1500,height=1500,res=150,bg="transparent")
plot.igraph(subGraph,main=paste0("Example families"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",vertex.color="black",edge.width=E(subGraph)$weight)

dev.off()


#### Compute family network connectedness...

# just compute over 3 or more families
connectivityNo3 = sapply(famsNo33,function(f){
    print(f)
    inds = names(famsNo3$membership)[famsNo3$membership==f]
    subGraph = igraph::induced_subgraph(networkNo3,inds,impl="copy_and_delete")
    con = edge_connectivity(subGraph)
    return(con)
})
names(connectivityNo3) = famsNo33

highConnectivityNo3 = as.numeric(names(connectivityNo3)[connectivityNo3>3])
#famsNo3$csize[highConnectivityNo3]

# plot the high connectivity families
inds = names(famsNo3$membership)[famsNo3$membership%in%highConnectivityNo3]
indsRels=unlist(sapply(inds,function(x) V(network)$name[neighbors(network, x)]))
subGraph1 = igraph::induced_subgraph(network,c(inds,indsRels),impl="copy_and_delete")
famIDs1 = getFamLabels(inds,famsNo3,subGraph1)
for(class in names(kinWeights)){
    # apply edge-weights
    E(subGraph1)$weight[E(subGraph1)$class%in%class] <- kinWeights[class]/3
}


png(paste0("plots/",OutputFile,"-families-no3degree-gt3-connectivity-graph-%02d.png"),width=1500,height=1500,res=150)
# with labels
plot.igraph(subGraph1,main=paste0("Families with connectivity 4 or greater \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=famIDs1,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight,vertex.label.color="black",vertex.label.family="Arial",vertex.label.dist=0.5)
# without labels
plot.igraph(subGraph1,main=paste0("Families with connectivity 4 or greater \n connected 1st or 2nd degree"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle",edge.width=E(subGraph1)$weight)

dev.off()




#### Plot a legend for the colours of lines in family plots
png(paste0("plots/",OutputFile,"-families-graph-LEGEND.png"),bg="transparent",width=1500,height=1500,res=150)
par(cex=3)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab=NA,ylab=NA)
ord = order(kinWeights,decreasing=TRUE)
legend("top",legend=names(colDef)[ord],col=colDef[ord],lwd=kinWeights[ord])
dev.off()


print( "DONE!" )
print( paste0("You can now run find-trios.R to create a list of trios and other families for phasing.") ) 
