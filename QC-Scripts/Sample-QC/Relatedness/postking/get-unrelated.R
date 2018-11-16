## Script to take king output and get maximal set of unrelated samples
if(!"args"%in%ls()) args = commandArgs(trailingOnly=T)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0","-ibs0","0.002", "-outdir" ,"postking","-outfile","")

library(igraph)

print(args)

h = args[-c(which(args%in%c("-in","-ibs0","-outdir","-exclude","-outfile")),1+which(args%in%c("-in","-ibs0","-outdir","-exclude","-outfile")))]
print(h)

for(helperScript in h){
    source(helperScript)
}

#quit()

KingOutFile = args[which(args=="-in")+1]
if("-outfile"%in%args) OutputFile = args[which(args=="-outfile")+1] else OutputFile = gsub(".kin0","",basename(KingOutFile))
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])
OutputDir = args[which(args=="-outdir")+1]


# get all sample names (default read.multiple.batch.info is without references)
otherInfo = read.multiple.batch.info(c("PIID","Pops"))

# read kinship file
kin = read.table(KingOutFile,header=TRUE,stringsAsFactors=FALSE)

# get kinship classes
kin$class = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)
kin$col = colDef[kin$class]
table(kin$class)

# unfiltered
#  2nd degree   3rd degree   dupe/twins parent/child         sibs 
#       11253        62491          890         6280        22737
# filtered
#  2nd degree   3rd degree   dupe/twins parent/child         sibs 
#       11191        61893          891         6318        22789


# exclude samples with 100s of 3rd-degree 'relatives'
nRelatives = table(c(kin$ID1[kin$class=="3rd degree"],kin$ID2[kin$class=="3rd degree"]))
extreme = names(nRelatives)[nRelatives > 100]

extremeFile = paste0(OutputDir,"/",OutputFile,"-gt100-relatives.txt")
print(paste0(length(extreme)," samples have more than 100 relatives. List of these samples saved to: ",extremeFile))
write.table(extreme,file=extremeFile,quote=FALSE,row.names=FALSE,col.names=FALSE)


# exclude any samples requested and add to extreme list
if("-exclude"%in%args) {
    exFiles = args[which(args=="-exclude")+1]
    exclusions = c()

    for(file in exFiles){
        print( paste0("Reading: ",file))
        e = read.table(file,stringsAsFactors=FALSE)[,1]
        print( paste0( length(e)," samples.") )
        exclusions = c(exclusions,e)
    }

    exclusions = unique(exclusions)
    
    print( paste0(sum(exclusions%in%extreme)," samples in exclude list already have more than 100 relatives.") )
    print( paste0(sum(!extreme%in%exclusions)," samples with more than 100 relatives that are not otherwise in the exclusion list." ) )
    print( paste0(length(exclusions)-sum(exclusions%in%extreme)," extra samples added to the samples to exclude.") )
    extreme = unique(c(exclusions,extreme))
}

#load(file=paste0("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-pruned.kinship.RData"),verbose=TRUE)

kin2 = kin[!((kin$ID1%in%extreme) | (kin$ID2%in%extreme)),c(-1,-3)]


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

# get list of samples not in filtered kinship table, nor in extreme
unrelatedEver = otherInfo$PIID[(!otherInfo$PIID%in%kin2$ID1) & (!otherInfo$PIID%in%kin2$ID2) & (!otherInfo$PIID%in%extreme)]

unrelatedEverNo3 = otherInfo$PIID[(!otherInfo$PIID%in%kin2No3$ID1) & (!otherInfo$PIID%in%kin2No3$ID2) & (!otherInfo$PIID%in%extreme)]


# compute maximal set on families less than size 40, and greater than size 2
# igraph uses algorithm described by S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm for generating all the maximal independent sets. SIAM J Computing, 6:505–517, 1977.
print("Computing maximal set of unrelated samples...")

# INCLUDING 3rd DEGREE
fams3= which(fams$csize > 2)
fam3inds = names(fams$membership)[fams$membership%in%fams3]
networkFam3 = igraph::induced_subgraph(network,fam3inds,impl="copy_and_delete")
h3 = clusters(networkFam3)

# get list of samples to keep in maximal independent set. Takes around 5mins to run
incl = sapply(1:h3$no,function(x){
    if(x%%100==0) print(round(100*x/h3$no,0) )
    f = igraph::induced_subgraph(networkFam3,names(h3$membership)[h3$membership==x],impl="auto")
    if(length(V(f)) < 40) {
        li = largest_ivs(f)
        if(length(li)>1) li=li[sample(1:length(li),1)] # randomly sample from one solution
        l = unlist(li)
    } else {
        print(paste0("too big... " ,x))
        l=NULL
    }
    return(l)
})

incl = names(unlist(incl))

# for the 2-person families include one at random
fam2inds = names(fams$membership)[!fams$membership%in%fams3]
this = kin2[(kin2$ID1%in%fam2inds)&(kin2$ID2%in%fam2inds),]
s = sample(c(TRUE,FALSE),nrow(this),replace=TRUE)
S = cbind(s,!s)
inclInds2 = c(this[S[,1],"ID1"],this[S[,2],"ID2"])

unrelated = unique(c(incl,inclInds2,unrelatedEver))


# check they really are unrelated
if( sum( (kin$ID1%in%unrelated) & (kin$ID2%in%unrelated) ) > 0 ) print("WARNING: chosen samples actually related!!")

# total samples we would have to exclude:
print(paste0("Total unrelated including 3rd degree: ",length(unrelated)))




# EXCLUDING 3rd DEGREE
fams3= which(famsNo3$csize > 2)
fam3inds = names(famsNo3$membership)[famsNo3$membership%in%fams3]
networkFam3 = igraph::induced_subgraph(network,fam3inds,impl="copy_and_delete")
h3 = clusters(networkFam3)

inclNo3 = sapply(1:h3$no,function(x){
    if(x%%100==0) print(round(100*x/h3$no,0) )
    f = igraph::induced_subgraph(networkFam3,names(h3$membership)[h3$membership==x],impl="auto")
    if(length(V(f)) < 40) {
        li = largest_ivs(f)
        if(length(li)>1) li=li[sample(1:length(li),1)] # randomly sample from one solution
        l = unlist(li)
    } else {
        print(paste0("too big... " ,x))
        l=NULL
    }
    return(l)
})

inclNo3 = names(unlist(inclNo3))

# for the 2-person families include one at random
fam2inds = names(famsNo3$membership)[!famsNo3$membership%in%fams3]
this = kin2[(kin2$ID1%in%fam2inds)&(kin2$ID2%in%fam2inds),]
s = sample(c(TRUE,FALSE),nrow(this),replace=TRUE)
S = cbind(s,!s)
inclInds2 = c(this[S[,1],"ID1"],this[S[,2],"ID2"])

unrelatedNo3 = unique(c(inclNo3,inclInds2,unrelatedEverNo3))

                                        # check they really are unrelated
if( sum( (kin$ID1[kin$class!="3rd degree"]%in%unrelatedNo3) & (kin$ID2[kin$class!="3rd degree"]%in%unrelatedNo3) ) > 0 ) print("WARNING: chosen samples actually related!!")


# total samples we would have to exclude:
print(paste0("Total unrelated ignoring 3rd degree: ",length(unrelatedNo3)))




# save list of unrelated samples
print( paste0("Saving lists of unrelated samples to ",OutputDir,"/",OutputFile,"-unrelated*.txt") )

write.table(unrelated,file=paste0(OutputDir,"/",OutputFile,"-unrelated.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(unrelatedNo3,file=paste0(OutputDir,"/",OutputFile,"-unrelated-ignoreDegree3.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
