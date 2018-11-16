# This script takes the output from filter-king.R and looks for trios. It outputs a table with these groups for use in imputation. It also checks the ages of samples in trios etc. as confirmation.

args = commandArgs(TRUE)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/ForRelease/b1__b11-b001__b095-sampleTable_v3_Kinship.txt","-ibs0","0.0012")

print(args)
h = args[-c(which(args%in%c("-in","-out","-ibs0","-exclude-dupes","-3degree")),1+which(args%in%c("-in","-out","-ibs0","-exclude-dupes","-3degree")))]
for(helperScript in h){
    source(helperScript)
}

library(igraph)
setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness"))

KingOutFile = args[which(args=="-in")+1]
OutputFile = gsub("\\.kin0|\\.txt","",basename(KingOutFile))
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])


## get kinship data
kin = read.table(KingOutFile,header=T,stringsAsFactors=F)
kin$class = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)
table(kin$class)
kin$col = colDef[kin$class]


if("-3degree"%in%args) {
    # exclude any pairs with kinship below given threshold
    
    deg3Thresh = as.numeric(args[which(args=="-3degree")+1])
    OutputFile = paste0(OutputFile,"-3deg-",deg3Thresh)
    print( paste0(sum(!kin$Kinship>=deg3Thresh)," pairs excluded because have kinship below ",deg3Thresh))
    nomore=kin[!kin$Kinship>=deg3Thresh,]
    print( paste0("This comprises ",length(unique(c(nomore$ID1,nomore$ID2)))," individuals.") )
    
    kin=kin[kin$Kinship>=deg3Thresh,]
        
}

kinshipNo3 = kin[kin$class!="3rd degree",]

## get other data
otherInfo = read.multiple.batch.info(batchInfoFields[c(1:59)])
otherInfo = tbl_df(otherInfo)


## get list of duplicates
dupFile = "/well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt"
dup = read.table(dupFile,header=FALSE)[,1]

if("-exclude-dupes" %in% args){
    print( "This run excludes duplicates" )
    otherInfo = filter(otherInfo,!PIID%in%dup)
    kin = kin[(!kin$ID1 %in% dup )&(!kin$ID2 %in% dup ),]
    OutputFile = paste0(OutputFile,"-nodupes")
    
    kinshipNo3 = kin[kin$class!="3rd degree",]
}


## get family connection data (this is created in filter-king.R)
load(paste0(OutputFile,"-families.RData"),verbose=TRUE)



########################
#### FIND CLEAN TRIOS (IGNORING 3rd DEGREE)
# write a list of clean trios - i.e have no other family members in the data
# "clean" trios are families that look like: 2 parent/child connections, and nothing else. Use age and sex information to confirm trios.

print("Getting list of trios...")

system(paste0("mkdir forPhasing"))


## Ignore 3rd degree relatives in trio
nConnections = colSums(famsNo33Configs)
trios = (nConnections==2)&(famsNo33Configs["parent/child",]==2) # parent/child families with only 2 connections indicates trio, because otherwise there would also be a connection for the two siblings (either full sibs or half sibs).
trioFams = as.numeric(colnames(famsNo33Configs))[trios]
nTrios = sum(trios)
# 724 trios (before validation)

# check they're all families of size 3
sum(famsNo3$csize[trioFams]!=3)

# write table of trios. One row for each, and a sample ID
#  Mum Dad Child
inds = names(famsNo3$membership)[famsNo3$membership%in%trioFams]
kins = kin[(kin$ID1%in%inds)&(kin$ID2%in%inds),]
nRel = table(c(kins$ID1,kins$ID2))
children = names(nRel)[nRel ==2]
aye1 = match(kins$ID1,otherInfo$PIID); aye2 = match(kins$ID2,otherInfo$PIID);
age1 = otherInfo$Age.when.attended.assessment.centre[aye1]
age2 = otherInfo$Age.when.attended.assessment.centre[aye2]
sex1 = as.character(otherInfo$Inferred.Gender[aye1])
sex2 = as.character(otherInfo$Inferred.Gender[aye2])

trioIDs = t(sapply(children, function(i){
    wrong=FALSE
    this = (kins$ID1==i)|(kins$ID2==i)
    ids = c(kins$ID1[this],kins$ID2[this])
    adult = which(ids!=i)
    A = c(age1[this],age2[this])
    S = c(sex1[this],sex2[this])
    if(sum(S[adult]=="F")!=1) {wrong=TRUE; print("Too few, or too many female 'parents'") }
    if(sum(S[adult]=="M")!=1) {wrong=TRUE; print("Too few, or too many male 'parents'") }
    if(sum( (A[adult] - A[-adult]) < 10 )  > 0 ) {wrong=TRUE; print("At least one 'adult' less than 10 years older than 'child'.") }
    Mum = ids[adult][S[adult]=="F"]
    Dad = ids[adult][S[adult]=="M"]
    if(wrong){
        print(i)
        return(c(i,NA,NA))
    } else {
        return(c(i,Mum,Dad))
    }
 }))

# Issues can occur only if real sibling connections aren't found, or where grand-parent/grandchild connections aren't found. Maybe these are called as 3rd-degree instead of 2nd degree?

# print list of weird trios
w = trioIDs[is.na(trioIDs[,2]),1]
weird = unique(c(kin$ID1[(kin$ID1%in%w)|(kin$ID2%in%w)],kin$ID2[(kin$ID1o%in%w)|(kin$ID2%in%w)])) # list individuals in pairs with the weird samples
print( paste0(length(weird)," individuals involved in odd trio configurations. Excluding these from trio list.") )

write.table(cbind(weird,weird),file=paste0(OutputFile,"-odd-trios.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)

# ====> run king separately on these odd trios

trioIDs = trioIDs[!is.na(trioIDs[,2]),] # remove any false ones
colnames(trioIDs)=c("Child","Mother","Father")
# 724 non-weird trios


# check M/F and Age again
if( sum(otherInfo$Inferred.Gender[match(trioIDs[,"Mother"],otherInfo$PIID)]!="F") !=0 ) print("Something wrong with trio set. Not all mothers are coded Female.")
if( sum(otherInfo$Inferred.Gender[match(trioIDs[,"Father"],otherInfo$PIID)]!="M") !=0 ) print("Something wrong with trio set. Not all mothers are coded Female.")
ageChild = otherInfo$Age.when.attended.assessment.centre[match(trioIDs[,"Child"],otherInfo$PIID)]
ageM = otherInfo$Age.when.attended.assessment.centre[match(trioIDs[,"Mother"],otherInfo$PIID)]
ageF = otherInfo$Age.when.attended.assessment.centre[match(trioIDs[,"Father"],otherInfo$PIID)]
if(min(ageM - ageChild) < 12) print("At least one mother with an age-difference less than 10 with her child...")
if(min(ageF - ageChild) < 12) print("At least one father with an age-difference less than 10 with her child...")

print(paste0(dim(trioIDs)[1]," clean trios (ignoring 3rd degree) found and validated with age and sex."))

write.table(trioIDs,file=paste0("forPhasing/",OutputFile,"-clean-trios-ignoreDegree3.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)



########################
#### FIND CLEAN TRIOS (ACCOUNTING FOR 3rd DEGREE)

# Trios with no 3rd degree attached (a subset of trioIDs)
nConnectionsAll = colSums(fams3Configs)
triosAll = (colSums(fams3Configs)==2)&(fams3Configs["parent/child",]==2)
trioAllFams = as.numeric(colnames(fams3Configs))[triosAll]
inds2 = names(fams$membership)[fams$membership%in%trioAllFams]

print( paste0( sum(weird%in%inds2)," weird trio samples still with no 3rd degree connections... "))

trioIDs2 = trioIDs[(trioIDs[,1]%in%inds2)&(trioIDs[,2]%in%inds2)&(trioIDs[,3]%in%inds2),]

write.table(trioIDs2,file=paste0("forPhasing/",OutputFile,"-clean-trios.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)


########################
#### ALL FAMILIES (IGNORING 3rd DEGREE)

# Other families:
# First excluding 3rd degree
kinshipNo3$familyID1 = famsNo3$membership[kinshipNo3$ID1]
kinshipNo3$familyID2 = famsNo3$membership[kinshipNo3$ID2]

if(sum(kinshipNo3$familyID1!=kinshipNo3$familyID2)!=0) print("Something wrong with kinship family ids...")

trioInds = unique(c(trioIDs[,1],trioIDs[,2],trioIDs[,3]))
inTrios = apply(kinshipNo3,1,function(x){
    #print(x)
    to = x[[1]];  from = x[[2]]
    if((to%in%trioInds) & (from%in%trioInds)) {
        
        testTo = c(which(trioIDs[,1]==to),which(trioIDs[,2]==to),which(trioIDs[,3]==to))
        testFrom = c(which(trioIDs[,1]==from),which(trioIDs[,2]==from),which(trioIDs[,3]==from))
        if(length(intersect(testTo,testFrom))>0) return(1) else return(0)
        
    } else {
        return(0)
    }        
})

# inTrio flag indicates any pair that are part of a clean trio, ignoring 3rd degree (trioIDs)
kinshipNo3$inTrio = inTrios

#write.table(kinshipNo3[,c("ID1","ID2","Kinship","class","familyID1","inTrio")],file=paste0("forPhasing/",OutputFile,"-kinship-ignoreDegree3.txt"),quote=FALSE,col.names=c("ID1","ID2","Kinship","Class","FamilyID","InTrio"),row.names=FALSE,sep="\t")


########################
#### ALL FAMILIES (ACCOUNTING FOR 3rd DEGREE)

# Then including 3rd degree
kin$familyID1 = fams$membership[kin$ID1]
kin$familyID2 = fams$membership[kin$ID2]

if(sum(kin$familyID1!=kin$familyID2)!=0) print("Something wrong with kinship family ids...")

trioInds = unique(c(trioIDs2[,1],trioIDs2[,2],trioIDs2[,3]))
inTrios = apply(kin,1,function(x){
    #print(x)
    to = x[[1]];  from = x[[2]]
    if((to%in%trioInds) & (from%in%trioInds)) {
        #print(which(kin$ID1==x[1]))
        testTo = c(which(trioIDs2[,1]==to),which(trioIDs2[,2]==to),which(trioIDs2[,3]==to))
        testFrom = c(which(trioIDs2[,1]==from),which(trioIDs2[,2]==from),which(trioIDs2[,3]==from))
        if(length(intersect(testTo,testFrom))>0) return(1) else return(0)
        
    } else {
        return(0)
    }        
})

# inTrio flag indicates any pair that are part of a clean trio, accounting for 3rd degree (trioIDs2)
kin$inTrio = inTrios

#write.table(kin[,c("ID1","ID2","Kinship","class","familyID1","inTrio")],file=paste0("forPhasing/",OutputFile,"-kinship.txt"),quote=FALSE,col.names=c("ID1","ID2","Kinship","Class","FamilyID","InTrio"),row.names=FALSE,sep="\t")



##########################################
###### FIND ALL NON-CLEAN and CLEAN TRIOs

## Find all non-clean trios (i.e allow other relatives involved)
# distinctive pattern is: at least 2 parent/child pairs in family, with
# for simplicity, just exclude 3rd-degree pairs, then double-check ages and sex

# 2 parent/child pairs:
pc2 = famsNo33Configs[,famsNo33Configs["parent/child",] >= 2] # 1240 groups with >=2 parent/child pairs. Maximum is 4 pairs

# For each parent/child pair, check if the child has two parents. The child is the younger of the two (assuming ages are correct)
inds = names(famsNo3$membership)[famsNo3$membership%in%as.numeric(colnames(pc2))]
pc2kins = kinshipNo3[(kinshipNo3$ID1%in%inds)|(kinshipNo3$ID2%in%inds),]
aye1 = match(pc2kins$ID1,otherInfo$PIID); aye2 = match(pc2kins$ID2,otherInfo$PIID);
age1 = otherInfo$Age.when.attended.assessment.centre[aye1]
age2 = otherInfo$Age.when.attended.assessment.centre[aye2]
sex1 = as.character(otherInfo$Inferred.Gender[aye1])
sex2 = as.character(otherInfo$Inferred.Gender[aye2])

# children in trios will always appear more than twice in the table in a parent/child pair

pcPairs = which(pc2kins$class=="parent/child")
pcPairsKin = pc2kins[pcPairs,]
nPairs = table(c(pcPairsKin$ID1,pcPairsKin$ID2))

toCheck = names(nPairs)[nPairs>=2]

trioList = t(sapply(toCheck,function(i){
    ourPairs = which(((pc2kins$ID1==i)|(pc2kins$ID2==i))&(pc2kins$class=="parent/child"))
    pids = c(pc2kins$ID1[ourPairs],pc2kins$ID2[ourPairs])
    pids = pids[pids!=i]
    ageChild = otherInfo$Age.when.attended.assessment.centre[otherInfo$PIID==i]
    ageParents = otherInfo$Age.when.attended.assessment.centre[match(pids,otherInfo$PIID)]
    sexParents = as.character(otherInfo$Inferred.Gender[match(pids,otherInfo$PIID)]    )
    ageDiff = ageParents - ageChild
    if(sum(ageDiff > 10) >= 2 ) {
        child = TRUE
        parents = pids[ageDiff>10]
        
        # are any of the 'parents' actually related or duplicates?
        relParents = pc2kins[(pc2kins$ID1%in%parents)&(pc2kins$ID2%in%parents),]
         
        if(dim(relParents)[1]>0){
            print(relParents)
            print("two 'parents' are related. Check this")
            print(i)
                                        # can we distinguish which are the two parents by looking at sex?
            #sex = otherInfo$Inferred.Gender[match(c(relParents$ID1,relParents$ID2),otherInfo$PIID)]
           # print(i,sex)
            #out = c(i,pids[sexParents=="F"],pids[sexParents=="M"])
            out = c(i,NA,"relParents")
        } else {
            out = c(i,pids[sexParents=="F"],pids[sexParents=="M"])
        }
            
    } else {
        print(i)
        print("not two parents")
        out = c(i,NA,NA)
    }
    return(out)
}))


# check that these are a superset of the other trio sets.
trioString = paste0(trioList[,1],trioList[,2],trioList[,3])
trioStringClean = paste0(trioIDs[,1],trioIDs[,2],trioIDs[,3])
trioStringClean2 = paste0(trioIDs2[,1],trioIDs2[,2],trioIDs2[,3])

if( length(intersect(trioString,trioStringClean)) != dim(trioIDs)[1] ) print("unclean trio list is not a superset of clean trio list!")
if( length(intersect(trioString,trioStringClean2)) != dim(trioIDs2)[1] ) print("unclean trio list is not a superset of clean (ignoring 3rd degree) trio list!")

weirdUnclean = trioList[which((is.na(trioList[,2])&(trioList[,3]=="relParents"))),1]

trioListSet = trioList[!is.na(trioList[,2]),]
colnames(trioListSet) = colnames(trioIDs)

# flag families with same parents
# check that no-one is 'father' and 'mother'
if( length(intersect(trioListSet[,2],trioListSet[,3])) != 0 ) print( "Someone is mother and father!!" )
# check sex and age
unique(as.character(otherInfo$Inferred.Gender[match(trioListSet[,2],otherInfo$PIID)]))
unique(as.character(otherInfo$Inferred.Gender[match(trioListSet[,3],otherInfo$PIID)]))

max( otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,1],otherInfo$PIID)] - otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,2],otherInfo$PIID)] )
max( otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,1],otherInfo$PIID)] - otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,3],otherInfo$PIID)] )

famID = famsNo3$membership[trioListSet[,1]]

sapply(names(which(table(famID)==2)),function(x) length(unique(trioListSet[famID==x,2:3])) )

print( paste0( nrow(trioListSet)," trios - any family group." ))
print( paste0( length(unique(famID))," unique trio family groups" ))


write.table(trioListSet,file=paste0("forPhasing/",OutputFile,"-all-trios.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)

write.table(weirdUnclean,file=paste0(OutputFile,"-all-trios-odd.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)


###### Investigate weirdunclean
out = c(0,0,0,0)
for (i in weirdUnclean) {
    print(i)
    ourPairs = which(((pc2kins$ID1==i)|(pc2kins$ID2==i))&(pc2kins$class=="parent/child"))
    pids = c(pc2kins$ID1[ourPairs],pc2kins$ID2[ourPairs])
    pids = pids[pids!=i]
    ageChild = otherInfo$Age.when.attended.assessment.centre[otherInfo$PIID==i]
    ageParents = otherInfo$Age.when.attended.assessment.centre[match(pids,otherInfo$PIID)]
    sexParents = as.character(otherInfo$Inferred.Gender[match(pids,otherInfo$PIID)]    )
    ageDiff = ageParents - ageChild
    print(pids)
    print(sexParents)
    print(ageParents)
    print(ageDiff)
    parents = pids[ageDiff>10]
                                        # are any of the 'parents' actually related or duplicates?
    relParents = pc2kins[(pc2kins$ID1%in%parents)&(pc2kins$ID2%in%parents),]
    print(relParents)
    if(length(parents)==3){
        moms = pids[sexParents=="F"]
        dads = pids[sexParents=="M"]
        if(length(dads)>1) print("WARNING: actually this person has twin dads, not mums!")
        out = rbind(out,c(i,moms,dads))
    }
}

out = out[-1,]; colnames(out) = c("Child","Mother1","Mother2","Father")

write.table(out,file=paste0("forPhasing/",OutputFile,"-all-trios-twinMothers.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)



##########################################
###### IDENTIFY PAIRS IN ANY TRIO

trioInds = unique(c(trioListSet[,1],trioListSet[,2],trioListSet[,3]))
inTrios = apply(kin,1,function(x){
                                        #    print(x)
    to = x[[1]];  from = x[[2]]
    if( (to%in%trioInds) & (from%in%trioInds) ) {
        testTo = c(which(trioListSet[,1]==to),which(trioListSet[,2]==to),which(trioListSet[,3]==to)) # all the trios we see ind1
        testFrom = c(which(trioListSet[,1]==from),which(trioListSet[,2]==from),which(trioListSet[,3]==from)) # all the trios we see ind2
        
        if(length(intersect(testTo,testFrom))>0) return(1) else return(0)
        
    } else {
        return(0)
    }        
})


##########################################
###### Families with 3 or more generations?

# child, parent, grandparent sets
# they will look like 'trios' but will have been filtered out because of age differences
toCheck3gens = trioList[which(is.na(trioList[,2])),]
toCheck3gensInds = unique(toCheck3gens[,1])
myAges = otherInfo$Age.when.attended.assessment.centre
names(myAges) = otherInfo$PIID
    
gens3List = sapply(toCheck3gensInds,function(ind){

    rels = kin[((kin$ID1==ind)|(kin$ID2==ind))&(kin$class=="parent/child"),c("ID1","ID2")]
    relsInds = unique(c(rels$ID1,rels$ID2)); relsInds=relsInds[relsInds!=ind]
    ages = myAges[c(ind,relsInds)]
    ageDiffs = ages[1] - ages[2:length(ages)]

    if((sum(ageDiffs>0)>0)&(sum(ageDiffs<0)>0)) {
        print(rels)        
        print(ageDiffs)
        parent = ind
        child = relsInds[ageDiffs>0]
        grandparent = relsInds[ageDiffs<0]
        out = c(child,parent,grandparent)
        print(out)
    } else {
        out = ages
        names(out)=c(ind,relsInds)
                     
    }
    return(out)
},simplify=TRUE)

# ====> There are none!


##########################################
###### FIND QUARTETS

quarts = trioListSet[duplicated(trioListSet[,2]),]
#apply(quarts,1,function(x) {
#    qs = trioListSet[(trioListSet[,2]==x[2])&(trioListSet[,3]==x[3]),]
#    print(dim(qs))
#})




##########################################
###### KINSHIP TABLE WITH INDICATOR FIELDS
# Clean trio (accounting for 3rd degree)
# Clean trio (ignoring for 3rd degree)
# Any trio (all types)
# Family size
# 3 generations

kinTrios = kin
kinTrios$in.trio = inTrios
kinTrios$in.clean.trio = kinTrios$inTrio
kinTrios$in.clean.trio.ignore3degree = 0
clean.trio.no3 = paste0(kinshipNo3$ID1[kinshipNo3$inTrio==1],"---",kinshipNo3$ID2[kinshipNo3$inTrio==1])
kinTrios$in.clean.trio.ignore3degree[paste0(kinTrios$ID1,"---",kinTrios$ID2)%in%clean.trio.no3]=1
kinTrios$family.id = fams$membership[kinTrios$ID1]
kinTrios$family.size = fams$csize[fams$membership[kinTrios$ID1]]

ig=kinTrios$class!="3rd degree"
kinTrios$family.id.ignore3degree = NA
kinTrios$family.size.ignore3degree = NA
kinTrios$family.id.ignore3degree[ig] = famsNo3$membership[kinTrios$ID1[ig]]
kinTrios$family.size.ignore3degree[ig] = famsNo3$csize[famsNo3$membership[kinTrios$ID1[ig]]]


toPrint = kinTrios[,c("ID1","ID2","IBS0","Kinship","class",
                   "family.id",
                   "family.size",
                   "family.id.ignore3degree",
                   "family.size.ignore3degree",
                   "in.clean.trio",
                   "in.clean.trio.ignore3degree","in.trio")]

write.table(toPrint,file=paste0("forPhasing/",OutputFile,"-kinship.txt"),quote=FALSE,row.names=FALSE,sep="\t")




# plot the unclean trios
toPlot = trioListSet[!trioListSet[,1]%in%trioIDs[,1],]

png(paste0("plots/",OutputFile,"-unclean-trios-graph-%02d.png"),width=1500,height=1500,res=150)
for(i in seq(1,ceiling(nrow(toPlot)),by=40) ){
    inds = toPlot[i:min((i+39),nrow(toPlot)),1]
    inds = c(kinshipNo3$ID1[(kinshipNo3$ID1%in%inds)|(kinshipNo3$ID2%in%inds)],kinshipNo3$ID2[(kinshipNo3$ID1%in%inds)|(kinshipNo3$ID2%in%inds)])
    subGraph = igraph::induced_subgraph(networkNo3,inds,impl="copy_and_delete")
    plot.igraph(subGraph,main=paste0("Unlcean trios and non 3rd-degree relatives"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle")
}
dev.off()

png(paste0("plots/",OutputFile,"-unclean-trios-and-relatives-graph-%02d.png"),width=1500,height=1500,res=150)
for(i in seq(1,ceiling(nrow(toPlot)),by=40) ){
    inds = toPlot[i:min((i+39),nrow(toPlot)),1]
    inds = c(kin$ID1[(kin$ID1%in%inds)|(kin$ID2%in%inds)],kin$ID2[(kin$ID1%in%inds)|(kin$ID2%in%inds)])
    subGraph = igraph::induced_subgraph(network,inds,impl="copy_and_delete")
    plot.igraph(subGraph,main=paste0("Unlcean trios and relatives"),vertex.frame.color=NA,vertex.label=NA,vertex.size=2,vertex.shape="circle")
}
dev.off()


# plot the odd sets
