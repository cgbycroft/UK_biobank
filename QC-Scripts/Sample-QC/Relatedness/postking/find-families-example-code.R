library(igraph)


###########
# INPUTS:


KinshipFile = "<Text file of pairwise kinships.  Must contain columns: ID1 and ID2 (the eids of individuals in each pair); Kinship; and IBS0 >"
AgeSexFile = "<Text file contain age and sex of all individuals.  Must contain the individual ID column to match to the kinship table.>"
OutputFile = "<Name of the file you want to write the trios to>"



###########
# Function to define classify related pairs

degrees = c(1/(2^(3/2)),1/(2^(5/2)),1/(2^(7/2)),1/(2^(9/2)))

get.kin.classes <- function(kin,ibs0Threshold=0.0012 ){
    
    ## get indexes of kinship classes
    classes = rep(NA,dim(kin)[1])
    p = kin$Kinship
    q = kin$IBS0
    classes[ p > degrees[1]  ] = "dupe/twins"
    classes[(p <= degrees[1] ) & ( q >= ibs0Threshold ) & (p > degrees[2] )] = "sibs"    
    classes[( p <= degrees[2] ) & ( p > degrees[3] )] = "2nd degree"
    classes[( p  <= degrees[3] ) & (p  > degrees[4] ) ] = "3rd degree"
    classes[(p <= degrees[1] ) & ( q < ibs0Threshold )] = "parent/child"

    if( sum(is.na(classes)) > 0 ) print("Some pairs unable to be classified...")
    return(classes)

}



###########
# Read UK Biobank kinship table & age/sex data

kin = read.table(KinshipFile,header=T,stringsAsFactors=F)

otherInfo = read.table(AgeSexFile,header=T,stringsAsFactors=F)




###########
# get relationship classes

kin$class = get.kin.classes(kin,ibs0Threshold=ibs0Threshold)
kin$class2 = factor(kin$class)

# exclude 3rd-degree pairs
kinNo3 = kin[kin$class!="3rd degree",]
kinNo3$class2 = factor(kinNo3$class)



###########
# Find and tabulate clusters of relatives ignoring 3rd-degree pairs. <-- Trios must contain at least 2 parent/child pairs

networkNo3 <- graph.data.frame(kinNo3, directed=F)
famsNo3 = clusters(networkNo3)

famsNo33 = which(famsNo3$csize > 2)  # index of which families contain 3 or more people.

# Create table of family configurations i.e. numbers of different relationship types in each family cluster.

famsNo33Configs = sapply(1:length(famsNo33),function(i){
        if(i%%100==0) print(i)
        x=famsNo33[i]
        inds=names(famsNo3$membership)[famsNo3$membership==x]
        config = table(kinNo3$class2[(kinNo3$ID1%in%inds)|(kinNo3$ID2%in%inds)]) # The count of pairwise connections.
    })
colnames(famsNo33Configs) = famsNo33



###########
# Find all trios

# The distinctive pattern is: at least 2 parent/child pairs in family, and children in trios will always appear more than twice in the kinship table as a parent/child pair.
# we then check ages and sex to figure out which of these clusters contain a child with two parents (as opposed to a parent with at least two children).

# NOTE:  The individual ID column in the 'otherInfo' table is called here "PIID". You might need to change this depending on your input file format. You'll also need the columns "Age.when.attended.assessment.centre" and "Inferred.Gender".


# Family groups with at least 2 parent/child pairs:
pc2 = famsNo33Configs[,famsNo33Configs["parent/child",] >= 2]


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

# For each parent/child pair, check if the child has two parents. The child is the younger of the two (assuming ages are correct)

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


# Remove those where the child doesn't have two parents or otherwise not actually a trio.
trioListSet = trioList[!is.na(trioList[,2]),]
colnames(trioListSet) = c("Child","Mother","Father")


#######
# flag families with same parents
# check that no-one is 'father' and 'mother'
if( length(intersect(trioListSet[,2],trioListSet[,3])) != 0 ) print( "Someone is mother and father!!" )

# re-check sex and age of parents and children
unique(as.character(otherInfo$Inferred.Gender[match(trioListSet[,2],otherInfo$PIID)]))
unique(as.character(otherInfo$Inferred.Gender[match(trioListSet[,3],otherInfo$PIID)]))

max( otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,1],otherInfo$PIID)] - otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,2],otherInfo$PIID)] )
max( otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,1],otherInfo$PIID)] - otherInfo$Age.when.attended.assessment.centre[match(trioListSet[,3],otherInfo$PIID)] )


#######
# Print counts

famID = famsNo3$membership[trioListSet[,1]]

print( paste0( nrow(trioListSet)," trios - any family group." )) # <--- Any combination of parents/child (quartets will be counted twice here)
print( paste0( length(unique(famID))," unique trio family groups" )) # <--- Family units will only be counted once here.


#######
# Write table where each line is a trio

write.table(trioListSet,file=OutputFile,quote=FALSE,col.names=TRUE,row.names=FALSE)


print("You're done!")
