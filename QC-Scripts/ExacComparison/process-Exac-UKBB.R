#########
# Script to process AF information from Exac and compare to UKBiobank
#########
library(hexbin)

GenotypesForReleaseFile="/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr"   

args = commandArg(trailingOnly=TRUE)

#args=c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R", "/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-exac","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Exac/ExAC.r1.sites.vep.UKBB-SNPs.aCounts_NFE.txt","-samples","b1__b11-b001__b095-sampleTable_v4")

print(args)
h = args[-c(which(args%in%c("-exac","-ukbb","-samples","-out")),1+which(args%in%c("-exac","-ukbb","-samples","-out")))]
for(helperScript in h){
    source(helperScript)
}


if("-ukbb"%in%args) snpFrequencyFiles = args[which(args=="-ukbb")+1]

releaseSampleQCPrefix = args[which(args=="-samples")+1]

if("-out"%in%args) OutputFile = args[which(args=="-out")+1] else OutputFile = paste0("Exac-",gsub("V2_QCed.export.|.chr","",basename(snpFrequencyFiles)),"-comparison")

exac = args[which(args=="-exac")+1]

print(paste0("using UKBB frequency files from: ",snpFrequencyFiles)," if they exist." )
print(paste0("using UKBB genotype files from: ",GenotypesForReleaseFile) )
                                        # this is defined in auxFunctions.R, but overridden by -ukbb if it exists
print(paste0("using UKBB sample table: ../../data/ForRelease/",releaseSampleQCPrefix,"_allColumns.RData"))


######## Get sample ethnicity and PC information
# Just use the set that were used for SNP-QC, as they're the 'CEU' cluster in 1000G.

load(paste0( "../../data/ForRelease/",releaseSampleQCPrefix,"_allColumns.RData"),verbose=TRUE)

# read in sample lists for SNP-QC

qcsamples = unlist( sapply(all.batches(),function(b) read.table( paste0( baseDataDir,"/",b,"/pca-1000G/",b,"-autosome-HapMap3.pcs.inliers"),stringsAsFactors=FALSE,header=FALSE)[,1] ) )

qcsamples=unique(qcsamples)

print(paste0("found " ,length(qcsamples)," samples used for snp qc."))

# exclude Finnish samples (using place of birth)
outTable$Country.of.birth = get.place.of.birth(outTable)

fins = outTable$PIID[outTable$Country.of.birth=="Finland"]

ukbbEur = qcsamples[!qcsamples%in%fins]

print( paste0("Excluding ",length(intersect(ukbbEur,fins))," Finnish samples used in snp qc"))
print( paste0(length(ukbbEur)," samples to use for UKBB European non-Finnish"))


# Actually NO finnish samples in snp qc. What countries are there?
print( sort(table( outTable$Country.of.birth[outTable$PIID%in%ukbbEur] )) )

write.table(cbind(ukbbEur,ukbbEur),file="ukb-autosome-HapMap3.pcs.inliers.no.fins.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

# Plot the ukbiobank only PCA results for these individuals

PCs = outTable[,grepl("PC",colnames(outTable))]

eth = outTable$Ethnic.background
pop = ethnicity2pop(eth)
Colors = ethnicity2col[pop]
Chars = ethnicity2char[pop]
# fade out samples we are not using (or make them gray?)
#Colors[!outTable$PIID%in%ukbbEur]= add.alpha(Colors[!outTable$PIID%in%ukbbEur],0.3)
Colors[!outTable$PIID%in%ukbbEur]= "gray"

Order = order.by.number.occurrences(Colors)


png(paste0('plots/',OutputFile,'-pca-UKB-NFE-%02d.png'),width=1500,height=1500,res=150)

for (pc in 1:4){
    
    if(pc%%2==0) next

    print(paste0('Plotting PCs',pc,' and ',pc + 1))

    x = PCs[[paste0("PC",pc)]]
    y = PCs[[paste0("PC",pc+1)]]

    plot(x[Order],y[Order],xlab=paste0('PC ',pc),ylab=paste0('PC ',pc+1),
         col=Colors[Order],pch=Chars[Order],axes=FALSE)    
    axis(1,cex.lab=1.3)
    axis(2,cex.lab=1.3)

}
dev.off()




######## Compute allele frequencies in UKBB European subset.

currentFiles = list.files(pattern=paste0(basename(snpFrequencyFiles),"-NFE.*.frqx"),path=dirname(exac))

for (i in c(1:22,"X","XY","Y","MT")){
    myFile = paste0(basename(snpFrequencyFiles),"-NFE",i,".freqx")
    
    #if( myFile %in% currentFiles )  next
    
    forSystem = paste0(plink," --bfile ",GenotypesForReleaseFile,i," --keep ukb-autosome-HapMap3.pcs.inliers.no.fins.txt --freqx --out ",dirname(exac),"/",basename(snpFrequencyFiles),i,"-NFE")

    system(forSystem)
}


####### Read in snps annotations for ukB
BB.ps2snp = rbind(ukbiobank.ps2snp("autosome"),ukbiobank.ps2snp("sexchrom"))
BL.ps2snp = rbind(ukbileve.ps2snp("autosome"),ukbileve.ps2snp("sexchrom"))
ps2snpBoth = unique(BL.ps2snp[BL.ps2snp$AffySNPID%in%BB.ps2snp$AffySNPID,])
ps2snpEither = unique(rbind(BL.ps2snp,BB.ps2snp))

####### Read in the UKBB frequencies

ukbbMaf = read.genotyped.maf(paste0(dirname(exac),"/",basename(snpFrequencyFiles),".*-NFE"))
ukbbBim = rbind_all(read.genotyped.snps("genome",genotypedSNPFile=GenotypesForRelease)) # read bims from GenotypesForRelease
ukbbMaf2 = left_join(ukbbMaf,ukbbBim[,c("SNP","BP","SNP2")],by=c("SNP"="SNP"))
ukbbMaf2$obs_A1 =  2*ukbbMaf2$C.HOM.A1. + ukbbMaf2$C.HET. + ukbbMaf2$C.HAP.A1
ukbbMaf2$obs_A2 =  2*ukbbMaf2$C.HOM.A2. + ukbbMaf2$C.HET. + ukbbMaf2$C.HAP.A2
ukbbMaf2$obs_chroms = ukbbMaf2$obs_A1 + ukbbMaf2$obs_A2

# how many observations (chromosomes) should there be altogether? 
ukbbMaf2$ukb_total_chroms = ukbbMaf2$obs_chroms + 2*ukbbMaf2$C.MISSING.

# non autosomes different. Only count haploid chromosomes once.

nMales = sum(outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="M")
nFemales = sum(outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="F")
ubl=outTable$Batch[outTable$PIID%in%ukbbEur]%in%ukbileve.batches()
ubb=outTable$Batch[outTable$PIID%in%ukbbEur]%in%ukbiobank.batches()

nM_ukbl = sum(( outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="M" )&( ubl ) )
nM_ukbb = sum(( outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="M" )&( ubb ) )
nF_ukbl = sum(( outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="F" )&( ubl ) )
nF_ukbb = sum(( outTable$Inferred.Gender[outTable$PIID%in%ukbbEur]=="F" )&( ubb ) )

ukbbSNPs = (ukbbMaf2$SNP%in%BB.ps2snp$dbSNPRSID)&(!ukbbMaf2$SNP%in%BL.ps2snp$dbSNPRSID)
ukblSNPs = (ukbbMaf2$SNP%in%BL.ps2snp$dbSNPRSID)&(!ukbbMaf2$SNP%in%BB.ps2snp$dbSNPRSID)
ukbothSNPs = (ukbbMaf2$SNP%in%BL.ps2snp$dbSNPRSID)&(ukbbMaf2$SNP%in%BB.ps2snp$dbSNPRSID)

table(ukbbSNPs+ukblSNPs+ukbothSNPs)


# both ukbileve and ukbiobank snps
ukbbMaf2$ukb_total_chroms[(ukbothSNPs)] = 2*(nMales + nFemales) # autosomes
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==23)&(ukbothSNPs)] = nMales + 2*nFemales # X nonPAR
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==24)&(ukbothSNPs)] = nMales # Y
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==26)&(ukbothSNPs)] = nMales + nFemales # MT

# only ukbiobank snps
ukbbMaf2$ukb_total_chroms[(ukbbSNPs)] = 2*(nM_ukbb + nF_ukbb) # autosomes
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==23)&(ukbbSNPs)] = nM_ukbb + 2*nF_ukbb # X nonPAR
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==24)&(ukbbSNPs)] = nM_ukbb # Y
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==26)&(ukbbSNPs)] = nM_ukbb + nF_ukbb # MT

# only ukbileve snps
ukbbMaf2$ukb_total_chroms[(ukblSNPs)] = 2*(nM_ukbl + nF_ukbl) # autosomes
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==23)&(ukblSNPs)] = nM_ukbl + 2*nF_ukbl # X nonPAR
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==24)&(ukblSNPs)] = nM_ukbl # Y
ukbbMaf2$ukb_total_chroms[(ukbbMaf2$CHR==26)&(ukblSNPs)] = nM_ukbl + nF_ukbl # MT



####### Read in the Exac frequencies

exacMaf = read.table(exac,header=FALSE,stringsAsFactors=FALSE)

exacMaf$SNP2 = altSNPID(exacMaf$V1,exacMaf$V2,exacMaf$V3,exacMaf$V4)
exacMaf$is.multi = 1*grepl(",",exacMaf$V4)
colnames(exacMaf)[5:6]= c("AC_NFE","AN_NFE")
                                        # parse the multi-allelic ones
alleles = str_split(exacMaf$V4,",") # the alleles
counts1 = str_split(exacMaf$AC_NFE,",") # the counts
counts2 = str_split(exacMaf$AN_NFE,",") # the counts

lengths=sapply(alleles,length)

ALELLES=sapply(1:max(lengths),function(l){
    sapply(1:length(alleles),function(x) if(lengths[x]>=l) alleles[[x]][l] else NA )
})
COUNTS1=sapply(1:max(lengths),function(l){
    sapply(1:length(counts1),function(x) if(lengths[x]>=l) counts1[[x]][l] else NA )
})



####### Join with exac data

MATCHES = sapply(1:max(lengths), function(c){
    s = altSNPID(exacMaf$V1,exacMaf$V2,exacMaf$V3,ALELLES[,c])
    ms = match(s,ukbbMaf2$SNP2)    
})

MATCHES_alt= sapply(1:max(lengths), function(c){
    s = altSNPID(exacMaf$V1,exacMaf$V2,ALELLES[,c],exacMaf$V3)
    ms = match(s,ukbbMaf2$SNP2)    
})


# this just makes sure that there is only one matching index per snp
test = unique( rowSums(is.na(MATCHES)) )
if( sum( !test%in%(ncol(MATCHES)-c(0,1)) )>0 ) print("ERROR!!! there is more than one matching index per snp!")

matches1 = rowSums(MATCHES,na.rm=TRUE); matches1[matches1==0]=NA

                                        # what have we missed?
test = unique( rowSums(is.na(MATCHES_alt)) )
if( sum( !test%in%(ncol(MATCHES_alt)-c(0,1)) )>0 ) print("ERROR!!! there is more than one matching index per snp!")

matches2 = rowSums(MATCHES_alt,na.rm=TRUE); matches2[matches2==0]=NA

matchesAll = matches1
matchesAll[is.na(matches1)] = matches2[is.na(matches1)]

chrPos = paste0(exacMaf$V1,exacMaf$V2)
chrPosMatches = match(chrPos,paste0(ukbbMaf2$CHR,ukbbMaf2$BP))
    
print( paste0( sum(!is.na(matchesAll))," markers found in both ExAC and UKBiobank"))


    
###### Do the join

exacUKB = cbind(exacMaf,ukbbMaf2[matchesAll,])
# get the correct alleles and counts from exac

Alts = sapply(1:nrow(exacUKB),function(i) {
    t = ALELLES[i,!is.na(MATCHES[i,]) ]
    if(length(t)==0){
        t = ALELLES[i,!is.na(MATCHES_alt[i,]) ]
    }
    if(length(t)==0){
        return(NA)
    } else {
        return(t)
    }
})

AC_NFE = sapply(1:nrow(exacUKB),function(i) {
    t = COUNTS1[i,!is.na(MATCHES[i,]) ]
    if(length(t)==0){
        t = COUNTS2[i,!is.na(MATCHES_alt[i,]) ]
    }
    if(length(t)==0){
        return(NA)
    } else {
        return(t)
    }
})



exacUKB$Alt = Alts
exacUKB$Ref = exacUKB$V3
exacUKB$exac_alt_count =  as.numeric(AC_NFE)
exacUKB$exac_n_count = exacUKB$AN_NFE
exacUKB$exac_ref_count = exacUKB$exac_n_count - exacUKB$exac_alt_count

exacUKB$freq_exac = exacUKB$exac_alt_count/exacUKB$exac_n_count


# align frequencies to the UKbiobank reference alleles
flipped =  (exacUKB$Ref!=exacUKB$A2)&(!is.na(exacUKB$A2))
sum((exacUKB$Ref[flipped]!=exacUKB$A1[flipped])&(!is.na(exacUKB$A2[flipped]))) # this should be zero, because all non-flipped alleles should have the same A1 alleles.

exacUKB$freq_exac_align = exacUKB$freq_exac
exacUKB$freq_exac_align[flipped] = 1-exacUKB$freq_exac[flipped]

exacUKB$MAF_exac_align = exacUKB$freq_exac_align
exacUKB$MAF_exac_align[(exacUKB$freq_exac_align>0.5)&(!is.na(exacUKB$freq_exac_align))] = 1 - exacUKB$freq_exac_align[(exacUKB$freq_exac_align>0.5)&(!is.na(exacUKB$freq_exac_align))]

max(exacUKB$MAF_exac_align - exacUKB$MAF,na.rm=TRUE)

o = exacUKB[order(exacUKB$MAF_exac_align - exacUKB$MAF),]


# Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data rs58571639 cluster_plots -b UKBiLEVEAX_b1,UKBiLEVEAX_b11,Batch_b001,Batch_b020,Batch_b060 -released &

###################
# compute a z-score for the differences in proportions of alleles (assume chromosomes are random within an individual)

zscore <- function(c1,c2,n1,n2){
                                        # for the H0 that p1!=p2
    p1 = c1/n1
    p2 = c2/n2
    p = (c1 + c2)/(n1 + n2)
    z = (p1 - p2)/(sqrt( p*(1-p) * ( (1/n1) + (1/n2))  ) )
    return(z)
}


# Desi's Fisher exact test
#pvalues = .Call("TestGenotypeFreqs",, 
#    as.integer(mcounts1), as.integer(fcounts2), as.integer(mcounts2))

# aligned counts, and treat like haploid calls in desi's algorithm

ukbbCounts = cbind(exacUKB$obs_A1,0,exacUKB$obs_A2,0)
exacCounts = cbind(exacUKB$exac_alt_count,0,exacUKB$exac_ref_count,0)
exacCounts[flipped,] = exacCounts[flipped,c(3,2,1,4)]

data = cbind(exacUKB$CHR,ukbbCounts,0,0,0,0,0,0,0,0,exacCounts)

fisherP = t(apply(data,1,CompareGenotypes))[,1]

zscores = apply(data,1,function(x) zscore(x[2],x[14],x[2]+x[4],x[14]+x[16]) )
zscoreP = 2*pnorm(-abs(zscores))


# Exclude anything that doesn't match to UKBB; and has actual non-sero read counts (AN_NFE) in ExAC
matched=( !is.na(exacUKB$MAF) )&( !is.na(exacUKB$MAF_exac_align ) )

print( paste0( sum(is.na(exacUKB$MAF_exac_align)), " markers have 0 read counts in ExAC. I.e AN_NFE=0, this leaving ",sum(matched)," matched markers for analysis." ) )

exacUKB = exacUKB[matched,]
zscoreP = zscoreP[matched]
fisherP = fisherP[matched]

###################
# save the output

save(exacUKB,fisherP,zscoreP,file=paste0(OutputFile,".RData"))




###################
# some plotting


apoeGene = c(45409006,45412652)
inapoe = (exacUKB$CHR==19)&(exacUKB$BP >= apoeGene[1])&(exacUKB$BP <= apoeGene[2])

Ymax=100
fisherP2 = fisherP; fisherP2[-log10(fisherP)>=Ymax] = 10^-Ymax

nsnpsFound = sum(!is.na(exacUKB$MAF))
mafDif = exacUKB$MAF - exacUKB$MAF_exac_align # differences between MAF in UKBiobank (with minor allele defined by biobank), and the frequency in exac.
relmafDif = mafDif/exacUKB$MAF

chisq1 = qchisq(fisherP,1,lower.tail=FALSE)
lambda1 = median(chisq1,na.rm=TRUE)/qchisq(0.5,1)


###### Some snp categories

# monomorphic in both
monoBoth = (exacUKB$MAF==0)&(exacUKB$MAF_exac_align==0)

# polymorphic in both
polyBoth = (exacUKB$MAF>0)&(exacUKB$MAF_exac_align>0)

# rare
rare = (exacUKB$MAF_exac_align<0.01)&(exacUKB$MAF<0.01)
vrare = (exacUKB$MAF_exac_align<0.001)&(exacUKB$MAF<0.001)
vrare2 = (exacUKB$MAF_exac_align<0.005)&(exacUKB$MAF<0.005)

# high call rates in both exac in ukbb
missThresh = 0.9
N_exac=33370*2 # This is from the NFE population in ExAC: http://exac.broadinstitute.org/faq
incl = ( exacUKB$exac_n_count > N_exac*missThresh )&( (exacUKB$obs_chroms/exacUKB$ukb_total_chroms) > missThresh )

# low fisher p-value
pthresh = 10^-12  # same as our other tests
lowP = fisherP < pthresh
    
######
###### Hex bin colors
my.hex.colors = colorRampPalette(colors=c("darkgray","blue2"))
######


# plot the differences
png(paste0('plots/',OutputFile,'-ExacComparison-histdiffs.png'),bg="transparent",width=1500,height=1500,res=150)
hist(mafDif,xlab="MAF.UKBiobank - MAF.ExAC",main=paste0(nsnpsFound," markers in ExAC matched to UK Biobank"),breaks=1000)
dev.off()

png(paste0('plots/',OutputFile,'-ExacComparison-qq.png'),bg="transparent",width=1500,height=1500,res=150)
qqman::qq(fisherP)
dev.off()


    
colors = rep("blue",length(mafDif))
colors[lowP] = "red"
#borders=colors
borders="transparent"
#borders[inapoe] = "black"

Order = order.by.number.occurrences(colors)

png(paste0('plots/',OutputFile,'-ExacComparison-scatter.png'),bg="transparent",width=1500,height=1500,res=150)
plot(exacUKB$MAF[Order],exacUKB$MAF_exac_align[Order],main=paste0(nsnpsFound," markers in exac matched to UK Biobank"),xlab="MAF.UKBiobank",ylab="MAF.ExAC",bg=add.alpha(colors[Order],0.2),col=borders[Order],pch=21)
dev.off()

hex1 = hexbin(exacUKB$MAF[Order],exacUKB$MAF_exac_align[Order],xbins=100)

png(paste0('plots/',OutputFile,'-ExacComparison-hexbin.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex1,main=paste0(nsnpsFound," markers in exac matched to UK Biobank"),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,10^seq(0:4)/sum( 10^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()

hex2 = hexbin(exacUKB$MAF[rare],exacUKB$MAF_exac_align[rare],xbins=100)

png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-zoom.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex2,main=paste0(sum(rare)," markers in exac matched to UK Biobank"),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,10^seq(0:4)/sum( 10^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()


# Filter exac + UKBiobank by some quality metric...
# First missing data. Compute Ukbiobank missing data by the two arrays separately
Order = order.by.number.occurrences(colors[incl])

png(paste0('plots/',OutputFile,'-ExacComparison-scatter-exacMiss',missThresh,'.png'),bg="transparent",width=1500,height=1500,res=150)
plot(exacUKB$MAF[incl][Order],exacUKB$MAF_exac_align[incl][Order],main=paste0(sum(incl,na.rm=TRUE)," markers in exac matched to UK Biobank\nand with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",bg=add.alpha(colors[incl][Order],0.2),col=borders[incl][Order],pch=21)
dev.off()

Order = order.by.number.occurrences(colors[(rare)&(incl)])
png(paste0('plots/',OutputFile,'-ExacComparison-scatter-exacMiss',missThresh,'-zoom.png'),bg="transparent",width=1500,height=1500,res=150)
plot(exacUKB$MAF[(rare)&(incl)][Order],exacUKB$MAF_exac_align[(rare)&(incl)][Order],main=paste0(sum((rare)&(incl),na.rm=TRUE)," markers in exac matched to UK Biobank\nand with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",bg=add.alpha(colors[(rare)&(incl)][Order],0.2),col=borders[(rare)&(incl)][Order],pch=21)
dev.off()

hex3 = hexbin(exacUKB$MAF[(rare)&(incl)],exacUKB$MAF_exac_align[(rare)&(incl)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'-zoom.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex3,main=paste0(sum((rare)&(incl),na.rm=TRUE)," rare markers"),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,2^seq(0:10)/sum( 2^seq(0:10)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()


hex3c = hexbin(exacUKB$MAF[(vrare2)&(incl)],exacUKB$MAF_exac_align[(vrare2)&(incl)],xbins=80)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'-zoom2.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex3c,main=paste0(sum((vrare2)&(incl),na.rm=TRUE)," markers in exac matched to UK Biobank\nand with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,2^seq(0:6)/sum( 2^seq(0:6)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()

hex3a = hexbin(exacUKB$MAF[(incl)],exacUKB$MAF_exac_align[(incl)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex3a,main=paste0(sum((incl),na.rm=TRUE)," markers in exac matched to UK Biobank\nand with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,5^seq(0:6)/sum( 5^seq(0:6)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()




# Manhattan plot!
png(paste0('plots/',OutputFile,'-ExacComparison-manhattan-exacMiss',missThresh,'.png'),width=61,height=12,units="in",res=150,bg="transparent")
myManhattan(cbind(fisherP,exacUKB),p="fisherP",ymax=Ymax,suggestiveline = FALSE,xpd=NA)
dev.off()


                                        # Fraction of SNPs in each bin, by UKBiobank MAF bin that have p-value within the threshold

mafBins = c(0,1/10000,1/1000,1/100,5/100,Inf)
myBinCols = c("yellow","red","purple","blue","green")
names(mafBins) = paste0("[",paste(mafBins[1:(length(mafBins)-1)],mafBins[2:length(mafBins)],sep=","),")")
names(mafBins)[length(mafBins)-1] = paste0("[",mafBins[length(mafBins)-1])

bins = getMafBins(exacUKB$MAF,mafBins=mafBins)
Bins = factor(bins)
levels(Bins) = names(mafBins)[as.numeric(levels(Bins))]

png(paste0('plots/',OutputFile,'-ExacComparison-pval-boxplot-byMaF.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(-log10(fisherP2)~Bins,ylab="-log(p), Fisher exact test",xlab="MAF in UK Biobank")
dev.off()

png(paste0('plots/',OutputFile,'-ExacComparison-pval-boxplot-byMaF-exacMiss',missThresh,'.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(-log10(fisherP2[incl])~Bins[incl],ylab="-log(p), Fisher exact test",xlab="MAF in UK Biobank")
dev.off()


png(paste0('plots/',OutputFile,'-ExacComparison-pval-hist-byMaF-exacMiss',missThresh,'.png'),bg="transparent",width=1500,height=1500,res=150)
bre=seq(0,Ymax,by=1)
for(i in 1:(length(mafBins)-1)){
    x = -log10(fisherP2[(incl)&(bins==i)])
    if(i==1) hist(x,xlab="-log(p), Fisher exact test",freq=FALSE,breaks=bre,col=add.alpha(myBinCols[i],0.5),border=NA ,xlim=c(5,Ymax),ylim=c(0,0.2)) else hist(x,freq=FALSE,breaks=bre ,add=TRUE,col=add.alpha(myBinCols[i],0.5),border=NA,xlim=c(5,Ymax),ylim=c(0,0.2))
}
dev.off()



##### Exclude monomorphic in both


hex6 = hexbin(exacUKB$MAF[(!monoBoth)],exacUKB$MAF_exac_align[(!monoBoth)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex6,main=paste0(sum((incl)&(!monoBoth),na.rm=TRUE)," markers in exac matched to UK Biobank\n and polymorphic in at least one study"),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,10^seq(0:4)/sum( 10^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()

# NOTE: non-monomorphic means at least one of ukb and exac is non-monomorphic
hex4 = hexbin(exacUKB$MAF[(incl)&(!monoBoth)],exacUKB$MAF_exac_align[(incl)&(!monoBoth)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex4,main=paste0(sum((incl)&(!monoBoth),na.rm=TRUE)," markers in exac matched to UK Biobank\nand polymorphic in at least one study; with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,10^seq(0:4)/sum( 10^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()


hex5 = hexbin(exacUKB$MAF[(rare)&(incl)&(!monoBoth)],exacUKB$MAF_exac_align[(rare)&(incl)&(!monoBoth)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'-polyboth-zoom.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex5,main=paste0(sum((rare)&(incl)&(!monoBoth),na.rm=TRUE)," rare markers in exac matched to UK Biobank\nand polymorphic in at least one study; with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,10^seq(0:4)/sum( 10^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()


hex7 = hexbin(exacUKB$MAF[(vrare)&(rare)&(incl)&(!monoBoth)],exacUKB$MAF_exac_align[(vrare)&(rare)&(incl)&(!monoBoth)],xbins=100)
png(paste0('plots/',OutputFile,'-ExacComparison-hexbin-exacMiss',missThresh,'-polyboth-zoom2.png'),bg="transparent",width=1500,height=1500,res=150)
hb = plot(hex7,main=paste0(sum((vrare)&(rare)&(incl)&(!monoBoth),na.rm=TRUE)," very rare markers in exac matched to UK Biobank\nand polymorphic in at least one study; with call rate > ",missThresh),xlab="MAF.UKBiobank",ylab="MAF.ExAC",colorcut=c(0,5^seq(0:4)/sum( 5^seq(0:4)),1),colramp=my.hex.colors)
hexVP.abline(hb$plot,0,1, col= "red",lty=3)
dev.off()


png(paste0('plots/',OutputFile,'-ExacComparison-pval-boxplot-byMaF-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(-log10(fisherP2[(!monoBoth)])~Bins[(!monoBoth)],ylab="-log(p), Fisher exact test",xlab="MAF in UK Biobank")
dev.off()

png(paste0('plots/',OutputFile,'-ExacComparison-pval-boxplot-byMaF-exacMiss',missThresh,'-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(-log10(fisherP2[(incl)&(!monoBoth)])~Bins[(incl)&(!monoBoth)],ylab="-log(p), Fisher exact test",xlab="MAF in UK Biobank")
dev.off()


png(paste0('plots/',OutputFile,'-ExacComparison-diff-boxplot-byMaF-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(mafDif[(!monoBoth)]~Bins[(!monoBoth)],ylab="MAF.UKBiobank - MAF.ExAC",xlab="MAF in UK Biobank")
dev.off()

png(paste0('plots/',OutputFile,'-ExacComparison-diff-boxplot-byMaF-exacMiss',missThresh,'-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(mafDif[(incl)&(!monoBoth)]~Bins[(incl)&(!monoBoth)],ylab="MAF.UKBiobank - MAF.ExAC",xlab="MAF in UK Biobank")
dev.off()


png(paste0('plots/',OutputFile,'-ExacComparison-reldiff-boxplot-byMaF-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(relmafDif[(!monoBoth)]~Bins[(!monoBoth)],ylab="(MAF.UKBiobank - MAF.ExAC)/MAF.UKBiobank",xlab="MAF in UK Biobank")
dev.off()

png(paste0('plots/',OutputFile,'-ExacComparison-reldiff-boxplot-byMaF-exacMiss',missThresh,'-polyboth.png'),bg="transparent",width=1500,height=1500,res=150)
boxplot(relmafDif[(incl)&(!monoBoth)]~Bins[(incl)&(!monoBoth)],ylab="(MAF.UKBiobank - MAF.ExAC)/MAF.UKBiobank",xlab="MAF in UK Biobank")
dev.off()



###### Cluster plots

# zero in one but not the other
zerosUKbb = (exacUKB$MAF==0)&(exacUKB$MAF_exac_align>0.0001)
zerosExac = (exacUKB$MAF>0.0001)&(exacUKB$MAF_exac_align==0)

# ukbb zero
snps1 = exacUKB$SNP[zerosUKbb&incl&(fisherP < pthresh)]
length(snps1)

#length(snps1)
#[1] 294

# exac zero
snps2 = exacUKB$SNP[zerosExac&incl&(fisherP < pthresh)]
length(snps2)

#length(snps2)
#[1] 59


system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(snps1[1:5],collapse=",")," plots/cluster_plots_mono_ukbb -b UKBiLEVEAX_b1,UKBiLEVEAX_b11,Batch_b001,Batch_b020,Batch_b060 -released -title random5 &") )


system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(snps2[1:5],collapse=",")," plots/cluster_plots_mono_exac -b UKBiLEVEAX_b1,UKBiLEVEAX_b11,Batch_b001,Batch_b020,Batch_b060 -released -title random5 &") )



ext = exacUKB$SNP[zerosUKbb&incl&(fisherP < pthresh)&(exacUKB$MAF_exac_align>0.01)]
system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(ext[1:5],collapse=",")," plots/cluster_plots_mono_ukbb -b UKBiLEVEAX_b1,UKBiLEVEAX_b11,Batch_b001,Batch_b020,Batch_b060 -released -title extreme5 &") )

ext = exacUKB$SNP[zerosExac&incl&(fisherP < pthresh)&(exacUKB$MAF>0.01)]
system( paste0( "Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(ext[1:5],collapse=",")," plots/cluster_plots_mono_exac -b UKBiLEVEAX_b1,UKBiLEVEAX_b11,Batch_b001,Batch_b020,Batch_b060 -released -title extreme5 &") )

