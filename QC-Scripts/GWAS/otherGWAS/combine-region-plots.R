#########
# Combine regions plots into one pdf file from two different runs of plot-GWAS-imputation-comparison.R

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)


#args <- commandArgs(trailingOnly=TRUE)

args = c("Standing.height-BOLT-LMM-v16.out.chr2.maf0.miss1.pruned-raw.v15-lreg-forASHG2","Standing.height-BOLT-LMM-v16.out.chr2.maf0.miss1.pruned-rawQChits-v19-lreg-forASHG2-region") # with v15


prefix1 = args[1]  # with v15
prefix2 = args[2]  # with v19

plot1 = list.files(path="plots/Region.plots/",pattern=paste0(prefix1,".*.png"))

regions = sapply(plot1,function(x) str_split(x,"region")[[1]][2])
pos = as.numeric(sapply(regions,function(x) str_split(x,"-")[[1]][2]))

regions = regions[order(pos)]
plot1 = plot1[order(pos)]

plot2 = sapply(regions,function(r){
    paste0(prefix2,r)
})

filename="plots/Region.plots/Standing.height-BOLT-LMM-v16.out.chr2.maf0.miss1.pruned-raw.v15-WITH-rawQChits.v19.pdf"

pdf(filename,width=41,height=12)
par(mar=c(0,0,0,0))
for(r in 1:length(regions)){
    print(r)    
    img1 <- readPNG(paste0("plots/Region.plots/",plot1[r])) # v15
    img2 <- readPNG(paste0("plots/Region.plots/",plot2[r])) # v19
                                        # plot each plot on an adjacent page
    plot.new()       
    rasterImage(img1,0,0,1,1)
    plot.new()
    rasterImage(img2,0,0,1,1)
}

dev.off()



#quit()




########## chrom X ( only v19s )

#args = c("Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.miss1.pruned-rawQChits.v19s-lreg-forASHG2-region") # Females only in LD computation
#args = c("Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.miss1.pruned-rawQChits.v19s-MandF-lreg-forASHG2-region")  # all samples in computation
args = c("Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.miss1.pruned-rawQChits.v19s-lreg-forASHG2-region")  # all samples in computation


prefix1 = args[1]  # with v15
#prefix2 = args[2]  # with v19


plot1 = list.files(path="plots/Region.plots/",pattern=paste0(prefix1,".*.png"))

regions = sapply(plot1,function(x) str_split(x,"region")[[1]][2])
pos = as.numeric(sapply(regions,function(x) str_split(x,"-")[[1]][2]))

regions = regions[order(pos)]
plot1 = plot1[order(pos)]
# 46 regions in v19s

print(length(regions))

#plot2 = sapply(regions,function(r){
#    paste0(prefix2,r)
#    })

filename="plots/Region.plots/Standing.height-BOLT-LMM-v16s-23.out.chr23.maf0.miss1.pruned-raw.WITH-rawQChits.v19s.pdf"

pdf(filename,width=41,height=12)
par(mar=c(0,0,0,0))
for(r in 1:length(regions)){
    print(r)    
    img1 <- readPNG(paste0("plots/Region.plots/",plot1[r]))
#    img2 <- readPNG(paste0("plots/Region.plots/",plot2[r]))
                                        # plot each plot on an adjacent page
    plot.new()       
    rasterImage(img1,0,0,1,1)
 #   plot.new()
  #  rasterImage(img2,0,0,1,1)
}

dev.off()



quit()




########## chrom PAR ( only v19s )

args = c("Standing.height-BOLT-LMM-v16s-25.out.chr25.maf0.miss1.pruned-rawQChits.v19s-lreg-forASHG2-region") # with v15


prefix1 = args[1]  # with v15
#prefix2 = args[2]  # with v19


plot1 = list.files(path="plots/Region.plots/",pattern=paste0(prefix1,".*.png"))

regions = sapply(plot1,function(x) str_split(x,"region")[[1]][2])
pos = as.numeric(sapply(regions,function(x) str_split(x,"-")[[1]][2]))

regions = regions[order(pos)]
plot1 = plot1[order(pos)]


print(length(regions))

#plot2 = sapply(regions,function(r){
#    paste0(prefix2,r)
#    })

filename="plots/Region.plots/Standing.height-BOLT-LMM-v16s-25.out.chr25.maf0.miss1.pruned-raw.WITH-rawQChits.v19s.pdf"

pdf(filename,width=41,height=12)
par(mar=c(0,0,0,0))
for(r in 1:length(regions)){
    print(r)    
    img1 <- readPNG(paste0("plots/Region.plots/",plot1[r])) # v15
#    img2 <- readPNG(paste0("plots/Region.plots/",plot2[r])) # v19
                                        # plot each plot on an adjacent page
    plot.new()       
    rasterImage(img1,0,0,1,1)
 #   plot.new()
  #  rasterImage(img2,0,0,1,1)
}

dev.off()



quit()

