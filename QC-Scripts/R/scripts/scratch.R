# ld statistics computation using bgen files (NOT FINISHED!!)

bgenix = '/apps/well/bgenix/20160708/bgenix'
edit.bgen='/tmp/for_clare/edit-bgen'  # this is on rescomp1

# to eventually go in auxFunctions.R
compute.ld.around.snp.bgen <- function(variants,bgenFile,bgenSampleFile,chrom,samplesFile,regions=NULL,DF){

    print(head(variants))
    
    if(is.null(bgenSampleFile)) {
        print("WARNING: must specify bgen sample file with argument -sample. Otherwise can't compute LD properly.")
        return(NULL)
    }

    # first subset bgen file using bgenx. make sure to subset based on the set of samples in samplesFile (e.g whatever set was used in the GWAS)
    rNumber = round(runif(1,1,100000),0)

    write.table(unique(variants),file = paste0("topVariants.",rNumber,".tmp"),quote=FALSE,col.names=FALSE,row.names=FALSE)

    
    DFchrom = DF[DF$CHR==chrom,]
    
     positions = sapply(1:length(variants),function(s){

        v = variants[s]
        print(v)
                                        # This won't necessarily match the plink '-ldwindow' argument, so just compute ld within the specified hit region (might contain a heap of SNPs!). If no region argument is supplied, just extract snps within 1MB of the centre SNP.
        
        if(!is.null(regions)){
            reg = regions[s,]
            starts = as.numeric(reg[1])
            ends = as.numeric(reg[2])
                                        # if region is smaller than 1MB either side of SNP, then take that region.

           # if( starts > ( DF$BP[DF$SNP==v] - 20e3 ) ) starts = DF$BP[DF$SNP==v] - 20e3
           # if( ends < ( DF$BP[DF$SNP==v] + 20e3 ) ) ends = DF$BP[DF$SNP==v] + 20e3
            
        } else {
            starts = DF$BP[DF$SNP==v] - 20e3; starts[starts<0]=0 
            ends = DF$BP[DF$SNP==v] + 20e3
        }

                                       
        #######  # TEMP FOR TESTING
        #starts = DF$BP[DF$SNP==v] - 10e3; starts[starts<0]=0 
        #ends = DF$BP[DF$SNP==v] + 10e3
        ########
        print(starts)
        print(ends)

        if(starts<0) starts=0
        if(ends>max(DFchrom$BP)) ends=max(DFchrom$BP)
        
        nSNPs = sum( (DFchrom$BP<=ends)&(DFchrom$BP>=starts) )
        if( nSNPs > 10000 ) print( paste0( nSNPs ," SNPs in chunk. May take a long time to process...") )
        
        chrs = sprintf("%02d",DF$CHR[DF$SNP==v])

        posString = paste0(chrs,":",starts,"-",ends)
        return(posString)
        
#        ldThisSNP = ld.bgen(snp,chrs,starts,ends,rNumber)        
    
    },simplify=FALSE)

    ld = ld.bgen(unlist(positions),rNumber)   

    #system(paste0("rm topVariants.",rNumber,".tmp.bgen"))
    #system(paste0("rm topVariants.",rNumber,".tmp.bed"))
    #system(paste0("rm topVariants.",rNumber,".tmp.bim"))
    #system(paste0("rm topVariants.",rNumber,".tmp.fam"))
    #system(paste0("rm topVariants.",rNumber,".tmp"))
    #system(paste0("rm taggedSNPs.",rNumber,".tmp.ld"))
    #system(paste0("rm taggedSNPs.",rNumber,".tmp.ld.log"))
    
#    ld = rbind_all(positions)
        
    return(ld)
}


#ld.bgen <- function(snp,chrs,starts,ends,rNumber){
ld.bgen <- function(positions,rNumber){
    
#    positions = paste0(chrs,":",starts,"-",ends)
    positions2 =  paste(positions,collapse=" ")

    # record set of positions used in this run
    write.table(positions,file=paste0("chr.",chrom,".genotypes.",rNumber,".tmp.positions.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)

                                        #run bgenx to get subset bgen file.   
    forSystem1 = paste0(bgenix," ",bgenFile," -incl-range ",positions2," > genotypes.",rNumber,".tmp.bgen")
    system(forSystem1)

                                        # The above creates a v1.1 bgen file with sample names in the header. Need to convert to something that plink can read!
    
#    forSystem2 = paste0(qctool," -g genotypes.",rNumber,".tmp.bgen -og genotypes.",rNumber,".tmp2.bgen")

 #   system(forSystem2)

    
                                        # ~1,000 variants in a typical region for imputed data!? Takes about 20secs to chunk up.
#    tempBgen = paste0("genotypes.",rNumber,".tmp2.bgen")
    tempBgen = paste0("genotypes.",rNumber,".tmp.bgen")

                                        # remove the sample names so plink can read it!
    forSystem2 = paste0(edit.bgen," -g ",tempBgen," -remove-sample-identifiers -really" )
    system(forSystem2)

                                        # convert to bed format (this can take a while (5mins or so). Depending on how big, and how many regions you have. But fast once you have bed format.)
    forSystem3 = paste0( plink," --bgen ",tempBgen," --sample ",bgenSampleFile," --list-duplicate-vars suppress-first --keep-allele-order --make-bed --out ",gsub(".bgen","",tempBgen))
    system(forSystem3)
    
    forSystem4 = paste0("awk '{print $4}' ",gsub(".bgen",".dupvar",tempBgen)," | tail -n +2 > ",gsub(".bgen",".dupvar.pos",tempBgen))
    system(forSystem4)
    
                                        # call plink to compute LD
                                        # NOTE: with genotype data we used --ld-window 1000 which means compute LD with only 1000 markers either side of main SNP. SNPs are much denser on imputed data, so will use a kb distance instead. 1000 markers in genotype data ~20000 BP = 20KB.
                                        # test = read.table(/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.bim)
                                        # markers = test[,4]
                                        # diffs1K = sapply(1:(length(markers)-1000), function(i) markers[i+1000]-markers[i])
                                        # mean(diffs1K)
    # we also didn't use a r2 threshold in the genotype data, but given that we only plot colours above 0.1, just use that value to avoid recording small r2 values. Also, allow the window to test any SNP within the regions we are testing.
                                        

    forSystem5 = paste0(plink," --bfile ",gsub(".bgen","",tempBgen)," --keep-allele-order --exclude ",gsub(".bgen",".dupvar.pos",tempBgen)," --keep ",samplesFile," --chr ",chrom," --r2 with-freqs dprime --ld-snp-list topVariants.",rNumber,".tmp --ld-window 100000 --ld-window-kb 1000 --ld-window-r2 0.1 --out taggedSNPs.",rNumber,".tmp")

    system(forSystem5)    
                                        # read back the LD results for this SNP
    
    ldThisSNP = read.table(paste0("taggedSNPs.",rNumber,".tmp.ld"),header=TRUE)
        
    return(ldThisSNP)    
    
}
