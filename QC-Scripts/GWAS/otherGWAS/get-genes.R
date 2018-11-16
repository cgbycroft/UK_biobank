# get latest GWAS catalogue data in the specified build. It will download it again if there is an updated version. For the time-stamp of the file, see:
# http://hgdownload.cse.ucsc.edu/goldenPath/<build>/database


args = commandArgs(TRUE)
print(args)

build = args[1]
dataFileDir = args[2]

system(paste0('mkdir ',dataFileDir))
setwd(dataFileDir)

fileName = paste0(dataFileDir,"/",build)
system(paste0('mkdir ',fileName))
setwd(fileName)


url=paste0("http://hgdownload.cse.ucsc.edu/goldenPath/",build,"/database/refGene.txt.gz")

# only download if the new version has been updated
if("refGene.txt.gz"%in%list.files()) system(paste0("wget -N ",url)) else system(paste0("wget -S ",url))

print("File saved: ")
print(fileName)
