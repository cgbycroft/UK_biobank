

args = commandArgs(trailingOnly = TRUE)
VarFile = args[1]


library(aberrant)
lambda = 20


X = read.table(VarFile,header=FALSE,stringsAsFactors=FALSE)
ID = X[,1]
x1 = X[,2]
x2 = X[,3]


set.seed(123456)
rslt = aberrant(data.frame(x1=x1,x2=x2),lambda,uncorr=TRUE,standardize=FALSE)


write.table(ID[rslt$inlier],
	    file=paste(VarFile,".inliers",sep=""),
	    quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(ID[rslt$outlier],
	    file=paste(VarFile,".outliers",sep=""),
	    quote=FALSE,row.names=FALSE,col.names=FALSE)
