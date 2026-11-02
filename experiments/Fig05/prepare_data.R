# Select the 15 MVA genes in the shared order and standardize all 118 observations.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
ROOT<-normalizePath(dirname(this));while(!file.exists(file.path(ROOT,'shared/settings.py'))) ROOT<-dirname(ROOT)
source(file.path(ROOT,'shared/pathway.R'))
output<-file.path(ROOT,'data/generated/microarray');dir.create(output,recursive=TRUE,showWarnings=FALSE)
path<-file.path(ROOT,'data/real/microarray/13059_2004_896_MOESM1_ESM.txt')
dat<-read.table(path,header=TRUE,skip=6,check.names=FALSE,stringsAsFactors=FALSE)
X<-t(as.matrix(dat[,paste0('c',1:118)]))
colnames(X)<-dat$Genename
X<-X[,genes,drop=FALSE]
stopifnot(identical(colnames(X),genes),identical(dim(X),c(118L,15L)))
write.table(X,file.path(output,'Arabidopsis_MVA_raw.csv'),sep=',',row.names=FALSE,col.names=FALSE)
write.table(scale(X),file.path(output,'Arabidopsis_MVA_standardized.csv'),sep=',',row.names=FALSE,col.names=FALSE)
