# Calculate the Figure 5 scores from saved microarray graph estimates.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
ROOT<-normalizePath(file.path(dirname(this),'..'))
source(file.path(ROOT,'shared/pathway.R'))
result_root<-Sys.getenv('JASA_RESULTS_ROOT',file.path(ROOT,'results'))
rows<-list()
for(l in 0:29) for(method in c('cam','npvar','mip','notears','globe','topic','eqvar','resit','gsges')) {
 f<-file.path(result_root,'Expression/subsamples',sprintf('dag_%s_%d.csv',method,l))
 if(!file.exists(f)) next
 first<-readLines(f,n=1);A<-as.matrix(read.table(f,header=FALSE,sep=if(grepl(',',first,fixed=TRUE)) ',' else ''))
 rows[[length(rows)+1]]<-cbind(data.frame(l=l,method=method),influence_metrics(A))
}
if(!length(rows)) stop('No microarray graph estimates found.')
dest<-file.path(result_root,'Expression');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,rows),file.path(dest,'influence_scores.csv'),row.names=FALSE)
cat('Calculated',length(rows),'microarray metric rows.\n')
