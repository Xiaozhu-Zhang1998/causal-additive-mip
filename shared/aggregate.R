# Aggregate headerless R/Python DAGs, retaining raw elapsed time and configuration counts.
args<-commandArgs(trailingOnly=TRUE)
root<-args[1];jobs<-read.csv(args[2],stringsAsFactors=FALSE);dest<-args[3]
source(file.path(root,'shared/pathway.R'))
read_matrix<-function(path) {
  first<-readLines(path,n=1,warn=FALSE)
  A<-as.matrix(read.table(path,header=FALSE,sep=if(grepl(',',first,fixed=TRUE)) ',' else ''))
  storage.mode(A)<-'double'
  if(any(!is.finite(A))) stop(paste('Non-finite matrix:',path))
  A
}
rows<-list();skipped<-list()
for(i in seq_len(nrow(jobs))) {
  j<-jobs[i,,drop=FALSE]
  timing_only<-j$method %in% c('pns','bootstrap')
  if(!timing_only && !file.exists(j$estimate)) {
    skipped[[length(skipped)+1]]<-data.frame(l=j$l,method=j$method,path=j$estimate,reason='missing estimate')
    next
  }
  if(!timing_only) {
    A<-read_matrix(j$estimate)
    if(nrow(A)!=ncol(A)||ncol(A)!=j$p) stop(paste('Wrong graph dimension:',j$estimate))
  }
  if(timing_only) {
    if(!file.exists(j$runtime)) {
      skipped[[length(skipped)+1]]<-data.frame(l=j$l,method=j$method,path=j$runtime,reason='missing timing-only file')
      next
    }
    row<-cbind(j[c('n','p','tau','edge_density','graph','df','gamma','zeta','l','method')],data.frame(shd=NA_real_,sid=NA_real_,f1=NA_real_))
  } else if(j$family=='microarray') {
    row<-cbind(j[c('l','method')],influence_metrics(A))
  } else {
    if(!file.exists(j$truth)) {
      skipped[[length(skipped)+1]]<-data.frame(l=j$l,method=j$method,path=j$truth,reason='missing truth')
      next
    }
    B<-read_matrix(j$truth)
    if(!identical(dim(A),dim(B))) stop('Truth/estimate dimensions differ')
    mask<-row(A)!=col(A);a<-A[mask]!=0;b<-B[mask]!=0
    tp<-sum(a&b);fp<-sum(a&!b);fn<-sum(!a&b)
    sid<-SID::structIntervDist(B,A)$sid
    row<-cbind(j[c('n','p','tau','edge_density','graph','df','gamma','zeta','l','method')],
      data.frame(shd=sum(abs(B-A)),sid=sid,f1=if(2*tp+fp+fn>0) 2*tp/(2*tp+fp+fn) else 0))
    if(j$family=='chamber') {
      row$trial<-j$l;row$directed_f1<-row$f1
      row$tdr<-if(tp+fn>0) tp/(tp+fn) else 0;row$fdr<-if(tp+fp>0) fp/(tp+fp) else 0
      labels<-c(cam='CAM',npvar='NPVAR',mip='MIP',notears='NoTears',globe='Globe',topic='TOPIC',eqvar='EqVar',resit='RESIT',gsges='GSGES')
      row$method<-labels[j$method]
    }
  }
  elapsed<-if(file.exists(j$runtime)) as.numeric(scan(j$runtime,quiet=TRUE,nmax=1)) else NA_real_
  if(is.na(elapsed)) skipped[[length(skipped)+1]]<-data.frame(l=j$l,method=j$method,path=j$runtime,reason='missing runtime; outcome retained')
  row$time<-elapsed;row$time_total<-elapsed
  row$tuning_configurations<-j$tuning_configurations
  row$time_per_configuration<-elapsed/j$tuning_configurations
  rows[[length(rows)+1]]<-row
}
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(length(rows)) write.csv(do.call(rbind,rows),file.path(dest,'RS.csv'),row.names=FALSE)
write.csv(if(length(skipped)) do.call(rbind,skipped) else data.frame(l=integer(),method=character(),path=character(),reason=character()),file.path(dest,'skipped_files.csv'),row.names=FALSE)
cat('Aggregated',length(rows),'outcome rows;',length(skipped),'missing-file records.\n')
