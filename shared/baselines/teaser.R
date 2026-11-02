# Figure 1 uses the original default CAM/NPVAR procedures, without common PNS.
args<-commandArgs(trailingOnly=TRUE)
root<-args[1];mode<-args[2];input<-args[3];dest<-args[4]
n<-as.integer(args[5]);l<-as.integer(args[6]);func<-args[7]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
h<-switch(func,'atan(x^2)'=function(x) atan(x^2),'exp(x)'=exp,
          '0.5x^3'=function(x) .5*x^3,'abs(x)'=abs,'sin(x)'=sin,stop('Unknown h'))
truth<-matrix(0,5,5);truth[1,2]<-truth[1,3]<-truth[2,3]<-1
if(mode=='generate') {
  for(noise_regime in c('equal','unequal')) {
    x1<-rnorm(n,0,.5);x2<-x1^2-mean(x1^2)+rnorm(n,0,if(noise_regime=='equal') .5 else .1)
    x3<-2*x1^2+h(x2)-mean(2*x1^2)-mean(h(x2))+rnorm(n,0,if(noise_regime=='equal') .5 else .3)
    X<-cbind(x1,x2,x3,rnorm(n,0,.5),rnorm(n,0,.5))
    write.table(X,file.path(dest,sprintf('%s_%d_%d.csv',noise_regime,n,l)),row.names=FALSE,col.names=FALSE)
  }
} else {
  source(file.path(root,'third_party/NPVAR/NPVAR.R'))
  source(file.path(root,'third_party/NPVAR/utils.R'))
  rows<-list()
  for(noise_regime in c('equal','unequal')) {
    X<-as.matrix(read.table(file.path(input,sprintf('%s_%d_%d.csv',noise_regime,n,l))))
    for(method in c('CAM','NPVAR')) {
      t<-Sys.time()
      A<-if(method=='CAM') CAM::CAM(X,pruning=TRUE)$Adj else prune(X,NPVAR(X))
      elapsed<-as.numeric(difftime(Sys.time(),t,units='secs'))
      write.table(A,file.path(dest,sprintf('dag_%s_%s_%d_%d.csv',method,noise_regime,n,l)),row.names=FALSE,col.names=FALSE)
      rows[[length(rows)+1]]<-data.frame(n=n,m=l,var=noise_regime,method=method,rec=all(A==truth),SHD=sum(abs(A-truth)),func=func,time=elapsed)
    }
  }
  write.csv(do.call(rbind,rows),file.path(dest,sprintf('R_%d_%d.csv',n,l)),row.names=FALSE)
}
