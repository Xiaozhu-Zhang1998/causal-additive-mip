# Shared R baseline execution for one saved dataset; called by every experiment adapter.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<7) stop('Internal interface: root mode X.csv output trial superstructure method')
root <- args[1]; mode <- args[2]; input <- args[3]; dest <- args[4]
l <- as.integer(args[5]); super_file <- args[6]; method <- args[7]
source(file.path(root,'third_party/NPVAR/NPVAR.R'))
source(file.path(root,'third_party/NPVAR/utils.R'))
for (f in c('DAG_from_Ordering.R','vselect.R','EqVarDAG_TD.R')) source(file.path(root,'third_party/EqVarDAG',f))
for (f in c('fitting.R','ICML.R','indtestAll.R','indtestHsic.R')) source(file.path(root,'third_party/RESIT',f))
source(file.path(root,'shared/baselines/r_common.R'))
X <- as.matrix(read.csv(input,header=FALSE)); storage.mode(X)<-'double'
colnames(X)<-paste0('X',seq_len(ncol(X)))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
save_fit <- function(name,A,elapsed) {
  write.table(A,file.path(dest,sprintf('dag_%s_%d.csv',name,l)),sep=' ',row.names=FALSE,col.names=FALSE)
  write(elapsed,file.path(dest,sprintf('time_%s_%d.txt',name,l)))
}
if (mode=='bootstrap') {
  selMat <- as.matrix(read.table(super_file,header=FALSE,sep=if(grepl(',',readLines(super_file,n=1))) ',' else ''))
  t <- Sys.time()
  result <- bootstrap_cam_npvar(X,fit_cam_once,fit_npvar_once,B=20,seed=NULL)
  elapsed <- as.numeric(difftime(Sys.time(),t,units='secs'))
  for (name in c('CAM','NPVAR')) write.table(result[[paste0(tolower(name),'_prob_adj')]],
      file.path(dest,sprintf('dag_%spartial_%d.csv',name,l)),sep=' ',row.names=FALSE,col.names=FALSE)
  write(elapsed,file.path(dest,sprintf('time_bootstrap_%d.txt',l)))
} else {
  t <- Sys.time(); pns <- suppressWarnings(cam_pns(X,at_least_selected=.05))
  selMat <- pns$candidate_parents
  save_fit('superstructure',selMat*1L,as.numeric(difftime(Sys.time(),t,units='secs')))
  file.rename(file.path(dest,sprintf('time_superstructure_%d.txt',l)),file.path(dest,sprintf('time_pns_%d.txt',l)))
  for (name in if(mode=='pns'||method=='pns') character() else if(method=='all') c('cam','npvar','eqvar','resit') else method) {
    t <- Sys.time()
    order <- switch(name,
      cam=suppressWarnings(cam_estimate_order(X,selMat,output=FALSE))$order,
      npvar=NPVAR(X,run.np=FALSE,run.mgcv=TRUE,layer.select=FALSE),
      eqvar=EqVarDAG_TD(X,mtd='cvlasso')$TO,
      resit=resit_estimate_order(X,model=train_gam,parsModel=list(numBasisFcts=10),
             indtest=indtestHsic,parsIndtest=list(method='ExactFastTrace'),output=FALSE)$order,
      stop('R methods: cam,npvar,eqvar,resit'))
    A <- gam_prune_from_order(X,order,selMat,cutoff=.001,num_basis_fcts=10,output=FALSE)
    save_fit(name,A,as.numeric(difftime(Sys.time(),t,units='secs')))
  }
}
