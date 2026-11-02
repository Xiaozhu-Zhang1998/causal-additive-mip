# MVA gene order and known influence graph used to evaluate Figure 5.
genes<-c('AACT1','AACT2','HMGS','HMGR1','HMGR2','MK','MPDC1','MPDC2','IPPI2',
         'FPPS1','FPPS2','DPPS1','DPPS3','GGPPS3','GGPPS4')
E0<-data.frame(from=c('AACT1','AACT2','HMGS','HMGS','HMGR1','HMGR2','MK','MK',
    'MPDC1','MPDC2','IPPI2','IPPI2','IPPI2','IPPI2','FPPS1','FPPS1','FPPS2','FPPS2'),
  to=c('HMGS','HMGS','HMGR1','HMGR2','MK','MK','MPDC1','MPDC2','IPPI2','IPPI2',
    'FPPS1','FPPS2','DPPS1','DPPS3','GGPPS3','GGPPS4','GGPPS3','GGPPS4'))
G0<-igraph::graph_from_data_frame(E0,directed=TRUE,vertices=genes)
D0<-igraph::distances(G0,v=genes,to=genes,mode='out')
R0<-1L*(is.finite(D0)&D0>0)
T0<-R0-t(R0)
# Normalize counts of known reachable pairs by the number of estimated edges.
influence_metrics<-function(A) {
  if(!identical(dim(A),c(length(genes),length(genes))) || any(!is.finite(A)))
    stop('Microarray adjacency must be a finite 15-by-15 matrix in the documented gene order.')
  A<-1L*(A!=0)
  rownames(A)<-colnames(A)<-genes
  G<-igraph::graph_from_adjacency_matrix(A,mode='directed',diag=FALSE)
  D<-igraph::distances(G,mode='out');R<-1L*(is.finite(D)&D>0)
  ind<-which(T0==1,arr.ind=TRUE);rev<-cbind(ind[,2],ind[,1])
  C<-sum(R[ind]);Rv<-sum(R[rev]);N<-nrow(ind);size<-sum(A)
  data.frame(size=size,correct=C,reversed=Rv,missing=N-C-Rv,
             net_influence_score=(C-Rv)/max(size,1),
             influence_score=C/max(size,1))
}
