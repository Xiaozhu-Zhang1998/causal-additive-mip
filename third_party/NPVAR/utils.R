# Selected NPVAR helpers; original Apache-2.0 license applies.
prune = function(x, est_order, cutoff = 0.001){
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  p = dim(x)[2]
  adj = matrix(0, p, p)
  for (i in 2:p) {
    node = est_order[i]
    ancestors = est_order[1:(i-1)]
    
    mgcvformula = "x[,node] ~ "
    for(a in ancestors){
      mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
    }
    mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
    mgcvformula = as.formula(mgcvformula)
    
    mod = mgcv::gam(mgcvformula)
    parents = ancestors[summary(mod)$s.pv < cutoff]
    adj[parents, node] = 1
  }
  return(adj)
}

### Calculate SHD given adjacency matrix
