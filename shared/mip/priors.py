"""Original heuristic CAM/NPVAR partial-order and stable-edge selection."""
import numpy as np
from shared.utils.graphs import isCyclic

def prior_edges(prob_adj, edge_superset, p):
    # Preserve the supplied thresholding, sorting, and prefix-filtering heuristic.
    par_edges=[]
    par_vals=np.array([])
    kappa_par=.95
    kappa_stab=1
    for i,j in edge_superset:
        if prob_adj[i,j]>=kappa_par and prob_adj[j,i]<=1-kappa_par:
            par_edges.append((i,j));par_vals=np.append(par_vals,prob_adj[i,j])
        if prob_adj[j,i]>=kappa_par and prob_adj[i,j]<=1-kappa_par:
            par_edges.append((j,i));par_vals=np.append(par_vals,prob_adj[j,i])
    idx=np.argsort(-par_vals);par_vals=par_vals[idx];par_edges=[par_edges[i] for i in idx]
    i=1
    while i<len(par_vals):
        if isCyclic(p,par_edges[0:i]):
            par_vals=np.delete(par_vals,i);par_edges.pop(i)
        else:i+=1
    stab_edges=[];stab_vals=np.array([])
    for i in range(p):
        for j in range(p):
            if prob_adj[i,j]>=kappa_stab and (i,j) in par_edges:
                stab_edges.append((i,j));stab_vals=np.append(stab_vals,prob_adj[i,j])
    idx=np.argsort(-stab_vals);stab_vals=stab_vals[idx];stab_edges=[stab_edges[i] for i in idx]
    i=1
    while i<len(stab_vals):
        if isCyclic(p,stab_edges[0:i]):
            stab_vals=np.delete(stab_vals,i);stab_edges.pop(i)
        else:i+=1
    return par_edges,stab_edges
