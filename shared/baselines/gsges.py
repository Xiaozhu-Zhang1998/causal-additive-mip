"""GES fitting with a generalized score and conversion of a CPDAG to a DAG."""
import numpy as np
from causallearn.search.ScoreBased.GES import ges
from causallearn.utils.PDAG2DAG import pdag2dag

def fit_gsges(X, kfold=5, lam=0.01, maxP=5, return_cpdag=False):
    X = np.asarray(X, dtype=float)
    result = ges(X, score_func='local_score_CV_general', maxP=maxP, parameters={'kfold': kfold, 'lambda': lam})
    cpdag = result['G']
    dag = pdag2dag(cpdag)
    p = X.shape[1]
    adj = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if dag.graph[i, j] == -1 and dag.graph[j, i] == 1:
                adj[i, j] = 1
    if return_cpdag:
        return (adj, cpdag)
    return adj
