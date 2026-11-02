"""Distinct sine/cosine SEM for supplementary Table 2.

"""
import numpy as np
from shared.utils.graphs import edge2adj, compute_caus_order

def gen_data(edges, n, p, func, sd=np.array([0.1])):
    if len(sd) == 1:
        sd = np.repeat(sd, p)
    adj_mat = edge2adj(edges, p)
    caus_order = compute_caus_order(adj_mat)
    X = np.zeros((n, p))
    for i in range(p):
        j = caus_order[i]
        parent = np.where(adj_mat[:, j] == 1)[0]
        if len(parent) == 0:
            X[:, j] = np.random.normal(scale=sd[j], size=n)
        else:
            Z = 0
            for pa in parent:
                Z = Z + func(X[:, pa])
            Z = Z - np.mean(Z, axis=0)
            X[:, j] = Z.reshape(-1) + np.random.normal(scale=sd[j], size=n)
    return X
