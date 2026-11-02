"""Adjacency differences and cycle checks; SHD counts a reversal twice.

"""
import numpy as np

def edge2adj(edges, p):
    adj = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            if (i, j) in edges:
                adj[i, j] = 1
    return adj

def edge_diff(edges_star, edges_hat, p):
    return np.sum(np.abs(edge2adj(edges_star, p) - edge2adj(edges_hat, p)))

def isCyclicUtil(adj, u, visited, recStack):
    if recStack[u]:
        return True
    if visited[u]:
        return False
    visited[u] = True
    recStack[u] = True
    for v in adj[u]:
        if isCyclicUtil(adj, v, visited, recStack):
            return True
    recStack[u] = False
    return False

def constructadj(V, edges):
    adj = [[] for _ in range(V)]
    for u, v in edges:
        adj[u].append(v)
    return adj

def isCyclic(V, edges):
    adj = constructadj(V, edges)
    visited = [False] * V
    recStack = [False] * V
    for i in range(V):
        if not visited[i] and isCyclicUtil(adj, i, visited, recStack):
            return True
    return False

# From TAM (Apache-2.0); see third_party/TAM_LICENSE.
def compute_caus_order(G):
    d = G.shape[0]
    remain = list(range(d))
    caus_order = np.empty(d, dtype=int)
    for i in range(d - 1):
        root = min(np.where(G.sum(axis=0) == 0)[0])
        caus_order[i] = remain[root]
        del remain[root]
        G = np.delete(G, root, axis=0)
        G = np.delete(G, root, axis=1)
    caus_order[d - 1] = remain[0]
    return caus_order
