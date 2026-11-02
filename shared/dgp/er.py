"""Erdos-Renyi DAGs with randomly permuted node order.

"""
import numpy as np

def generate_er_dag(p, edge_prob, seed=None, return_order=False):
    """
    Generate a random Erdős-Rényi DAG.

    Parameters
    ----------
    p : int
        Number of nodes.
    edge_prob : float
        Probability of including each admissible forward edge.
        Must be between 0 and 1.
    seed : int or None
        Random seed.
    return_order : bool
        Whether to return the randomly generated topological order.

    Returns
    -------
    adj : ndarray, shape (p, p)
        Directed adjacency matrix, where adj[i, j] = 1 means i -> j.
    edges : list[tuple[int, int]]
        Directed edge list.
    order : ndarray, optional
        Random topological ordering. Returned only if
        return_order=True.
    """
    if not isinstance(p, (int, np.integer)) or p < 1:
        raise ValueError('p must be a positive integer.')
    if not 0 <= edge_prob <= 1:
        raise ValueError('edge_prob must lie between 0 and 1.')
    rng = np.random.default_rng(seed)
    order = rng.permutation(p)
    upper_adj = np.zeros((p, p), dtype=int)
    upper_indices = np.triu_indices(p, k=1)
    upper_adj[upper_indices] = (rng.random(len(upper_indices[0])) < edge_prob).astype(int)
    adj = np.zeros((p, p), dtype=int)
    adj[np.ix_(order, order)] = upper_adj
    edges = [(int(parent), int(child)) for parent, child in np.argwhere(adj == 1)]
    if return_order:
        return (adj, edges, order)
    return (adj, edges)
