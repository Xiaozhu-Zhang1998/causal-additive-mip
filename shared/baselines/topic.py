"""TOPIC fitting and adjacency output for baseline comparisons."""
from __future__ import annotations
from typing import Any
import numpy as np
import networkx as nx
from topic.topic.topic import Topic
from topic.topic.scoring.fitting import DataType, ScoreType

def fit_topic(X: np.ndarray, *, extra_refinement: bool=True, verbose: int=0) -> dict[str, Any]:
    """
    Fit TOPIC to i.i.d. continuous data.

    Parameters
    ----------
    X
        Data matrix of shape (n_samples, n_variables).
    extra_refinement
        Whether to run TOPIC's additional parent-set refinement phase.
        The repository default is True.
    verbose
        TOPIC verbosity level. Use 0 for no routine output.

    Returns
    -------
    result
        Dictionary containing the directed adjacency matrix, edge list,
        topological ordering, fitted graph, and raw TOPIC object.

    Convention
    ----------
    adjacency[i, j] = 1 means i -> j.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError('X must be a two-dimensional array.')
    n, p = X.shape
    if n < 2 or p < 2:
        raise ValueError('X must contain at least two observations and two variables.')
    if not np.isfinite(X).all():
        raise ValueError('X contains NaN or infinite values.')
    model = Topic(data_type=DataType.CONTINUOUS, score_type=ScoreType.GAM, extra_refinement=bool(extra_refinement), vb=int(verbose))
    graph, order = model.fit(X)
    if not isinstance(graph, nx.DiGraph):
        raise TypeError(f'TOPIC returned an unexpected graph type: {type(graph).__name__}.')
    if not nx.is_directed_acyclic_graph(graph):
        raise RuntimeError('TOPIC returned a graph containing a cycle.')
    order = [int(node) for node in order]
    if len(order) != p or set(order) != set(range(p)):
        raise RuntimeError('TOPIC returned an invalid topological ordering.')
    adjacency = np.zeros((p, p), dtype=int)
    for parent, child in graph.edges():
        adjacency[int(parent), int(child)] = 1
    np.fill_diagonal(adjacency, 0)
    edges = [(int(parent), int(child)) for parent, child in graph.edges()]
    return {'adjacency': adjacency, 'edges': edges, 'order': order, 'graph': graph, 'model': model}
