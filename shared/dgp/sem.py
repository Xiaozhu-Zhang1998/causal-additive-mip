"""Nonlinear additive Gaussian and standardized Student-t SEMs.

"""
import numpy as np
from shared.utils.graphs import compute_caus_order

def gen_data_from_adj(adj_mat, n, p, func, sd=np.array([0.1]), seed=None, return_edge_functions=False):
    """
    Generate data from a nonlinear additive SEM.

    Each generated column is centered to have sample mean zero.

    `func` can be either:
      1. one callable used for every edge; or
      2. a collection of callables, with one function sampled per edge.
    """
    rng = np.random.default_rng(seed)
    sd = np.asarray(sd, dtype=float).reshape(-1)
    if len(sd) == 1:
        sd = np.repeat(sd, p)
    if len(sd) != p:
        raise ValueError('sd must have length 1 or p.')
    caus_order = compute_caus_order(adj_mat)
    if callable(func):
        function_pool = None
    else:
        function_pool = list(func)
        if len(function_pool) == 0:
            raise ValueError('The function pool cannot be empty.')
        if not all((callable(f) for f in function_pool)):
            raise TypeError('Every entry in func must be callable.')
    edge_functions = {}
    for parent, child in np.argwhere(adj_mat == 1):
        parent = int(parent)
        child = int(child)
        if function_pool is None:
            selected_func = func
        else:
            selected_func = function_pool[rng.integers(len(function_pool))]
        edge_functions[parent, child] = selected_func
    X = np.zeros((n, p), dtype=float)
    for i in range(p):
        j = int(caus_order[i])
        parents = np.where(adj_mat[:, j] == 1)[0]
        if len(parents) == 0:
            X[:, j] = rng.normal(scale=sd[j], size=n)
        else:
            Z = np.zeros(n, dtype=float)
            for parent in parents:
                parent = int(parent)
                edge_func = edge_functions[parent, j]
                contribution = np.asarray(edge_func(X[:, parent]), dtype=float).reshape(-1)
                if len(contribution) != n:
                    raise ValueError(f'The function assigned to edge {parent} -> {j} returned the wrong length.')
                Z += contribution
            noise = rng.normal(scale=sd[j], size=n)
            X[:, j] = Z + noise
        X[:, j] -= np.mean(X[:, j])
    if return_edge_functions:
        return (X, edge_functions)
    return X

def linear(x):
    return x

def sine(x):
    return np.sin(x)

def cosine(x):
    return np.cos(x)

def tanh(x):
    return np.tanh(x)

def arctan(x):
    return np.arctan(x)

def soft_cubic(x):
    return np.tanh(x ** 3)

def gen_data_from_adj_tdisn(adj_mat, n, p, func, sd=np.array([0.1]), df=None, seed=None, return_edge_functions=False):
    """
    Generate data from a nonlinear additive SEM.

    Each generated column is centered to have sample mean zero.

    Parameters
    ----------
    adj_mat
        Directed adjacency matrix.
    n
        Sample size.
    p
        Number of variables.
    func
        Either:
        1. one callable used for every edge; or
        2. a collection of callables, with one function sampled per edge.
    sd
        Noise standard deviations. Must have length 1 or p.
    df
        Degrees of freedom of the Student-t noise.
        If None, Gaussian noise is used.
        If specified, must satisfy df > 2. The Student-t noise is
        standardized to variance one before multiplication by sd.
    seed
        Random seed.
    return_edge_functions
        Whether to return the sampled edge-specific functions.
    """
    rng = np.random.default_rng(seed)
    sd = np.asarray(sd, dtype=float).reshape(-1)
    if len(sd) == 1:
        sd = np.repeat(sd, p)
    if len(sd) != p:
        raise ValueError('sd must have length 1 or p.')
    if df is not None and df <= 2:
        raise ValueError('df must be greater than 2 for finite-variance t noise.')

    def sample_noise(j):
        if df is None:
            return rng.normal(scale=sd[j], size=n)
        z = rng.standard_t(df, size=n)
        z *= np.sqrt((df - 2) / df)
        return sd[j] * z
    caus_order = compute_caus_order(adj_mat)
    if callable(func):
        function_pool = None
    else:
        function_pool = list(func)
        if len(function_pool) == 0:
            raise ValueError('The function pool cannot be empty.')
        if not all((callable(f) for f in function_pool)):
            raise TypeError('Every entry in func must be callable.')
    edge_functions = {}
    for parent, child in np.argwhere(adj_mat == 1):
        parent = int(parent)
        child = int(child)
        if function_pool is None:
            selected_func = func
        else:
            selected_func = function_pool[rng.integers(len(function_pool))]
        edge_functions[parent, child] = selected_func
    X = np.zeros((n, p), dtype=float)
    for i in range(p):
        j = int(caus_order[i])
        parents = np.where(adj_mat[:, j] == 1)[0]
        if len(parents) == 0:
            X[:, j] = sample_noise(j)
        else:
            Z = np.zeros(n, dtype=float)
            for parent in parents:
                parent = int(parent)
                edge_func = edge_functions[parent, j]
                contribution = np.asarray(edge_func(X[:, parent]), dtype=float).reshape(-1)
                if len(contribution) != n:
                    raise ValueError(f'The function assigned to edge {parent} -> {j} returned the wrong length.')
                Z += contribution
            X[:, j] = Z + sample_noise(j)
        X[:, j] -= np.mean(X[:, j])
    if return_edge_functions:
        return (X, edge_functions)
    return X

function_pool = [sine, cosine, tanh, arctan, soft_cubic]
