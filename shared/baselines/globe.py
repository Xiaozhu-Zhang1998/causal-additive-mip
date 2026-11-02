"""GLOBE fitting and adjacency output for baseline comparisons."""
import numpy as np
import pandas as pd
from pathlib import Path
import tempfile
from Globe.globeWrapper import GlobeWrapper

def fit_globe(X: np.ndarray, *, max_interactions: int=2, log_results: bool=False, verbose: bool=False) -> np.ndarray:
    """
    Fit GLOBE and return a directed adjacency matrix.

    Parameters
    ----------
    X
        Data matrix with shape (n_samples, n_variables).
    max_interactions
        Maximum interaction degree used by the MARS regression.
        The repository example uses 2.
    log_results
        Whether GLOBE writes its internal results to files.
    verbose
        Whether GLOBE prints progress information.

    Returns
    -------
    adjacency
        A p x p NumPy integer matrix. This wrapper assumes
        adjacency[i, j] = 1 means i -> j.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError('X must be a two-dimensional array.')
    n, p = X.shape
    if n < 2 or p < 2:
        raise ValueError('X must contain at least two rows and two columns.')
    if not np.isfinite(X).all():
        raise ValueError('X contains NaN or infinite values.')
    if max_interactions < 1:
        raise ValueError('max_interactions must be at least 1.')
    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = Path(temp_dir) / 'globe_input.csv'
        pd.DataFrame(X, columns=[f'X{j}' for j in range(p)]).to_csv(input_path, index=False)
        globe = GlobeWrapper(max_interactions, log_results, verbose)
        globe.loadData(str(input_path))
        adjacency = globe.run()
    adjacency = np.asarray(adjacency, dtype=int)
    if adjacency.shape != (p, p):
        raise RuntimeError(f'GLOBE returned shape {adjacency.shape}; expected {(p, p)}.')
    adjacency = (adjacency != 0).astype(int)
    np.fill_diagonal(adjacency, 0)
    return adjacency
