from __future__ import annotations

import re
from typing import Any

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import default_converter, numpy2ri
from rpy2.robjects.packages import importr


# R package used for multivariate adaptive regression splines.
MARS = importr("earth")

# Converter used only inside explicit `with` blocks.
NUMPY_CONVERTER = default_converter + numpy2ri.converter


def _to_numpy(obj: Any) -> np.ndarray:
    """Convert an R object to a NumPy array using a local converter."""
    with NUMPY_CONVERTER.context():
        return np.asarray(ro.conversion.rpy2py(obj))


def REarth(
    X: np.ndarray,
    Y: np.ndarray,
    M: int = 1,
):
    """
    Fit an R `earth` model and extract the quantities expected by GLOBE.

    Parameters
    ----------
    X
        Predictor matrix with shape (n_samples, n_predictors).
    Y
        Response vector or one-column matrix.
    M
        Maximum interaction degree supplied to `earth`.

    Returns
    -------
    tuple
        (sse, [coeffs], np.array([no_of_terms]), interactions)

    Notes
    -----
    The return format is preserved from the original GLOBE implementation.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    if X.ndim != 2:
        raise ValueError("X must be a two-dimensional array.")

    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    elif Y.ndim != 2 or Y.shape[1] != 1:
        raise ValueError("Y must be a vector or a one-column matrix.")

    if X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must have the same number of rows.")

    if not np.isfinite(X).all() or not np.isfinite(Y).all():
        raise ValueError("X and Y must contain only finite values.")

    if M < 1:
        raise ValueError("M must be at least 1.")

    n_rows, n_cols = X.shape

    # Convert NumPy arrays to R matrices locally.
    with NUMPY_CONVERTER.context():
        rX = ro.r.matrix(
            ro.conversion.py2rpy(X),
            nrow=n_rows,
            ncol=n_cols,
            byrow=False,
        )
        rY = ro.r.matrix(
            ro.conversion.py2rpy(Y),
            nrow=n_rows,
            ncol=1,
            byrow=False,
        )

        try:
            rearth = MARS.earth(
                x=rX,
                y=rY,
                degree=int(M),
            )
        except Exception:
            print(
                "Singular fit encountered; "
                "retrying with Max Interactions = 1."
            )
            rearth = MARS.earth(
                x=rX,
                y=rY,
                degree=1,
            )

    # Indices in the R earth object used by the original GLOBE code.
    RSS_INDEX = 0
    DIRS_INDEX = 5
    CUTS_INDEX = 6
    SELECTED_INDEX = 7
    COEFF_INDEX = 11

    selected = _to_numpy(rearth[SELECTED_INDEX]).reshape(-1)
    no_of_terms = selected.size

    # R uses one-based indices.
    working_index = selected.astype(int) - 1

    dirs_all = _to_numpy(rearth[DIRS_INDEX])
    cuts_all = _to_numpy(rearth[CUTS_INDEX])
    regression_coefficients = _to_numpy(
        rearth[COEFF_INDEX]
    ).reshape(-1, 1)

    dir_rows = dirs_all[working_index, :]
    cut_rows = cuts_all[working_index, :]

    # Retained for compatibility with the original implementation,
    # even though it is not returned.
    dirs = np.sum(np.abs(dir_rows), axis=0)
    unused = (
        len(np.flatnonzero(dirs)) + 1
    ) < X.shape[1]
    del unused

    interactions: list[int] = []

    for row in dir_rows:
        interactions.append(
            int(np.sum(np.asarray(row) != 0))
        )

    coeffs: list[float] = []

    for row in cut_rows:
        nonzero_indices = np.flatnonzero(row)

        for index in nonzero_indices:
            coeffs.append(float(row[index]))

    for value in regression_coefficients[:, 0]:
        coeffs.append(float(value))

    rss = _to_numpy(rearth[RSS_INDEX]).reshape(-1)
    sse = float(rss[0])

    return (
        sse,
        [coeffs],
        np.array([no_of_terms]),
        interactions,
    )


# ------------------------------------------------------------------
# The functions below were commented out in the original file.
# They are preserved as comments for reference.
# ------------------------------------------------------------------

"""
SID_ = importr("SID")
PCALG = importr("pcalg")


def getStructuralDistances(G, H):
    G = np.asarray(G)
    H = np.asarray(H)

    rows, cols = G.shape

    with NUMPY_CONVERTER.context():
        rG = ro.r.matrix(
            ro.conversion.py2rpy(G.reshape(-1)),
            ncol=cols,
            byrow=True,
        )
        rH = ro.r.matrix(
            ro.conversion.py2rpy(H.reshape(-1)),
            ncol=cols,
            byrow=True,
        )

        sid_result = SID_.structIntervDist(rG, rH)
        shd_result = SID_.hammingDist(rG, rH)

    return (
        sid_result[2],
        sid_result[1],
        shd_result[0],
    )
"""