"""Centered basis matrices, OLS graph scores and vanilla sufficient statistics.

"""
import numpy as np
import patsy
from shared.utils.graphs import edge2adj

def gen_basis_feat(X, basis_setting):
    basis_type = basis_setting['type']
    if basis_type == 'linear':
        _, p = X.shape
        basis = {}
        for k in range(p):
            B = X[:, k:k + 1].copy()
            B -= B.mean(axis=0)
            basis[k] = B
    elif basis_type == 'spline':
        n_knots = basis_setting['n_knots']
        degree = basis_setting['degree']
        include_intercept = basis_setting['include_intercept']
        space = 1 / (n_knots + 1)
        Q = space * np.arange(1, n_knots + 1)

        def gen_basis(x):
            return patsy.bs(x, knots=np.quantile(x, Q), degree=degree, include_intercept=include_intercept)
        _, p = X.shape
        basis = {}
        for k in range(p):
            B = np.asarray(gen_basis(X[:, k]))
            B -= B.mean(axis=0)
            basis[k] = B
    elif basis_type == 'rbf':
        n_centers = basis_setting['n_centers']
        bandwidth = basis_setting.get('bandwidth', None)
        Q = np.arange(1, n_centers + 1) / (n_centers + 1)
        _, p = X.shape
        basis = {}
        for k in range(p):
            x = X[:, k]
            centers = np.quantile(x, Q)
            if bandwidth is None:
                gaps = np.diff(centers)
                gaps = gaps[gaps > 0]
                h = np.median(gaps) if len(gaps) else np.std(x)
                h = h if h > 0 else 1.0
            else:
                h = float(bandwidth)
            B = np.exp(-0.5 * ((x[:, None] - centers[None, :]) / h) ** 2)
            B -= B.mean(axis=0)
            basis[k] = B
    else:
        raise NotImplementedError(f'Unknown basis type: {basis_type}')
    R = basis[0].shape[1]
    return (basis, R)

def loss_from_edges(edges, X, basis_setting):
    n, p = X.shape
    basis_type = basis_setting['type']
    if basis_type == 'linear':

        def gen_basis(x):
            return x[:, None]
    elif basis_type == 'spline':
        n_knots = basis_setting['n_knots']
        space = 1 / (n_knots + 1)
        Q = space * np.arange(1, n_knots + 1)
        degree = basis_setting['degree']
        include_intercept = basis_setting['include_intercept']

        def gen_basis(x):
            return patsy.bs(x, knots=np.quantile(x, Q), degree=degree, include_intercept=include_intercept)
    elif basis_type == 'rbf':
        n_centers = basis_setting['n_centers']
        bandwidth = basis_setting.get('bandwidth', None)
        Q = np.arange(1, n_centers + 1) / (n_centers + 1)

        def gen_basis(x):
            centers = np.quantile(x, Q)
            if bandwidth is None:
                gaps = np.diff(centers)
                gaps = gaps[gaps > 0]
                h = np.median(gaps) if len(gaps) else np.std(x)
                h = h if h > 0 else 1.0
            else:
                h = bandwidth
            return np.exp(-0.5 * ((x[:, None] - centers) / h) ** 2)
    else:
        raise NotImplementedError(f'Unknown basis type: {basis_type}')
    Adj_mat = edge2adj(edges, p)
    var_hat = np.array([])
    for i in range(p):
        idx = np.nonzero(Adj_mat[:, i])[0]
        if idx.size != 0:
            y_model = X[:, i]
            X_model = X[:, idx]
            Z_model = np.zeros((n, 0))
            for j in range(X_model.shape[1]):
                A_model = gen_basis(X_model[:, j])
                A_model -= np.mean(A_model, axis=0)
                Z_model = np.append(Z_model, A_model, axis=1)
            r = y_model - Z_model @ np.linalg.lstsq(Z_model, y_model, rcond=None)[0]
        else:
            r = X[:, i]
        var_hat = np.append(var_hat, np.square(np.linalg.norm(r, ord=2)) / n)
    return np.sum(np.log(var_hat))

def loss_from_adj(Adj_mat, X, basis_setting):
    n, p = X.shape
    basis_type = basis_setting['type']
    if basis_type == 'spline':
        n_knots = basis_setting['n_knots']
        space = 1 / (n_knots + 1)
        Q = space * np.arange(1, n_knots + 1)
        degree = basis_setting['degree']
        include_intercept = basis_setting['include_intercept']

        def gen_basis(x):
            return patsy.bs(x, knots=np.quantile(x, Q), degree=degree, include_intercept=include_intercept)
    var_hat = np.array([])
    for i in range(p):
        idx = np.nonzero(Adj_mat[:, i])[0]
        if idx.size != 0:
            y_model = X[:, i]
            X_model = X[:, idx]
            p_ = X_model.shape[1]
            Z_model = np.zeros((n, 0))
            for j in range(p_):
                A_model = gen_basis(X_model[:, j])
                A_model = A_model - np.mean(A_model, axis=0)
                Z_model = np.append(Z_model, A_model, axis=1)
            r = y_model - Z_model @ np.linalg.inv(Z_model.T @ Z_model) @ Z_model.T @ y_model
        else:
            r = X[:, i]
        var_hat = np.append(var_hat, np.square(np.linalg.norm(r, ord=2)) / n)
    loss_val = np.sum(np.log(var_hat))
    return loss_val

def gen_suff_stat(X, basis_setting):
    basis_type = basis_setting['type']
    if basis_type == 'spline':
        n_knots = basis_setting['n_knots']
        space = 1 / (n_knots + 1)
        Q = space * np.arange(1, n_knots + 1)
        degree = basis_setting['degree']
        include_intercept = basis_setting['include_intercept']

        def gen_basis(x):
            return patsy.bs(x, knots=np.quantile(x, Q), degree=degree, include_intercept=include_intercept)
    Xsp = {}
    n, p = X.shape
    for j in range(p):
        A = gen_basis(X[:, j])
        A = A - np.mean(A, axis=0)
        Xsp[j] = A
    R = A.shape[1]
    Z = X
    for r in range(R):
        for j in range(p):
            Z = np.append(Z, Xsp[j][:, r].reshape(-1, 1), axis=1)
    Sigma = Z.T @ Z / n
    return {'Sigma': Sigma, 'R': R}
