"""NOTEARS fitting and BIC evaluation for baseline comparisons."""
import copy
import numpy as np
import torch

def notears_bic(fitted_model, X, W_est, max_iter=1000, tolerance_grad=1e-09, tolerance_change=1e-12, include_constant_parameters=False):
    """
    Compute post-selection BIC for a fitted one-hidden-layer NotearsMLP.

    Parameters
    ----------
    fitted_model : NotearsMLP
        Model already fitted by notears_nonlinear().
    X : numpy.ndarray, shape (n, p)
        Data used to fit NOTEARS.
    W_est : numpy.ndarray, shape (p, p)
        Thresholded adjacency-weight matrix returned by notears_nonlinear().
        W_est[i, j] != 0 represents edge i -> j.
    include_constant_parameters : bool
        If False, count only graph-dependent parameters: m * |E|.
        If True, also count biases, output weights, and error variances.

    Returns
    -------
    result : dict
        Contains BIC, RSS values, parameter count, and refitted model.
    """
    X = np.asarray(X)
    n, p = X.shape
    if fitted_model.dims[0] != p:
        raise ValueError('X and fitted_model have incompatible dimensions.')
    if len(fitted_model.dims) != 3 or fitted_model.dims[-1] != 1:
        raise ValueError('This function is written for dims=[p, hidden_units, 1].')
    adjacency = np.asarray(W_est) != 0
    np.fill_diagonal(adjacency, False)
    number_edges = int(adjacency.sum())
    hidden_units = fitted_model.dims[1]
    refitted_model = copy.deepcopy(fitted_model)
    refitted_model.train()
    parameter = next(refitted_model.parameters())
    X_tensor = torch.as_tensor(X, dtype=parameter.dtype, device=parameter.device)
    mask = torch.as_tensor(adjacency.T, dtype=parameter.dtype, device=parameter.device).unsqueeze(1).repeat(1, hidden_units, 1).reshape(p * hidden_units, p)

    def enforce_graph():
        """Keep all excluded first-layer connections equal to zero."""
        with torch.no_grad():
            refitted_model.fc1_pos.weight.mul_(mask)
            refitted_model.fc1_neg.weight.mul_(mask)
    enforce_graph()
    optimizer = torch.optim.LBFGS(refitted_model.parameters(), max_iter=max_iter, tolerance_grad=tolerance_grad, tolerance_change=tolerance_change, line_search_fn='strong_wolfe')

    def closure():
        optimizer.zero_grad()
        enforce_graph()
        X_hat = refitted_model(X_tensor)
        loss = 0.5 / n * torch.sum((X_hat - X_tensor) ** 2)
        loss.backward()
        refitted_model.fc1_pos.weight.grad.mul_(mask)
        refitted_model.fc1_neg.weight.grad.mul_(mask)
        return loss
    optimizer.step(closure)
    enforce_graph()
    refitted_model.eval()
    with torch.no_grad():
        X_hat = refitted_model(X_tensor)
        residuals = X_tensor - X_hat
        rss_node = torch.sum(residuals ** 2, dim=0)
    rss_node = rss_node.detach().cpu().numpy()
    rss_total = float(rss_node.sum())
    tiny = np.finfo(float).tiny
    rss_node_safe = np.maximum(rss_node, tiny)
    negative_twice_loglik = float(n * np.sum(np.log(rss_node_safe / n)))
    graph_df = hidden_units * number_edges
    if include_constant_parameters:
        model_df = graph_df + 2 * p * hidden_units + 2 * p
    else:
        model_df = graph_df
    bic = negative_twice_loglik + model_df * np.log(n)
    return {'bic': float(bic), 'negative_twice_loglik': negative_twice_loglik, 'rss_node': rss_node, 'rss_total': rss_total, 'number_edges': number_edges, 'graph_df': graph_df, 'model_df': model_df, 'adjacency': adjacency}
