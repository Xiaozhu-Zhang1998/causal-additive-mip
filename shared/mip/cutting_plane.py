"""Projected MIP outer approximation; main (5)-(7), Algorithm 1, supplement C.

"""
import numpy as np

def solve_sem_DAG_OA_withinitial1(X, basis, edge_superset, gamma, lam2, partial_order=None, known_edges=None, mip_gap=0.001, eps=1e-06, n_jobs=-1, K_init=10, time_limit=300, verbose=True):
    """Fit one DAG. gamma is inverse ridge strength; lam2 is the edge penalty.

    mip_gap retains the original relative Gurobi MIPGap convention.
    time_limit=None leaves Gurobi unlimited; otherwise applies per optimize call.
    Returns the edge indicators, layer order, and solver objective.
    """
    import numpy as np
    import gurobipy as gp
    from gurobipy import GRB
    from joblib import Parallel, delayed
    partial_order = partial_order or []
    known_edges = set(known_edges or [])
    n, p = X.shape
    BkTXj = {(k, j): basis[k].T @ X[:, j] / n for k in basis for j in range(p)}
    BkTBl = {(k, l): basis[k].T @ basis[l] / n for k in basis for l in basis}
    cache = {j: {} for j in range(p)}

    def compute_node(j, parents):
        key = frozenset(parents)
        if key in cache[j]:
            return cache[j][key]
        Xj = X[:, j]
        if len(parents) == 0:
            alpha_j = Xj / np.sqrt(n)
        else:
            R = basis[parents[0]].shape[1]
            d = len(parents)
            G = np.zeros((d * R, d * R))
            rhs = np.zeros(d * R)
            for a, k in enumerate(parents):
                rhs[a * R:(a + 1) * R] = BkTXj[k, j]
            for a, k in enumerate(parents):
                for b, l in enumerate(parents):
                    G[a * R:(a + 1) * R, b * R:(b + 1) * R] = BkTBl[k, l]
            A = G + np.eye(d * R) / gamma
            beta = np.linalg.solve(A, rhs)
            Zbeta = sum((basis[k] @ beta[a * R:(a + 1) * R] for a, k in enumerate(parents)))
            alpha_j = Xj / np.sqrt(n) - Zbeta / np.sqrt(n)
        val = 0.5 * (Xj @ alpha_j) / np.sqrt(n)
        grad_norms = {k: np.linalg.norm(basis[k].T @ alpha_j) ** 2 for k in basis}
        cache[j][key] = (val, grad_norms)
        return (val, grad_norms)

    def evaluate_g(g_dict):
        parents_by_j = {j: [k2 for (k2, j2), v in g_dict.items() if j2 == j and v == 1] for j in range(p)}
        results = Parallel(n_jobs=n_jobs)((delayed(compute_node)(j, parents_by_j[j]) for j in range(p)))
        fbar = 0.0
        fjs = {}
        for j, (val, _) in enumerate(results):
            fj = max(val, 1e-12)
            fjs[j] = fj
            fbar += np.log(fj)
        fbar += lam2 * sum(g_dict.values())
        grad = {}
        for j, (_, grad_norms) in enumerate(results):
            fj = fjs[j]
            for k2, j2 in g_dict:
                if j2 == j:
                    grad[k2, j] = -(gamma / (2 * n)) * grad_norms[k2] / fj + lam2
        return (fbar, grad)

    def compute_node_fractional(j, g_dict, edge_superset_j, tol=0):
        Xj = X[:, j]
        parents = [k for k in edge_superset_j if g_dict[k, j] > tol]
        weights = {k: g_dict[k, j] for k in parents}
        if len(parents) == 0:
            alpha_j = Xj / np.sqrt(n)
        else:
            R = basis[parents[0]].shape[1]
            d = len(parents)
            G = np.zeros((d * R, d * R))
            rhs = np.zeros(d * R)
            for a, k in enumerate(parents):
                gk = np.sqrt(weights[k])
                rhs[a * R:(a + 1) * R] = gk * BkTXj[k, j]
            for a, k in enumerate(parents):
                gk = np.sqrt(weights[k])
                for b, l in enumerate(parents):
                    gl = np.sqrt(weights[l])
                    G[a * R:(a + 1) * R, b * R:(b + 1) * R] = gk * gl * BkTBl[k, l]
            A = G + np.eye(d * R) / gamma
            beta = np.linalg.solve(A, rhs)
            Zbeta = sum((np.sqrt(weights[k]) * basis[k] @ beta[a * R:(a + 1) * R] for a, k in enumerate(parents)))
            alpha_j = Xj / np.sqrt(n) - Zbeta / np.sqrt(n)
        val = 0.5 * (Xj @ alpha_j) / np.sqrt(n)
        grad_norms = {k: np.linalg.norm(basis[k].T @ alpha_j) ** 2 for k in edge_superset_j}
        return (val, grad_norms)

    def evaluate_g_fractional(g_dict, edge_superset_by_j):
        results = Parallel(n_jobs=n_jobs)((delayed(compute_node_fractional)(j, g_dict, edge_superset_by_j[j]) for j in range(p)))
        fbar = 0.0
        fjs = {}
        for j, (val, _) in enumerate(results):
            fj = max(val, 1e-12)
            fjs[j] = fj
            fbar += np.log(fj)
        fbar += lam2 * sum(g_dict.values())
        grad = {}
        for j, (_, grad_norms) in enumerate(results):
            fj = fjs[j]
            for k in edge_superset_by_j[j]:
                grad[k, j] = -(gamma / (2 * n)) * grad_norms[k] / fj + lam2
        return (fbar, grad)
    model = gp.Model('SEM_DAG_OA')
    g = {e: model.addVar(vtype=GRB.BINARY) for e in edge_superset}
    psi = model.addVars(p, lb=1, ub=p)
    t = model.addVar(lb=-GRB.INFINITY)
    model.setObjective(t, GRB.MINIMIZE)
    for k, j in known_edges:
        if (k, j) in g:
            model.addConstr(g[k, j] == 1)
    for k, j in edge_superset:
        if (j, k) in g:
            model.addConstr(g[k, j] + g[j, k] <= 1)
    for k, j in edge_superset:
        model.addConstr(1 - p + p * g[k, j] <= psi[j] - psi[k])
    for i, j in partial_order:
        model.addConstr(psi[i] <= psi[j])
    model.setObjective(t, GRB.MINIMIZE)
    model.Params.OutputFlag = int(verbose)
    model.Params.LazyConstraints = 1
    model.Params.MIPGap = mip_gap
    if time_limit is not None:
        model.Params.TimeLimit = time_limit
    edge_superset_by_j = {j: [] for j in range(p)}
    for k, j in edge_superset:
        edge_superset_by_j[j].append(k)
    g_val = {e: 1.0 if e in known_edges else np.random.rand() for e in g}
    for e in g:
        g[e].VType = GRB.CONTINUOUS
        g[e].LB = 0
        g[e].UB = 1
    eta = 0
    for _ in range(K_init):
        f_val, grad = evaluate_g_fractional(g_val, edge_superset_by_j)
        expr = f_val
        for e in g:
            if e in known_edges:
                continue
            expr += grad[e] * (g[e] - g_val[e])
        model.addConstr(t >= expr)
        rand_dir = {e: np.random.randn() for e in g if e not in known_edges}
        model.setObjective(t + eta * gp.quicksum((rand_dir[e] * g[e] for e in rand_dir)), GRB.MINIMIZE)
        model.optimize()
        for e in g:
            g_val[e] = 1.0 if e in known_edges else g[e].X
    for e in g:
        g[e].VType = GRB.BINARY

    def oa_callback(model, where):
        if where == GRB.Callback.MIPSOL:
            g_sol = {e: int(round(model.cbGetSolution(g[e]))) for e in g}
            t_sol = model.cbGetSolution(t)
            fval, grad = evaluate_g(g_sol)
            if t_sol < fval - eps:
                expr = fval + gp.quicksum((grad[e] * (g[e] - g_sol[e]) for e in g))
                model.cbLazy(t >= expr)
    model.optimize(oa_callback)
    g_opt = {e: int(round(g[e].X)) for e in g}
    psi_opt = np.array([psi[j].X for j in range(p)])
    return (g_opt, psi_opt, model.ObjVal)
