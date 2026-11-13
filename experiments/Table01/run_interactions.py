"""Spline-based DAG learning with optional hierarchical pairwise interactions.

gamma_main and gamma_interaction multiply the quadratic ridge penalties
directly: gamma_main * ||beta||^2 and gamma_interaction * ||zeta||^2.
Both ridge penalty parameters default to one. 
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np
import pandas as pd
import time
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix, vstack
import math
import json
import os
import tempfile
from contextlib import contextmanager
import argparse
import uuid
import signal
import hashlib
ROOT = Path(__file__).resolve().parents[2]

Edge = Tuple[int, int]
Interaction = Tuple[int, int, int]

def _sorted_unique_edges(edges: Iterable[Edge], p: int) -> List[Edge]:
    result = sorted({(int(k), int(j)) for k, j in edges})
    invalid = [(k, j) for k, j in result if k == j or not (0 <= k < p and 0 <= j < p)]
    if invalid:
        raise ValueError(f'Invalid directed edges: {invalid[:5]}')
    return result

def interaction_superset(edge_superset: Sequence[Edge], p: int) -> List[Interaction]:
    """All parent pairs that are allowed by a directed edge superstructure."""
    parents = {j: [] for j in range(p)}
    for k, j in _sorted_unique_edges(edge_superset, p):
        parents[j].append(k)
    return [(l, m, j) for j in range(p) for l, m in combinations(sorted(parents[j]), 2)]

def generate_spline_basis(X: np.ndarray, n_knots: int=3, degree: int=3, include_intercept: bool=False, normalize_columns: bool=True, eps: float=1e-12) -> Tuple[Dict[int, np.ndarray], int]:
    """Generate centered B-spline features for every column of X."""
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError('X must be a two-dimensional array')
    n, p = X.shape
    quantiles = np.arange(1, n_knots + 1) / (n_knots + 1)
    basis: Dict[int, np.ndarray] = {}
    for k in range(p):
        x = X[:, k]
        knots = np.unique(np.quantile(x, quantiles))
        import patsy
        B = np.asarray(patsy.bs(x, knots=knots, degree=degree, include_intercept=include_intercept))
        B = B - B.mean(axis=0, keepdims=True)
        if normalize_columns:
            norms = np.linalg.norm(B, axis=0)
            keep = norms > eps
            B = B[:, keep]
            norms = norms[keep]
            B = B * (np.sqrt(n) / norms)
        basis[k] = B
    dimensions = {B.shape[1] for B in basis.values()}
    if len(dimensions) != 1:
        raise ValueError(f'Variables produced unequal basis dimensions: {dimensions}')
    return (basis, dimensions.pop())

@dataclass
class OracleEvaluation:
    objective: float
    gradient_g: Dict[Edge, float]
    gradient_w: Dict[Interaction, float]
    node_scores: np.ndarray

class InteractionProjectedOracle:
    """Projected ridge score and gradients for main and interaction blocks."""

    def __init__(self, X: np.ndarray, basis: Mapping[int, np.ndarray], edge_superset: Sequence[Edge], gamma_main: float, gamma_interaction: Optional[float]=None, lambda_edge: float=0.0, lambda_pair: float=0.0, fit_intercept: bool=True, center_interactions: bool=True, normalize_interaction_columns: bool=True, include_interactions: bool=True, interaction_design_overrides: Optional[Mapping[Tuple[int, int], np.ndarray]]=None, zero_tol: float=1e-12) -> None:
        self.X = np.asarray(X, dtype=float)
        if self.X.ndim != 2:
            raise ValueError('X must be two dimensional')
        self.n, self.p = self.X.shape
        self.basis = {int(k): np.asarray(v, dtype=float) for k, v in basis.items()}
        if set(self.basis) != set(range(self.p)):
            raise ValueError('basis must contain exactly one matrix for every variable')
        if gamma_main <= 0:
            raise ValueError('gamma_main must be positive')
        if gamma_interaction is not None and gamma_interaction <= 0:
            raise ValueError('gamma_interaction must be positive')
        self.gamma_main = float(gamma_main)
        self.gamma_interaction = float(gamma_main if gamma_interaction is None else gamma_interaction)
        self.lambda_edge = float(lambda_edge)
        self.lambda_pair = float(lambda_pair)
        self.fit_intercept = bool(fit_intercept)
        self.response = self.X - self.X.mean(axis=0, keepdims=True) if self.fit_intercept else self.X
        self.center_interactions = bool(center_interactions)
        self.normalize_interaction_columns = bool(normalize_interaction_columns)
        self.include_interactions = bool(include_interactions)
        self.interaction_design_overrides = {(min(int(pair[0]), int(pair[1])), max(int(pair[0]), int(pair[1]))): np.asarray(design, dtype=float) for pair, design in (interaction_design_overrides or {}).items()}
        for pair, design in self.interaction_design_overrides.items():
            if design.ndim == 1:
                self.interaction_design_overrides[pair] = design[:, None]
            elif design.ndim != 2:
                raise ValueError(f'interaction override {pair} must be one or two dimensional')
            if self.interaction_design_overrides[pair].shape[0] != self.n:
                raise ValueError(f'interaction override {pair} has the wrong number of rows')
        self.zero_tol = float(zero_tol)
        self.edges = _sorted_unique_edges(edge_superset, self.p)
        self.interactions = interaction_superset(self.edges, self.p) if self.include_interactions else []
        self.parents_by_child = {j: [] for j in range(self.p)}
        self.pairs_by_child = {j: [] for j in range(self.p)}
        for k, j in self.edges:
            self.parents_by_child[j].append(k)
        for l, m, j in self.interactions:
            self.pairs_by_child[j].append((l, m))
        self._interaction_design_cache: Dict[Tuple[int, int], np.ndarray] = {}
        self._integer_node_cache: Dict[Tuple[int, frozenset[int]], Tuple[float, Dict[int, float], Dict[Tuple[int, int], float]]] = {}

    def interaction_design(self, l: int, m: int) -> np.ndarray:
        """Return the centered tensor-product basis for variables l and m."""
        key = (min(l, m), max(l, m))
        if key in self._interaction_design_cache:
            return self._interaction_design_cache[key]
        if key in self.interaction_design_overrides:
            C = self.interaction_design_overrides[key].copy()
        else:
            Bl, Bm = (self.basis[key[0]], self.basis[key[1]])
            C = np.einsum('ir,is->irs', Bl, Bm, optimize=True).reshape(self.n, -1)
        if self.center_interactions:
            C = C - C.mean(axis=0, keepdims=True)
        if self.normalize_interaction_columns:
            norms = np.linalg.norm(C, axis=0)
            keep = norms > self.zero_tol
            C = C[:, keep]
            C = C * (np.sqrt(self.n) / norms[keep])
        self._interaction_design_cache[key] = C
        return C

    def _solve_node(self, j: int, main_weights: Mapping[int, float], pair_weights: Mapping[Tuple[int, int], float], allow_cache: bool) -> Tuple[float, Dict[int, float], Dict[Tuple[int, int], float]]:
        candidate_parents = self.parents_by_child[j]
        candidate_pairs = self.pairs_by_child[j]
        binary_main = all((abs(float(main_weights.get(k, 0.0)) - round(float(main_weights.get(k, 0.0)))) <= self.zero_tol for k in candidate_parents))
        active_parent_set = frozenset((k for k in candidate_parents if float(main_weights.get(k, 0.0)) > 0.5))
        expected_pairs = {pair for pair in candidate_pairs if pair[0] in active_parent_set and pair[1] in active_parent_set}
        binary_and_pairs = all((abs(float(pair_weights.get(pair, 0.0)) - float(pair in expected_pairs)) <= self.zero_tol for pair in candidate_pairs))
        cache_key = (j, active_parent_set)
        if allow_cache and binary_main and binary_and_pairs and (cache_key in self._integer_node_cache):
            return self._integer_node_cache[cache_key]
        blocks: List[np.ndarray] = []
        for k in candidate_parents:
            weight = float(main_weights.get(k, 0.0))
            if weight > self.zero_tol:
                blocks.append(np.sqrt(weight / self.gamma_main) * self.basis[k])
        for pair in candidate_pairs:
            weight = float(pair_weights.get(pair, 0.0))
            if weight > self.zero_tol:
                blocks.append(np.sqrt(weight / self.gamma_interaction) * self.interaction_design(*pair))
        y = self.response[:, j]
        if blocks:
            Z = np.column_stack(blocks)
            qdim = Z.shape[1]
            if qdim <= self.n:
                A = Z.T @ Z / self.n + np.eye(qdim)
                rhs = Z.T @ y / self.n
                coef = np.linalg.solve(A, rhs)
                residual = y - Z @ coef
            else:
                H = np.eye(self.n) + Z @ Z.T / self.n
                residual = np.linalg.solve(H, y)
        else:
            residual = y.copy()
        alpha = residual / np.sqrt(self.n)
        q = float(0.5 * y @ residual / self.n)
        if not np.isfinite(q) or q <= 0:
            raise FloatingPointError(f'Nonpositive projected score for node {j}: {q}')
        main_norms = {k: float(np.linalg.norm(self.basis[k].T @ alpha) ** 2) for k in candidate_parents}
        pair_norms = {pair: float(np.linalg.norm(self.interaction_design(*pair).T @ alpha) ** 2) for pair in candidate_pairs}
        result = (q, main_norms, pair_norms)
        if allow_cache and binary_main and binary_and_pairs:
            self._integer_node_cache[cache_key] = result
        return result

    def evaluate(self, g_values: Mapping[Edge, float], w_values: Mapping[Interaction, float], allow_cache: bool=True) -> OracleEvaluation:
        """Evaluate the convex score and its gradient at (g, w)."""
        gradient_g: Dict[Edge, float] = {}
        gradient_w: Dict[Interaction, float] = {}
        node_scores = np.zeros(self.p)
        objective = 0.0
        for j in range(self.p):
            main_weights = {k: float(g_values.get((k, j), 0.0)) for k in self.parents_by_child[j]}
            pair_weights = {(l, m): float(w_values.get((l, m, j), 0.0)) for l, m in self.pairs_by_child[j]}
            q, main_norms, pair_norms = self._solve_node(j, main_weights, pair_weights, allow_cache=allow_cache)
            node_scores[j] = q
            objective += np.log(q)
            for k in self.parents_by_child[j]:
                gradient_g[k, j] = -main_norms[k] / self.gamma_main / (2 * self.n * q) + self.lambda_edge
            for l, m in self.pairs_by_child[j]:
                gradient_w[l, m, j] = -pair_norms[l, m] / self.gamma_interaction / (2 * self.n * q) + self.lambda_pair
        objective += self.lambda_edge * sum((float(g_values.get(e, 0.0)) for e in self.edges))
        objective += self.lambda_pair * sum((float(w_values.get(a, 0.0)) for a in self.interactions))
        return OracleEvaluation(objective, gradient_g, gradient_w, node_scores)

def all_pairs(edges):
    parents = {}
    for k, j in edges:
        parents.setdefault(j, []).append(k)
    return sorted(((k, l, j) for j, ps in parents.items() for k, l in combinations(sorted(ps), 2)))

def screen_support(edges, pairs, points, threshold):
    if len(points) != 1:
        raise ValueError('per-penalty screening requires exactly one relaxation solution')
    kept_e = {e for point in points for e, v in point['g'].items() if v >= threshold}
    kept_w = {a for point in points for a, v in point['w'].items() if v >= threshold}
    for k, l, j in kept_w:
        kept_e.update([(k, j), (l, j)])
    return (sorted(kept_e), sorted(kept_w))

class Master:

    def __init__(self, X, basis, edges, pairs=None, forced=False, relaxed=False, gamma_main=1.0, gamma_interaction=1.0, backend='scipy'):
        self.X, self.basis = (X, basis)
        self.n, self.p = X.shape
        self.edges = sorted(set(map(tuple, edges)))
        self.pairs = all_pairs(self.edges) if pairs is None else sorted(set(map(tuple, pairs)))
        self.forced, self.relaxed, self.backend = (forced, relaxed, backend)
        self.E, self.W = (len(self.edges), len(self.pairs))
        self.size = self.E + self.W + self.p + 1
        self.t = self.size - 1
        self.ei = {e: i for i, e in enumerate(self.edges)}
        self.wi = {a: self.E + i for i, a in enumerate(self.pairs)}
        self.lb = np.zeros(self.size)
        self.ub = np.ones(self.size)
        self.lb[self.E + self.W:self.t] = 1
        self.ub[self.E + self.W:self.t] = self.p
        self.lb[self.t] = -np.inf
        self.ub[self.t] = np.inf
        self.integrality = np.zeros(self.size, dtype=int)
        if not relaxed:
            self.integrality[:self.E] = 1
            if not forced:
                self.integrality[self.E:self.E + self.W] = 1
        self.rows = []
        self.lower = []
        self.upper = []
        for k, j in self.edges:
            self.constraint({self.ei[k, j]: self.p, self.E + self.W + j: -1, self.E + self.W + k: 1}, -np.inf, self.p - 1)
            if k < j and (j, k) in self.ei:
                self.constraint({self.ei[k, j]: 1, self.ei[j, k]: 1}, -np.inf, 1)
        for k, l, j in self.pairs:
            if (k, j) not in self.ei or (l, j) not in self.ei:
                raise ValueError('interaction superset must contain both parent edges')
            wi = self.wi[k, l, j]
            self.constraint({wi: 1, self.ei[k, j]: -1}, -np.inf, 0)
            self.constraint({wi: 1, self.ei[l, j]: -1}, -np.inf, 0)
            if forced:
                self.constraint({wi: -1, self.ei[k, j]: 1, self.ei[l, j]: 1}, -np.inf, 1)
        self.oracle = InteractionProjectedOracle(X, basis, self.edges, gamma_main=gamma_main, gamma_interaction=gamma_interaction, include_interactions=bool(self.pairs))
        self.oracle.interactions = list(self.pairs)
        self.oracle.pairs_by_child = {j: [] for j in range(self.p)}
        for k, l, j in self.pairs:
            self.oracle.pairs_by_child[j].append((k, l))
        self.cuts = 0
        self.add_cut(np.zeros(self.E + self.W))
        if backend == 'gurobi':
            import gurobipy as gp
            self.gp = gp
            self.model = gp.Model('persistent_interaction_oa')
            self.model.Params.OutputFlag = 0
            self.variables = []
            for i in range(self.size):
                self.variables.append(self.model.addVar(lb=float(self.lb[i]), ub=float(self.ub[i]), vtype=gp.GRB.BINARY if self.integrality[i] else gp.GRB.CONTINUOUS))
            self.installed = 0

    def constraint(self, entries, lo, hi):
        row = csr_matrix((list(entries.values()), ([0] * len(entries), list(entries))), shape=(1, self.size))
        self.rows.append(row)
        self.lower.append(lo)
        self.upper.append(hi)

    def evaluate(self, x):
        g = {e: float(x[i]) for e, i in self.ei.items()}
        w = {a: float(x[i]) for a, i in self.wi.items()}
        return self.oracle.evaluate(g, w, allow_cache=False)

    def add_cut(self, x):
        ev = self.evaluate(x)
        gradient = np.zeros(self.size)
        for e, i in self.ei.items():
            gradient[i] = ev.gradient_g[e]
        for a, i in self.wi.items():
            gradient[i] = ev.gradient_w[a]
        intercept = ev.objective - gradient[:self.E + self.W] @ x[:self.E + self.W]
        gradient[self.t] = -1
        self.rows.append(csr_matrix(gradient.reshape(1, -1)))
        self.lower.append(-np.inf)
        self.upper.append(-intercept)
        self.cuts += 1
        return ev

    def linear_solve(self, c, budget, gap):
        if self.backend == 'scipy':
            res = milp(c, integrality=self.integrality, bounds=Bounds(self.lb, self.ub), constraints=LinearConstraint(vstack(self.rows, format='csr'), self.lower, self.upper), options={'time_limit': max(0.01, budget), 'mip_rel_gap': gap})
            x = res.x
            lower = getattr(res, 'mip_dual_bound', None)
            if lower is None and res.status == 0:
                lower = res.fun
            return (x, lower, res.status)
        gp = self.gp
        m = self.model
        while self.installed < len(self.rows):
            i = self.installed
            row = self.rows[i]
            expr = gp.LinExpr([(float(v), self.variables[k]) for k, v in zip(row.indices, row.data) if v])
            if np.isfinite(self.lower[i]):
                m.addConstr(expr >= float(self.lower[i]))
            if np.isfinite(self.upper[i]):
                m.addConstr(expr <= float(self.upper[i]))
            self.installed += 1
        m.setObjective(gp.LinExpr([(float(v), self.variables[i]) for i, v in enumerate(c) if v]))
        m.Params.TimeLimit = max(0.01, budget)
        m.Params.MIPGap = gap
        m.optimize()
        return ([v.X for v in self.variables] if m.SolCount else None, m.ObjBound, m.Status)

    def solve(self, le, li, time_limit=120, tolerance=0.001, max_iterations=1000, initial_support=None):
        start = time.monotonic()
        c = np.zeros(self.size)
        c[:self.E] = le
        c[self.E:self.E + self.W] = li
        c[self.t] = 1
        baseline = np.zeros(self.size)
        baseline[self.E + self.W:self.t] = np.arange(1, self.p + 1)
        baseline[self.t] = self.evaluate(baseline).objective
        best = baseline.copy()
        upper = float(c @ best)
        lower = -np.inf
        if initial_support is not None:
            if self.relaxed:
                raise ValueError('initial_support is for integer solves')
            edges, pairs = map(set, initial_support)
            if not edges <= set(self.ei) or not pairs <= set(self.wi):
                raise ValueError('initial support lies outside the candidate set')
            if any(((k, j) not in edges or (l, j) not in edges for k, l, j in pairs)):
                raise ValueError('initial interactions require both parent edges')
            if self.forced and pairs != {a for a in self.pairs if (a[0], a[2]) in edges and (a[1], a[2]) in edges}:
                raise ValueError('initial support violates forced interaction constraints')
            parents = {j: {k for k, child in edges if child == j} for j in range(self.p)}
            order = []
            while len(order) < self.p:
                ready = sorted((j for j in range(self.p) if j not in order and parents[j] <= set(order)))
                if not ready:
                    raise ValueError('initial edge support contains a cycle')
                order.extend(ready)
            guess = baseline.copy()
            for e in edges:
                guess[self.ei[e]] = 1
            for a in pairs:
                guess[self.wi[a]] = 1
            for rank, j in enumerate(order, 1):
                guess[self.E + self.W + j] = rank
            ev = self.add_cut(guess)
            guess[self.t] = ev.objective
            value = float(c @ guess)
            if value < upper:
                best = guess.copy()
                upper = value
            if self.backend == 'gurobi':
                for variable, value in zip(self.variables, guess):
                    variable.Start = float(value)
        status = 'time_limit'
        iterations = 0
        for iterations in range(1, max_iterations + 1):
            remaining = time_limit - (time.monotonic() - start)
            if remaining <= 0:
                break
            x, bound, _ = self.linear_solve(c, remaining, tolerance * 0.2)
            if bound is not None and np.isfinite(bound):
                lower = max(lower, float(bound))
            if x is None:
                break
            x = np.array(x)
            if not self.relaxed:
                x[:self.E] = np.rint(x[:self.E])
                x[self.E:self.E + self.W] = np.rint(x[self.E:self.E + self.W])
            ev = self.add_cut(x)
            actual = ev.objective + le * x[:self.E].sum() + li * x[self.E:self.E + self.W].sum()
            if actual < upper:
                best = x.copy()
                best[self.t] = ev.objective
                upper = actual
            certificate = max(0, upper - lower) / max(1, abs(upper))
            if certificate <= tolerance:
                status = 'optimal'
                break
        if iterations == max_iterations and status != 'optimal':
            status = 'iteration_limit'
        return {'g': {e: float(best[i]) for e, i in self.ei.items()}, 'w': {a: float(best[i]) for a, i in self.wi.items()}, 'objective': upper, 'lower_bound': lower, 'gap': max(0, upper - lower) / max(1, abs(upper)), 'status': status, 'seconds': time.monotonic() - start, 'iterations': iterations, 'cuts': self.cuts}

def score_graph(X, basis, edges, pairs, gm, gi, only_nominal=False):
    oracle = InteractionProjectedOracle(X, basis, edges, gamma_main=gm, gamma_interaction=gi, include_interactions=bool(pairs))
    n, p = X.shape
    rb = ob = edf = dim = 0.0
    admissible = True
    for j in range(p):
        blocks = [basis[k] for k, c in edges if c == j]
        blocks += [oracle.interaction_design(k, l) for k, l, c in pairs if c == j]
        scales = [gm for k, c in edges if c == j] + [gi for k, l, c in pairs if c == j]
        y = X[:, j] - X[:, j].mean()
        if blocks:
            H = np.column_stack(blocks)
            d = H.shape[1]
            if d > 0.8 * (n - 1) and only_nominal:
                admissible = False
                ob = np.inf
                dim += d
                continue
            Z = np.column_stack([B / np.sqrt(s) for s, B in zip(scales, blocks)])
            if only_nominal:
                r = y
                df = 0.0
            elif Z.shape[1] <= n:
                gram = Z.T @ Z / n
                system = np.eye(Z.shape[1]) + gram
                a = np.linalg.solve(system, Z.T @ y / n)
                r = y - Z @ a
                df = float(np.trace(np.linalg.solve(system, gram)))
            else:
                A = np.eye(n) + Z @ Z.T / n
                inv = np.linalg.solve(A, np.eye(n))
                r = inv @ y
                df = float(n - np.trace(inv))
            if d <= 0.8 * (n - 1):
                coef, _, rank, _ = np.linalg.lstsq(H, y, rcond=None)
                ro = y - H @ coef
                admissible &= rank == d
            else:
                admissible = False
                ro = y
                ob = np.inf
        else:
            d = df = 0
            r = ro = y
        rb += n * np.log(max(r @ r / n, np.finfo(float).tiny)) + df * np.log(n)
        ob += n * np.log(max(ro @ ro / n, np.finfo(float).tiny)) + d * np.log(n)
        edf += df
        dim += d
    return {'ridge_bic': np.nan if only_nominal else rb, 'ols_bic': ob, 'effective_df': np.nan if only_nominal else edf, 'nominal_parameters': dim, 'ols_admissible': bool(admissible)}

def provisional_dag(edges, pairs, weights):
    """Keep high-weight edges first, skipping those that would create a cycle."""
    children = {}
    kept = set()
    for k, j in sorted(edges, key=lambda e: (-weights.get(e, 0.0), e)):
        pending = [j]
        seen = set()
        while pending:
            node = pending.pop()
            if node in seen:
                continue
            seen.add(node)
            pending.extend(children.get(node, ()))
        if k in seen:
            continue
        kept.add((k, j))
        children.setdefault(k, set()).add(j)
    interactions = [a for a in pairs if (a[0], a[2]) in kept and (a[1], a[2]) in kept]
    return (sorted(kept), sorted(interactions))

def rank_candidates(candidates, size):
    """Union of top nominal-BIC and effective-df-BIC models. Fit each once."""
    for c in candidates:
        c.update(shortlist_rank=None, shortlisted=False, nominal_shortlist_rank=None, effective_df_shortlist_rank=None, shortlisted_by_nominal=False, shortlisted_by_effective_df=False)
    for column, rank_key, flag, nominal in [('ols_bic', 'nominal_shortlist_rank', 'shortlisted_by_nominal', True), ('ridge_bic', 'effective_df_shortlist_rank', 'shortlisted_by_effective_df', False)]:
        eligible = [c for c in candidates if math.isfinite(c[column]) and (not nominal or c['ols_admissible'])]
        eligible.sort(key=lambda c: (c[column], c['screened_interactions'], c['screened_edges'], -c['lambda_edge'], c['rho']))
        for rank, c in enumerate(eligible, 1):
            c[rank_key] = rank
            if rank <= size:
                c[flag] = True
                c['shortlisted'] = True
                c['shortlist_rank'] = min(c['shortlist_rank'] or rank, rank)
    return sorted((c for c in candidates if c['shortlisted']), key=lambda c: (c['shortlist_rank'], -c['lambda_edge'], c['rho']))

def iter_jsonl(path):
    path = Path(path)
    if not path.exists():
        return
    with path.open('r+b') as handle:
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                if handle.read(1):
                    raise
                handle.seek(offset)
                handle.truncate()
                print(f'Removed incomplete final record from {path.name}', flush=True)
                break
            if not line.endswith(b'\n'):
                handle.write(b'\n')
            yield record

def read_jsonl(path):
    return list(iter_jsonl(path))

def append(path, record):
    path = Path(path)
    if path.exists() and path.stat().st_size:
        with path.open('rb') as tail:
            tail.seek(-1, 2)
            needs_newline = tail.read(1) != b'\n'
        if needs_newline:
            with path.open('ab') as handle:
                handle.write(b'\n')
    with Path(path).open('a', encoding='utf-8') as handle:
        handle.write(json.dumps(record) + '\n')
        handle.flush()
        os.fsync(handle.fileno())

def _atomic_refresh(path, write, checkpoint_index=False):
    """Derived-file I/O must not abort fits already saved in the raw journals."""
    path = Path(path)
    temporary = None
    descriptor = None
    try:
        descriptor, name = tempfile.mkstemp(prefix=path.name + '.tmp.', dir=path.parent)
        temporary = Path(name)
        handle = os.fdopen(descriptor, 'w', encoding='utf-8', newline='')
        descriptor = None
        with handle:
            write(handle)
        os.replace(temporary, path)
        return True
    except OSError as error:
        if checkpoint_index:
            print(f'Cannot refresh {path.name}: {error}. Completed fits remain saved in checkpoint.jsonl and its index will be rebuilt on resume.', flush=True)
        else:
            print(f'Cannot refresh {path.name}: {error}. All results remain saved in the checkpoints. Continuing the study. Use --export-only after resolving the file access problem.', flush=True)
        return False
    finally:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except OSError:
                pass
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass

def atomic_csv(path, rows):
    return _atomic_refresh(path, lambda handle: pd.DataFrame(rows).to_csv(handle, index=False))

@contextmanager
def phase(path, run_id, condition_id, replicate, stage):
    start = time.perf_counter()
    status = 'complete'
    try:
        yield
    except BaseException:
        status = 'interrupted'
        raise
    finally:
        append(path, dict(run_id=run_id, condition_id=condition_id, replicate=replicate, stage=stage, seconds=time.perf_counter() - start, status=status))

def encode(point):
    result = {k: v for k, v in point.items() if k not in ('g', 'w')}
    result['g'] = [[list(e), v] for e, v in point['g'].items()]
    result['w'] = [[list(a), v] for a, v in point['w'].items()]
    return result

def decode(point):
    return dict(point, g={tuple(e): v for e, v in point['g']}, w={tuple(a): v for a, v in point['w']})

class PointCache:
    """Index raw points on disk instead of keeping all fractional arrays in RAM."""

    def __init__(self, out):
        self.path = Path(out) / 'checkpoint.jsonl'
        self.index_path = Path(out) / 'checkpoint_index.json'
        self.offsets = {}
        start = 0
        if self.index_path.exists():
            index = json.loads(self.index_path.read_text())
            if self.path.exists() and index['bytes'] <= self.path.stat().st_size:
                self.offsets = index['offsets']
                start = index['bytes']
        if self.path.exists():
            with self.path.open('r+b') as handle:
                handle.seek(start)
                while True:
                    offset = handle.tell()
                    line = handle.readline()
                    if not line:
                        break
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError:
                        if handle.read(1):
                            raise
                        handle.seek(offset)
                        handle.truncate()
                        break
                    if not line.endswith(b'\n'):
                        handle.write(b'\n')
                        line += b'\n'
                    self.offsets[record['key']] = [offset, len(line)]

    def __contains__(self, key):
        return key in self.offsets

    def get(self, key):
        offset, size = self.offsets[key]
        with self.path.open('rb') as handle:
            handle.seek(offset)
            return decode(json.loads(handle.read(size)))

    def put(self, key, point):
        data = (json.dumps(encode(point)) + '\n').encode('utf-8')
        with self.path.open('ab') as handle:
            offset = handle.tell()
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        self.offsets[key] = [offset, len(data)]

    def flush(self):
        payload = json.dumps(dict(bytes=self.path.stat().st_size if self.path.exists() else 0, offsets=self.offsets))
        return _atomic_refresh(self.index_path, lambda handle: handle.write(payload), checkpoint_index=True)

def export(out):
    out = Path(out)
    candidates = {r['key']: r for r in iter_jsonl(out / 'rows_checkpoint.jsonl')}
    screens = {r['key']: r for r in iter_jsonl(out / 'screening_checkpoint.jsonl')}
    shortlisted = {r['final_key'] for r in screens.values() if r.get('shortlisted')}
    rows = [r for r in candidates.values() if r['key'] in shortlisted]
    atomic_csv(out / 'candidates.csv', rows)
    atomic_csv(out / 'screening.csv', list(screens.values()))
    selected = []
    if rows:
        data = pd.DataFrame(rows)
        for criterion, column in [('ols_nominal_parameters', 'ols_bic'), ('ridge_effective_df', 'ridge_bic')]:
            eligible = data[data.ols_admissible] if criterion == 'ols_nominal_parameters' else data
            for _, group in eligible.groupby(['condition_id', 'replicate', 'approach', 'threshold'], dropna=False):
                row = group.loc[group[column].idxmin()].to_dict()
                row.update(selection_criterion=criterion, selection_bic=row[column])
                selected.append(row)
    atomic_csv(out / 'selected.csv', selected)
    if selected:
        data = pd.DataFrame(selected)
        metrics = ['sid', 'shd', 'shd_skeleton', 'exact_dag', 'edge_tpr', 'edge_fdr', 'interaction_tpr', 'interaction_fdr', 'exact_interaction_support', 'true_interactions_recovered', 'false_interactions', 'all_true_interactions_recovered', 'exact_motif_parent_sets', 'all_motif_parents_recovered', 'motif_parent_recall', 'motif_reversed_edges', 'motif_missing_adjacencies', 'exact_target_parents', 'target_parent_recall', 'third_parent_selected', 'third_parent_reversed', 'third_parent_screen_retained', 'unnecessary_third_parent_interactions', 'estimated_edges', 'estimated_interactions', 'screened_edges', 'screened_interactions', 'gap', 'seconds', 'nominal_parameters', 'effective_df', 'lambda_edge', 'rho']
        groups = ['condition_id', 'block', 'p', 'n', 'interaction_count', 'interaction_strength', 'weak_strength', 'third_strength', 'setting_number', 'main_form', 'noise_tau', 'approach', 'threshold', 'selection_criterion']
        summary = data.groupby(groups, dropna=False)[metrics].mean()
        summary['replicates'] = data.groupby(groups, dropna=False).size()
        atomic_csv(out / 'summary.csv', summary.reset_index().to_dict('records'))
        paired = []
        for (cid, rep, criterion, threshold), group in data.groupby(['condition_id', 'replicate', 'selection_criterion', 'threshold'], dropna=False):
            baseline = group[group.approach == 'main_effects_only']
            if baseline.empty:
                continue
            baseline = baseline.iloc[0]
            for _, row in group[group.approach == 'optional_interactions'].iterrows():
                paired.append(dict(condition_id=cid, replicate=rep, selection_criterion=criterion, threshold=row.threshold, baseline_shd=baseline.shd, interaction_shd=row.shd, delta_shd=row.shd - baseline.shd, delta_bic=row.selection_bic - baseline.selection_bic, baseline_exact_dag=baseline.exact_dag, interaction_exact_dag=row.exact_dag))
        atomic_csv(out / 'paired_comparison.csv', paired)
    events = read_jsonl(out / 'timing_events.jsonl')
    if events:
        times = pd.DataFrame(events).groupby(['condition_id', 'replicate', 'stage']).seconds.sum().unstack(fill_value=0)
        times.columns = [c + '_seconds' for c in times.columns]
        atomic_csv(out / 'timings.csv', times.reset_index().to_dict('records'))
VERSION = 'interaction_study_v1'
BLOCKS = ['interaction', 'additive', 'pure', 'three_parent']

@dataclass(frozen=True)
class Condition:
    condition_id: str
    block: str
    p: int
    n: int
    interaction_count: int
    interaction_strength: float
    weak_strength: float = 0.5
    third_strength: float = 0.0
    setting_number: int = 0
    main_form: str = 'nonlinear'
    noise_tau: float = 0.1

    def metadata(self):
        return asdict(self)

def token(value):
    return format(float(value), '.8g').replace('.', 'p').replace('-', 'm')

def make_conditions(args):
    specs = [('interaction', 1, e, 0.5, 0.0, 'nonlinear') for e in (0.5, 1.0, 2.0)]
    specs += [('interaction', 3, e, 0.5, 0.0, 'nonlinear') for e in (0.5, 1.0, 2.0)]
    specs += [('interaction', 3, e, 0.5, 0.0, 'linear') for e in (0.5, 2.0)]
    specs += [('additive', 0, 0.0, 0.5, 0.0, 'nonlinear'), ('pure', 1, 1.0, 0.5, 0.0, 'nonlinear')]
    specs += [('three_parent', 1, 1.0, 0.5, b, 'nonlinear') for b in (0.2, 0.5)]
    specs += [('interaction', 3, 1.0, 0.5, 0.0, 'strong_cos')]
    result = []
    for tau in args.noise_taus:
        for number, (block, k, eta, b, third, form) in enumerate(specs):
            name = f's{number:02d}_{block}_p8_n300_k{k}_eta{token(eta)}_{form}_tau{token(tau)}'
            if block == 'three_parent':
                name += f'_b{token(third)}'
            result.append(Condition(name, block, 8, 300, k, eta, b, third, number, form, tau))
    if len({c.condition_id for c in result}) != len(result):
        raise ValueError('duplicate noise levels')
    if args.condition_ids:
        unknown = set(args.condition_ids) - {c.condition_id for c in result}
        if unknown:
            raise ValueError('unknown condition IDs: ' + ', '.join(sorted(unknown)))
        result = [c for c in result if c.condition_id in args.condition_ids]
    return result

def standardize(v):
    v = np.asarray(v, dtype=float)
    sd = v.std(ddof=0)
    if sd <= 1e-12:
        raise ValueError('constant generated component')
    return (v - v.mean()) / sd
FUNCTIONS = ('sin', 'cos', 'tanh', 'arctan', 'tanh_cube')

def raw_function(v, name):
    if name == 'sin':
        return np.sin(v)
    if name == 'cos':
        return np.cos(v)
    if name == 'tanh':
        return np.tanh(v)
    if name == 'arctan':
        return np.arctan(v)
    if name == 'tanh_cube':
        return np.tanh(v ** 3)
    if name == 'linear':
        return v
    raise ValueError(name)
INTERACTION_FUNCTIONS = ('product', 'cos_sum', 'sin_product', 'smooth_ratio')

def interaction_choice(condition, rep, seed, k, l, j):
    """Sample uniformly from four surfaces using independent seeded streams."""
    if condition.block == 'pure':
        return 'product'
    tag = k * condition.p * condition.p + l * condition.p + j
    slot = int(rng(seed, rep, 107, tag).integers(5))
    index = int(rng(seed, rep, 109, tag).integers(4)) if slot == 3 else {0: 0, 1: 1, 2: 2, 4: 3}[slot]
    return INTERACTION_FUNCTIONS[index]

def raw_interaction(u, v, name, reverse=False):
    if name == 'product':
        return u * v
    if name == 'cos_sum':
        return np.cos(u + v)
    if name == 'sin_product':
        return np.sin(u * v)
    if name == 'smooth_ratio':
        if reverse:
            u, v = (v, u)
        return u / np.hypot(v, 0.5)
    raise ValueError(name)

def nonlinear(v):
    return standardize(np.cos(v))

def rng(seed, rep, tag, node=0):
    return np.random.default_rng(np.random.SeedSequence([seed, rep, tag, node]))

def simulate(condition, rep, seed=1, noise_sd=1.0):
    p, n = (condition.p, condition.n)
    log_scale = condition.noise_tau * rng(seed, rep, 105).normal(size=p)
    raw = np.exp(log_scale - log_scale.max())
    sigma = noise_sd * raw / np.sqrt(np.mean(raw ** 2))
    X = np.zeros((n, p))
    X[:, :3] = rng(seed, rep, 101).normal(size=(n, 3)) * sigma[:3]
    errors = {j: rng(seed, rep, 102, j).normal(size=n) * sigma[j] for j in range(3, p)}
    parent_sets = {3: [0, 1], 4: [0, 2], 5: [1, 2]}
    if condition.block == 'three_parent':
        parent_sets[3] = [0, 1, 2]
    truth = []
    true_interactions = []
    components = []
    function_assignments = {}
    diagnostics = []
    interaction_assignments = []

    def component(k, j, linear=False):
        if linear:
            name = 'linear'
        elif condition.main_form == 'strong_cos':
            name = 'cos'
        else:
            name = FUNCTIONS[int(rng(seed, rep, 106, k * p + j).integers(len(FUNCTIONS)))]
        value = standardize(raw_function(X[:, k], name))
        x = X[:, k] - X[:, k].mean()
        r2 = float((x @ value) ** 2 / (x @ x * (value @ value)))
        function_assignments[f'{k}->{j}'] = name
        diagnostics.append(dict(parent=k, child=j, function=name, linear_R2=r2))
        return value
    for index, j in enumerate([3, 4, 5]):
        parents = parent_sets[j]
        active = index < condition.interaction_count
        eta = condition.interaction_strength if active else 0.0
        mains = [0.5, 0.5]
        if condition.block == 'additive' and j == 3:
            mains = [1.0, condition.weak_strength]
        if condition.block == 'pure' and j == 3:
            mains = [0.0, 0.0]
        form = 'linear' if active and condition.main_form == 'linear' else 'nonlinear'
        main = sum((a * component(k, j, form == 'linear') for a, k in zip(mains, parents[:2]) if a != 0.0))
        if len(parents) == 3:
            main += condition.third_strength * component(parents[2], j)
        interaction = np.zeros(n)
        if active:
            k, l = parents[:2]
            tag = k * p * p + l * p + j
            name = interaction_choice(condition, rep, seed, k, l, j)
            reverse = bool(rng(seed, rep, 108, tag).integers(2)) if name == 'smooth_ratio' else False
            surface = raw_interaction(standardize(X[:, k]), standardize(X[:, l]), name, reverse)
            interaction = eta * standardize(surface)
            interaction_assignments.append(dict(parents=[k, l], child=j, function=name, reverse_inputs=reverse, numerator_parent=(l if reverse else k) if name == 'smooth_ratio' else None, denominator_parent=(k if reverse else l) if name == 'smooth_ratio' else None))
        X[:, j] = main + interaction + errors[j]
        truth.extend(((k, j) for k in parents))
        if active:
            true_interactions.append((parents[0], parents[1], j))
        components.append(dict(node=j, parents=parents, main_strengths=mains, main_form=form, third_main_strength=condition.third_strength if len(parents) == 3 else 0.0, third_main_form='nonlinear' if len(parents) == 3 else None, interaction_strength=eta, main_sd=float(np.std(main)), interaction_sd=float(np.std(interaction))))
    for j in range(6, p):
        parent = int(rng(seed, rep, 103, j).integers(3, j))
        X[:, j] = component(parent, j) + errors[j]
        truth.append((parent, j))
        parent_sets[j] = [parent]
    edges = sorted(truth)
    pairs = sorted(true_interactions)
    third = (2, 3) if condition.block == 'three_parent' else None
    metadata = dict(condition=condition.metadata(), replicate=rep, seed=seed, noise_sd=noise_sd, noise_tau=condition.noise_tau, noise_standard_deviations=sigma.tolist(), function_assignments=function_assignments, nonlinearity_diagnostics=diagnostics, interaction_assignments=interaction_assignments, interaction_pool=list(INTERACTION_FUNCTIONS), dgp_version=VERSION, column_order='X0,...,X7', graph_index_base=0, original_to_observed=list(range(p)), true_edges=edges, true_interactions=pairs, host_parents={str(j): parent_sets[j] for j in [3, 4, 5]}, target_node=3, third_parent_edge=third, unused_third_parent_interactions=[(0, 2, 3), (1, 2, 3)] if third else [], components=components)
    return (X, metadata)

def recovery(metadata, edges, pairs, screen_edges, screen_pairs):
    truth = set(map(tuple, metadata['true_edges']))
    ti = set(map(tuple, metadata['true_interactions']))
    edges = set(edges)
    pairs = set(pairs)
    ke = set(screen_edges)
    kw = set(screen_pairs)
    skeleton = lambda es: {tuple(sorted(e)) for e in es}
    hosts = {int(j): set(ps) for j, ps in metadata['host_parents'].items()}
    motif = {(k, j) for j, ps in hosts.items() for k in ps}
    target = metadata['target_node']
    target_parents = {k for k, j in edges if j == target}
    result = dict(shd=len(truth ^ edges), shd_skeleton=len(skeleton(truth) ^ skeleton(edges)), exact_dag=int(edges == truth), edge_tpr=len(truth & edges) / len(truth), edge_fdr=len(edges - truth) / len(edges) if edges else 0.0, interaction_tpr=len(ti & pairs) / len(ti) if ti else np.nan, interaction_fdr=len(pairs - ti) / len(pairs) if pairs else 0.0, exact_interaction_support=int(pairs == ti), true_interactions_recovered=len(ti & pairs), false_interactions=len(pairs - ti), all_true_interactions_recovered=int(ti <= pairs) if ti else np.nan, exact_motif_parent_sets=sum(({k for k, j2 in edges if j2 == j} == ps for j, ps in hosts.items())) / len(hosts), all_motif_parents_recovered=int(motif <= edges), motif_parent_recall=len(motif & edges) / len(motif), motif_reversed_edges=sum(((j, k) in edges for k, j in motif)), motif_missing_adjacencies=sum(((k, j) not in edges and (j, k) not in edges for k, j in motif)), exact_target_parents=int(target_parents == hosts[target]), target_parent_recall=len(target_parents & hosts[target]) / len(hosts[target]), target_true_parent_count=len(hosts[target]))
    third = metadata['third_parent_edge']
    result.update(third_parent_selected=int(tuple(third) in edges) if third else np.nan, third_parent_reversed=int(tuple(third[::-1]) in edges) if third else np.nan, third_parent_screen_retained=int(tuple(third) in ke) if third else np.nan, unnecessary_third_parent_interactions=len(pairs & set(map(tuple, metadata['unused_third_parent_interactions']))) if third else np.nan)
    return result
PIPELINE_VERSION = 'dual_bic_shortlist_union_v2'
METHODS = ('main_effects_only', 'optional_interactions')

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summarize-results', type=Path, help='Read graph estimates and recompute their metrics in a separate directory')
    p.add_argument('--noise-taus', nargs='+', type=float, default=[0.1, 0.2])
    p.add_argument('--condition-ids', nargs='+')
    p.add_argument('--settings', nargs='+', type=int, help='Setting numbers, 0 through 12')
    trials = p.add_mutually_exclusive_group()
    trials.add_argument('--replicates', type=int, default=20)
    trials.add_argument('--trial', type=int, choices=range(20), metavar='0..19', help='Run one trial in each selected condition')
    p.add_argument('--save-data-only', action='store_true', help='Generate observations and ground truth without fitting')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--noise-sd', type=float, default=1.0)
    p.add_argument('--lambda-values', nargs='+', type=float, help='Optional absolute grid overriding the dimension-scaled default')
    p.add_argument('--lambda-multipliers', nargs='+', type=float, default=[1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125], help='Default lambda_E = 2*p*log(p)^2/n times these multipliers')
    p.add_argument('--rho-values', nargs='+', type=float, default=[1.0, 2.0])
    p.add_argument('--shortlist-size', type=int, default=2, help='Top candidates PER BIC and method. Solve the union once (up to twice this number).')
    p.add_argument('--thresholds', nargs='+', type=float, default=[0.3])
    p.add_argument('--gamma-main', type=float, default=1.0)
    p.add_argument('--gamma-interaction', type=float, default=1.0)
    p.add_argument('--n-knots', type=int, default=3)
    p.add_argument('--degree', type=int, default=3)
    p.add_argument('--time-limit', type=float, default=300.0)
    p.add_argument('--relax-time-limit', type=float, default=120.0)
    p.add_argument('--gap', type=float, default=0.001)
    p.add_argument('--relax-gap', type=float, default=0.001, help='OA relative-gap tolerance for BOTH continuous screening methods')
    p.add_argument('--backend', choices=['gurobi', 'scipy'], default='gurobi')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--output-dir', type=Path, help='Directory for generated datasets or results')
    p.add_argument('--export-every', type=int, default=1, help='Refresh CSV tables every this many datasets; checkpoints remain per fit')
    p.add_argument('--list-conditions', action='store_true')
    p.add_argument('--validate-only', action='store_true')
    p.add_argument('--export-only', action='store_true')
    return p

def save_dataset(out, condition, rep, X, metadata):
    folder = out / 'data'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f'{condition.condition_id}_run{rep:02d}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            if not np.array_equal(saved['X'], X):
                raise ValueError(f'Existing observations differ: {path}')
        return
    edges = np.array(metadata['true_edges'], dtype=int).reshape(-1, 2)
    pairs = np.array(metadata['true_interactions'], dtype=int).reshape(-1, 3)
    A = np.zeros((X.shape[1], X.shape[1]), dtype=np.int8)
    for k, j in edges:
        A[k, j] = 1
    tmp = path.with_suffix('.npz.tmp')
    with tmp.open('wb') as handle:
        np.savez_compressed(handle, X=X, true_edges=edges, true_interactions=pairs, adjacency=A, noise_standard_deviations=metadata['noise_standard_deviations'], metadata_json=np.array(json.dumps(metadata)), sha256_float64=np.array(hashlib.sha256(np.ascontiguousarray(X, dtype='<f8').tobytes()).hexdigest()))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)

def dispose(master):
    if getattr(master, 'backend', None) == 'gurobi':
        master.model.dispose()

def process(condition, rep, args, out, run_id, cache, metadata_keys):
    cid = condition.condition_id
    prefix = f'{cid}/{rep}'
    event = out / 'timing_events.jsonl'

    def timed(stage):
        return phase(event, run_id, cid, rep, stage)
    with timed('generation'):
        X, metadata = simulate(condition, rep, args.seed, args.noise_sd)
        save_dataset(out, condition, rep, X, metadata)
        key = f'{cid}/{rep}'
        if key not in metadata_keys:
            append(out / 'datasets.jsonl', dict(metadata, key=key))
            metadata_keys.add(key)
    with timed('basis'):
        basis, R = generate_spline_basis(X, n_knots=args.n_knots, degree=args.degree, include_intercept=False, normalize_columns=True)
    n, num_nodes = X.shape
    edges = [(k, j) for k in range(num_nodes) for j in range(num_nodes) if k != j]
    pairs = all_pairs(edges)
    lambdas = list(dict.fromkeys(args.lambda_values or [2 * num_nodes * np.log(num_nodes) ** 2 / n * v for v in args.lambda_multipliers]))
    grid = list(dict.fromkeys(((float(le), float(rho)) for rho in args.rho_values for le in lambdas)))

    def master(es, ws=None, relaxed=False):
        m = Master(X, basis, es, ws, forced=False, relaxed=relaxed, gamma_main=args.gamma_main, gamma_interaction=args.gamma_interaction, backend=args.backend)
        if args.backend == 'gurobi':
            m.model.Params.Threads = args.threads
        return m

    def run(m, key, le, rho, relaxed=False, start=None):
        if key in cache:
            return dict(cache.get(key), cache_hit=True)
        point = m.solve(le, rho * le, args.relax_time_limit if relaxed else args.time_limit, args.relax_gap if relaxed else args.gap, initial_support=start)
        point.update(key=key, lambda_edge=le, lambda_interaction=rho * le, rho=rho, cache_hit=False)
        cache.put(key, point)
        print(f"{key}: {point['status']}, gap={point['gap']:.3g}, {point['seconds']:.1f}s", flush=True)
        return point

    def collect(point, approach, threshold, kept_e, kept_w):
        se = sorted((e for e, v in point['g'].items() if v > 0.5))
        sw = sorted((a for a, v in point['w'].items() if v > 0.5))
        row = {k: v for k, v in point.items() if k not in ('g', 'w')}
        row.update(condition.metadata(), configuration=cid, replicate=rep, approach=approach, pipeline_version=PIPELINE_VERSION, dgp_version=VERSION, threshold=threshold, R=R, true_edge_count=len(metadata['true_edges']), true_interaction_count=len(metadata['true_interactions']), selected_edges=json.dumps(np.asarray(se, dtype=int).reshape(-1, 2).tolist()), selected_interactions=json.dumps(np.asarray(sw, dtype=int).reshape(-1, 3).tolist()), estimated_edges=len(se), estimated_interactions=len(sw), screened_edges=len(kept_e), screened_interactions=len(kept_w), lambda_grid_min=min(lambdas), lambda_grid_max=max(lambdas), lambda_at_grid_boundary=point['lambda_edge'] in (min(lambdas), max(lambdas)))
        row.update(score_graph(X, basis, se, sw, args.gamma_main, args.gamma_interaction))
        truth = np.zeros((num_nodes, num_nodes), dtype=int)
        estimated = truth.copy()
        for k, j in metadata['true_edges']:
            truth[k, j] = 1
        for k, j in se:
            estimated[k, j] = 1
        row['sid'] = structural_intervention_distance(truth, estimated)
        row.update(recovery(metadata, se, sw, kept_e, kept_w))
        append(out / 'rows_checkpoint.jsonl', row)
    for approach in METHODS:
        method_pairs = [] if approach == 'main_effects_only' else pairs
        method_grid = [(float(le), 0.0) for le in lambdas] if approach == 'main_effects_only' else grid
        points = []
        with timed(approach + '_screening'):
            relaxation = master(edges, method_pairs, relaxed=True)
            try:
                for le, rho in method_grid:
                    key = f'{prefix}/{approach}/relax/{le:.17g}/{rho:.17g}'
                    points.append(run(relaxation, key, le, rho, True))
            finally:
                dispose(relaxation)
        by_threshold = {t: [] for t in dict.fromkeys(args.thresholds)}
        with timed(approach + '_provisional_refit'):
            for point in points:
                le, rho = (point['lambda_edge'], point['rho'])
                for threshold, candidates in by_threshold.items():
                    kept_e, kept_w = screen_support(edges, method_pairs, [point], threshold)
                    dag_e, dag_w = provisional_dag(kept_e, kept_w, point['g'])
                    final_key = f'{prefix}/{approach}/screened/{threshold:.17g}/{le:.17g}/{rho:.17g}'
                    candidate = dict(condition.metadata(), key=final_key, final_key=final_key, approach=approach, pipeline_version=PIPELINE_VERSION, replicate=rep, threshold=threshold, lambda_edge=le, lambda_interaction=rho * le, rho=rho, R=R, relaxation_status=point['status'], relaxation_gap=point['gap'], relaxation_seconds=point['seconds'], screened_edges=len(kept_e), screened_interactions=len(kept_w), provisional_edge_count=len(dag_e), provisional_interaction_count=len(dag_w), provisional_edges=json.dumps(dag_e), provisional_interactions=json.dumps(dag_w), screened_edge_support=json.dumps(kept_e), screened_interaction_support=json.dumps(kept_w), provisional_shd=recovery(metadata, dag_e, dag_w, kept_e, kept_w)['shd'])
                    candidate.update(score_graph(X, basis, dag_e, dag_w, args.gamma_main, args.gamma_interaction, only_nominal=False))
                    if metadata['third_parent_edge']:
                        candidate['third_parent_screen_retained'] = int(tuple(metadata['third_parent_edge']) in kept_e)
                    candidates.append(candidate)
            for threshold, candidates in by_threshold.items():
                shortlist = rank_candidates(candidates, args.shortlist_size)
                for candidate in candidates:
                    append(out / 'screening_checkpoint.jsonl', candidate)
                if not shortlist:
                    raise ValueError(f'{cid}/{approach}: no finite provisional BIC model for either criterion. Inspect screening.csv; use larger penalties or a higher threshold.')
        del points
        restricted = {}
        with timed(approach + '_final_mip'):
            try:
                for threshold, candidates in by_threshold.items():
                    finalists = [c for c in candidates if c['shortlisted']]
                    print(f'{cid}/{rep}/{approach}: threshold={threshold:g}, {len(finalists)}/{len(method_grid)} penalty points shortlisted', flush=True)
                    for candidate in finalists:
                        es = [tuple(e) for e in json.loads(candidate['screened_edge_support'])]
                        ws = [tuple(a) for a in json.loads(candidate['screened_interaction_support'])]
                        ds = [tuple(e) for e in json.loads(candidate['provisional_edges'])]
                        di = [tuple(a) for a in json.loads(candidate['provisional_interactions'])]
                        signature = (tuple(es), tuple(ws))
                        if signature not in restricted:
                            restricted[signature] = master(es, ws)
                        point = run(restricted[signature], candidate['final_key'], candidate['lambda_edge'], candidate['rho'], start=(ds, di))
                        point = dict(point, shortlist_rank=candidate['shortlist_rank'], provisional_ols_bic=candidate['ols_bic'], provisional_ridge_bic=candidate['ridge_bic'], shortlisted_by_nominal=candidate['shortlisted_by_nominal'], shortlisted_by_effective_df=candidate['shortlisted_by_effective_df'], nominal_shortlist_rank=candidate['nominal_shortlist_rank'], effective_df_shortlist_rank=candidate['effective_df_shortlist_rank'], relaxation_status=candidate['relaxation_status'], relaxation_gap=candidate['relaxation_gap'], relaxation_seconds=candidate['relaxation_seconds'])
                        collect(point, approach, threshold, es, ws)
            finally:
                for m in restricted.values():
                    dispose(m)

def main():

    def stop_for_checkpoint(*_):
        raise KeyboardInterrupt('Termination requested; checkpointing')
    signal.signal(signal.SIGTERM, stop_for_checkpoint)
    cli = parser()
    args = cli.parse_args()
    if args.summarize_results:
        out = args.output_dir or ROOT / 'results/reaggregated/interaction'
    elif args.save_data_only:
        out = args.output_dir or ROOT / 'data/generated/interaction'
    else:
        out = args.output_dir or ROOT / 'results/rerun/interaction/mip'
    out = out.expanduser().resolve()
    protected = [ROOT / 'data/simulated/interaction', ROOT / 'results/interaction']
    if any(out == path or path in out.parents for path in protected):
        cli.error('Choose an output directory outside the supplied interaction data and results.')
    if args.summarize_results:
        source = args.summarize_results.expanduser().resolve()
        if out == source:
            cli.error('Choose an output directory separate from the results being summarized.')
        summarize_saved_results(source, out)
        return
    if args.export_only:
        if not out.exists():
            cli.error('output directory does not exist')
        export(out)
        print('CSV exports refreshed:', out.resolve())
        return
    if any((not np.isfinite(t) or t < 0 for t in args.noise_taus)):
        cli.error('noise taus must be finite and nonnegative')
    if args.replicates < 1 or args.shortlist_size < 1:
        cli.error('replicates >=1 and shortlist-size >=1 required')
    if args.seed < 0 or args.threads < 1 or args.n_knots < 1 or (args.degree < 1) or (args.export_every < 1):
        cli.error('seed >=0 and positive threads/basis settings required')
    if any((not np.isfinite(v) or v <= 0 for v in [args.noise_sd, args.gamma_main, args.gamma_interaction, args.time_limit, args.relax_time_limit, args.gap, args.relax_gap, *args.lambda_multipliers])):
        cli.error('strengths and scale/budget parameters must be finite and positive')
    if any((not np.isfinite(v) or v < 0 for v in args.rho_values + (args.lambda_values or []))):
        cli.error('penalty grids must be finite and nonnegative')
    if any((not 0 < t <= 1 for t in args.thresholds)):
        cli.error('thresholds must be in (0,1]')
    try:
        conditions = make_conditions(args)
    except ValueError as error:
        cli.error(str(error))
    try:
        if args.settings:
            if not set(args.settings) <= set(range(13)):
                raise ValueError('settings must be between 0 and 12')
            conditions = [c for c in conditions if c.setting_number in args.settings]
    except ValueError as error:
        cli.error(str(error))
    if not conditions:
        cli.error('no conditions selected')

    def run_reps(c):
        return range(args.replicates) if args.trial is None else (args.trial,)
    dataset_count = sum((len(run_reps(c)) for c in conditions))
    if args.list_conditions:
        for c in conditions:
            print(json.dumps(c.metadata()))
        print(len(conditions), 'conditions')
        return
    if args.save_data_only:
        out.mkdir(parents=True, exist_ok=True)
        (out / 'generation_design.json').write_text(json.dumps(dict(seed=args.seed, noise_sd=args.noise_sd, dgp_version=VERSION, replicates=args.replicates, trial=args.trial, conditions=[c.metadata() for c in conditions]), indent=2))
        for c in conditions:
            for rep in run_reps(c):
                X, meta = simulate(c, rep, args.seed, args.noise_sd)
                save_dataset(out, c, rep, X, meta)
        print('Saved data and ground truth:', out.resolve())
        return
    if args.validate_only:
        for c in conditions:
            rep = 0 if args.trial is None else args.trial
            X, meta = simulate(c, rep, args.seed, args.noise_sd)
            assert np.isfinite(X).all() and len(meta['true_interactions']) == c.interaction_count
            print(c.condition_id, 'shape', X.shape, 'true edges', len(meta['true_edges']), 'true interactions', len(meta['true_interactions']), 'candidate edges', c.p * (c.p - 1), 'candidate interactions', c.p * (c.p - 1) * (c.p - 2) // 2)
        print('Validated', len(conditions), 'conditions')
        return
    out.mkdir(parents=True, exist_ok=True)
    excluded = {'replicates', 'trial', 'shortlist_size', 'blocks', 'sample_sizes', 'base_p', 'dimensions', 'interaction_counts', 'interaction_strengths', 'weak_strengths', 'third_parent_strengths', 'settings', 'condition_ids', 'output_dir', 'list_conditions', 'validate_only', 'export_only', 'export_every', 'run_part', 'save_data_only'}
    settings = {k: v for k, v in vars(args).items() if k not in excluded}
    settings['dgp_version'] = VERSION
    settings['pipeline_version'] = PIPELINE_VERSION
    design_path = out / 'design.json'
    if design_path.exists() and json.loads(design_path.read_text())['settings'] != settings:
        cli.error('data seed, fitting or screening settings changed: use a new output directory')
    design_path.write_text(json.dumps(dict(settings=settings, replicates=args.replicates, trial=args.trial, shortlist_size=args.shortlist_size, conditions=[c.metadata() for c in conditions]), indent=2))
    atomic_csv(out / 'conditions.csv', [c.metadata() for c in conditions])
    cache = PointCache(out)
    metadata_keys = {r['key'] for r in read_jsonl(out / 'datasets.jsonl')}
    completed_keys = {r['key'] for r in read_jsonl(out / 'completed_datasets.jsonl')}
    run_id = uuid.uuid4().hex
    start = time.perf_counter()
    status = 'complete'
    print(f'{len(conditions)} conditions, {dataset_count} scheduled datasets; shortlist-size={args.shortlist_size} PER BIC AND METHOD; both methods screened', flush=True)
    try:
        try:
            from threadpoolctl import threadpool_limits
            limiter = threadpool_limits(limits=args.threads)
        except ImportError:
            limiter = None
        completed = 0
        for c in conditions:
            for rep in run_reps(c):
                dataset_key = f'{c.condition_id}/{rep}'
                if dataset_key in completed_keys:
                    print(f'Already complete: {dataset_key}', flush=True)
                    continue
                print(f'Starting {c.condition_id}, replicate {rep}', flush=True)
                with phase(out / 'timing_events.jsonl', run_id, c.condition_id, rep, 'dataset_wall'):
                    process(c, rep, args, out, run_id, cache, metadata_keys)
                    cache.flush()
                    completed += 1
                    append(out / 'completed_datasets.jsonl', dict(key=dataset_key))
                    completed_keys.add(dataset_key)
                    if completed % args.export_every == 0:
                        with phase(out / 'timing_events.jsonl', run_id, c.condition_id, rep, 'export'):
                            export(out)
    except BaseException:
        status = 'interrupted'
        raise
    finally:
        cache.flush()
        export(out)
        append(out / 'runs.jsonl', dict(run_id=run_id, status=status, seconds=time.perf_counter() - start, conditions=len(conditions), replicates=args.replicates, shortlist_size=args.shortlist_size))
        if 'limiter' in locals() and limiter is not None:
            limiter.restore_original_limits()
    print('Saved', out.resolve(), 'elapsed seconds', round(time.perf_counter() - start, 2))

def structural_intervention_distance(true_adjacency, estimated_adjacency):
    """Count incorrect ordered intervention pairs for two DAGs.

    Implements Lemma 1 and Proposition 2 of Peters and Buehlmann (2015),
    Structural Intervention Distance for Evaluating Causal Graphs.
    A[parent, child] = 1. All simple skeleton paths are enumerated, so this
    implementation is intended for small graphs such as the eight-node study.
    """
    G = np.asarray(true_adjacency, dtype=bool)
    H = np.asarray(estimated_adjacency, dtype=bool)
    if G.ndim != 2 or G.shape[0] != G.shape[1] or H.shape != G.shape:
        raise ValueError('Expected equally sized square adjacency matrices')
    p = len(G)

    def reachability(A):
        R = A.copy()
        for k in range(p):
            R |= R[:, k, None] & R[None, k, :]
        if np.diag(R).any():
            raise ValueError('SID requires acyclic graphs')
        return R
    reach = reachability(G)
    reachability(H)
    descendants = reach | np.eye(p, dtype=bool)
    neighbors = [np.flatnonzero(G[i] | G[:, i]).tolist() for i in range(p)]
    mistakes = 0
    for i in range(p):
        Z = set(np.flatnonzero(H[:, i]))
        collider_open = descendants[:, list(Z)].any(axis=1) if Z else np.zeros(p, bool)
        for j in range(p):
            if i == j:
                continue
            if j in Z:
                mistakes += int(reach[i, j])
                continue
            causal_nodes = np.flatnonzero(reach[i] & descendants[:, j])
            forbidden = descendants[causal_nodes].any(axis=0) if len(causal_nodes) else np.zeros(p, bool)
            if any((forbidden[z] for z in Z)):
                mistakes += 1
                continue

            def open_nondirected_path(path):
                last = path[-1]
                if last == j:
                    if all((G[a, b] for a, b in zip(path, path[1:]))):
                        return False
                    for a, b, c in zip(path, path[1:], path[2:]):
                        collider = G[a, b] and G[c, b]
                        if collider and (not collider_open[b]) or (not collider and b in Z):
                            return False
                    return True
                return any((open_nondirected_path(path + [n]) for n in neighbors[last] if n not in path))
            mistakes += int(open_nondirected_path([i]))
    return mistakes

def fit_dag(X, *, interactions=True, lambda_edge=0.05, rho=1.0, gamma_main=1.0, gamma_interaction=1.0, threshold=0.3, relax_time_limit=120.0, time_limit=300.0, gap=0.001, threads=8, backend='gurobi'):
    """Fit one penalty configuration, with continuous screening and an integer solve.

    Returns zero-based edge and interaction lists, adjacency, optimization
    diagnostics, BIC scores, and separate screening/solve times. Use the CLI
    to search and shortlist the full penalty grid instead.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or not np.isfinite(X).all():
        raise ValueError('X must be a finite observations-by-variables matrix')
    basis, _ = generate_spline_basis(X)
    edges = [(k, j) for k in range(X.shape[1]) for j in range(X.shape[1]) if k != j]
    pairs = all_pairs(edges) if interactions else []
    options = dict(gamma_main=gamma_main, gamma_interaction=gamma_interaction, backend=backend)
    relaxation = Master(X, basis, edges, pairs, relaxed=True, **options)
    try:
        if backend == 'gurobi':
            relaxation.model.Params.Threads = threads
        point = relaxation.solve(lambda_edge, rho * lambda_edge, relax_time_limit, gap)
    finally:
        dispose(relaxation)
    es, ws = screen_support(edges, pairs, [point], threshold)
    start = provisional_dag(es, ws, point['g'])
    integer = Master(X, basis, es, ws, **options)
    try:
        if backend == 'gurobi':
            integer.model.Params.Threads = threads
        fitted = integer.solve(lambda_edge, rho * lambda_edge, time_limit, gap, initial_support=start)
    finally:
        dispose(integer)
    se = sorted((e for e, value in fitted['g'].items() if value > 0.5))
    sw = sorted((a for a, value in fitted['w'].items() if value > 0.5))
    A = np.zeros((X.shape[1], X.shape[1]), dtype=int)
    for k, j in se:
        A[k, j] = 1
    return dict(edges=np.asarray(se, int).reshape(-1, 2).tolist(), interactions=np.asarray(sw, int).reshape(-1, 3).tolist(), adjacency=A, screening_seconds=point['seconds'], solve_seconds=fitted['seconds'], gap=fitted['gap'], status=fitted['status'], **score_graph(X, basis, se, sw, gamma_main, gamma_interaction))

def generate_dataset(setting=0, run=0, noise_tau=0.1, seed=1):
    """Generate one dataset; setting (0..12) and run (0..19) are zero-based.

    Returns X, a true adjacency matrix, zero-based true interactions and metadata.
    The low-level metadata uses zero-based graph indices, as in the solver.
    """
    if setting not in range(13) or run < 0:
        raise ValueError('Expected setting in 0..12 and a nonnegative run number')
    args = argparse.Namespace(noise_taus=[noise_tau], condition_ids=None)
    condition = next((c for c in make_conditions(args) if c.setting_number == setting))
    X, metadata = simulate(condition, run, seed)
    A = np.zeros((8, 8), dtype=int)
    for k, j in metadata['true_edges']:
        A[k, j] = 1
    return dict(X=X, adjacency=A, interactions=np.asarray(metadata['true_interactions'], int).reshape(-1, 3).tolist(), metadata=metadata)

def summarize_saved_results(directory, output_dir):
    """Recompute graph metrics and aggregates from the portable results package."""
    directory = Path(directory)
    models = pd.read_csv(directory / 'mip_results.csv')
    with np.load(directory / 'mip_dags.npz', allow_pickle=False) as graphs:
        estimated = graphs['estimated_adjacency']
        truth = graphs['true_adjacency']
    if len(models) != len(estimated) or estimated.shape != truth.shape:
        raise ValueError('Model table and graph arrays do not match')
    indices = models.graph_number.to_numpy(dtype=int)
    if set(indices) != set(range(len(models))):
        raise ValueError('graph_number must enumerate every graph exactly once')
    models['shd'] = [int(np.count_nonzero(truth[i] != estimated[i])) for i in indices]
    models['sid'] = [structural_intervention_distance(truth[i], estimated[i]) for i in indices]
    records = []
    keys = ['setting', 'noise_tau', 'setting_label', 'method', 'bic']
    for values, group in models.groupby(keys, sort=False):
        row = dict(zip(keys, values))
        row['completed_runs'] = len(group)
        for col in ['shd', 'sid', 'total_method_seconds']:
            row['mean_' + col] = group[col].mean()
            row['se_' + col] = group[col].std(ddof=1) / np.sqrt(len(group))
        row['time_limit_rate'] = (group.status == 'time_limit').mean()
        records.append(row)
    output_dir.mkdir(parents=True, exist_ok=True)
    models.to_csv(output_dir / 'mip_results.csv', index=False)
    pd.DataFrame(records).to_csv(output_dir / 'mip_summary.csv', index=False)
    print(f'Summarized {len(models)} models into {len(records)} groups')
if __name__ == '__main__':
    main()
