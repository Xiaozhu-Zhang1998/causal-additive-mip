"""Sufficient-statistic vanilla formulation; supplement B and Table 2.

"""
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def Mips_DAG_diffvar(Sigma_hat, p, lam, gamma, R, moral, time_limit=300, known_edges=[], partial_order=[], tau=0.0, verbose=1):
    list_edges = []
    for edge in moral:
        list_edges.append((edge[0], edge[1]))
        list_edges.append((edge[1], edge[0]))
    m = gp.Model()
    Gamma = {}
    for i in range(p):
        for j in range(p):
            if i == j:
                Gamma[i, j] = m.addVar(lb=1e-05, vtype=GRB.CONTINUOUS, name='Gamma_%s_%s_%s' % (i, j, 0))
            else:
                Gamma[i, j] = 0
    for r in range(1, R + 1):
        for j in range(p):
            for i in range(p):
                if (i, j) in list_edges:
                    Gamma[j, i + p * r] = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='Gamma_%s_%s_%s' % (i, j, r))
                else:
                    Gamma[j, i + p * r] = 0
    psi = m.addMVar((p, 1), lb=1, ub=p, vtype=GRB.CONTINUOUS, name='psi')
    g = {}
    for i, j in list_edges:
        if (i, j) in known_edges:
            g[i, j] = 1
        else:
            g[i, j] = m.addVar(vtype=GRB.BINARY, name='g_%s_%s' % (i, j))
    T = {}
    for i in range(p):
        T[i] = m.addVar(lb=-10, ub=100, vtype=GRB.CONTINUOUS, name='T_%s' % i)
    for i in range(p):
        m.addGenConstrLog(Gamma[i, i], T[i])
    log_term = gp.LinExpr()
    for i in range(p):
        log_term += -2 * T[i]
    trace = gp.QuadExpr()
    for k in range((R + 1) * p):
        for i in range((R + 1) * p):
            for j in range(p):
                trace += Gamma[j, k] * Gamma[j, i] * Sigma_hat[i, k]
    penalty = gp.LinExpr()
    for i, j in list_edges:
        penalty += lam / 2 * g[i, j]
    ridge_bar = gp.QuadExpr()
    for j in range(p):
        for k in range(p, (R + 1) * p):
            ridge_bar += Gamma[j, k] * Gamma[j, k]
    ridge_bar *= 1 / gamma
    m.setObjective(log_term + trace + ridge_bar + penalty, GRB.MINIMIZE)
    m.Params.OutputFlag = verbose
    m.Params.TimeLimit = time_limit
    m.update()
    m.optimize()
    big_M = 0
    for i in range(p):
        big_M = max(big_M, abs(Gamma[i, i].x))
    for r in range(1, R + 1):
        for i, j in list_edges:
            big_M = max(big_M, abs(Gamma[j, i + p * r].x))
    M = 2 * big_M
    m.addConstrs((Gamma[i, i] <= M for i in range(p)))
    for r in range(1, R + 1):
        for i, j in list_edges:
            m.addConstr(M * g[i, j] >= Gamma[j, i + p * r])
            m.addConstr(-M * g[i, j] <= Gamma[j, i + p * r])
    for i, j in list_edges:
        m.addConstr(1 - p + p * g[i, j] <= psi[j] - psi[i])
    for i, j in partial_order:
        m.addConstr(psi[i] <= psi[j])
    m.setObjective(log_term + trace + ridge_bar + penalty, GRB.MINIMIZE)
    m.Params.OutputFlag = verbose
    m.Params.TimeLimit = time_limit
    if tau > 0:
        m.Params.MIPGapAbs = tau
    m.update()
    m.optimize()
    Gamma_val = np.zeros((p, p * (R + 1)))
    for var in m.getVars():
        if 'Gamma' in var.varName:
            name = var.varName.split('_')[1:]
            i = int(name[0])
            j = int(name[1])
            r = int(name[2])
            Gamma_val[j, i + p * r] = var.X
    g_val = []
    for var in m.getVars():
        if 'g' in var.varName and var.X == 1:
            name = var.varName.split('_')[1:]
            g_val.append((int(name[0]), int(name[1])))
    for i, j in list_edges:
        if (i, j) in known_edges:
            g_val.append((i, j))
    psi_val = {}
    for var in m.getVars():
        if 'psi' in var.varName:
            name = var.varName.split('[')[1].split(']')[0].split(',')
            psi_val[int(name[0])] = var.X
    varhat = 1 / np.diag(Gamma_val) ** 2
    Beta_val = Gamma_val * np.sqrt(varhat)[:, np.newaxis]
    loss_val = np.sum(np.log(varhat)) + p
    return (g_val, m.ObjVal)
