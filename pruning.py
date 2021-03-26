import pressure as Pr 
from build import F0, F1, z0, z1, F_mult, dt, c1
from config import Fth_prun, dth_prun, z_prun
import numpy as np

def final_pruning(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult):
    Fmax = 0
    for i,e in enumerate(reg_reg_edges + reg_something_edges + in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F > Fmax:
            Fmax = F
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 0.001
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 0.001
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 0.001
        in_edges[i] = (n1, n2, d, l)
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    print (pmatrix)
    pnow = Pr.solve_equation(pmatrix, presult)
    G1 = Pr.update_network(G1, reg_reg_edges, reg_something_edges, pnow)
    return G1


def final_pruning_save(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult):
    qmax = max([edge[2] for edge in G1.edges(data='q')])
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:
                d = 0.001
            reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:
                d = 0.001
            reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:
                d = 0.001
            in_edges[i] = (n1, n2, d, l)
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    print (pmatrix)
    pnow = Pr.solve_equation(pmatrix, presult)
    G1 = Pr.update_network(G1, reg_reg_edges, reg_something_edges, pnow)
    return G1

def d_update(F):
    #zmiana średnicy pod względem siły F
    result = 0
    if (F > F0):
        if (F < F1):
            result = z0+(F-F0)*(z1-z0)/(F1-F0)
        else:
            result = z1
    else:
        result = z0
    return result * dt

def update_graph(G1, pnow, reg_reg_edges, reg_something_edges, in_edges):
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        if d <= 0:
            d = 0.001
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        if d <= 0:
            d = 0.001
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        if d <= 0:
            d = 0.001
        in_edges[i] = (n1, n2, d, l)

    return reg_reg_edges, reg_something_edges, in_edges