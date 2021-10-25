import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import numpy as np
from build import nkw, F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, D, k, dth, in_nodes_ox, out_nodes_ox



def solve_equation(matrix, oxresult):
    return sprlin.spsolve(matrix, oxresult)

def create_vector():
    oxresult = np.zeros(nkw)
    for node in in_nodes_ox:
        oxresult[node] = 1
    #for node in out_nodes_ox:
    #    oxresult[node] = 1
    return oxresult


def update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges):
    data, row, col = [], [], []

    diag = np.ones(nkw) * (-k)
    for n1, n2, d, l in reg_reg_edges+reg_something_edges+other_edges:
        if oxresult[n1] == 1:
            if oxresult[n2] == 1:
                diag[n1] = 1
                diag[n2] = 1
            else:
                diag[n1] = 1
                res = D / l
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res
        elif oxresult[n2] == 1:
            diag[n2] = 1
            res = D / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
        else:
            res = D / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n1] -= res
            diag[n2] -= res

    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw))


def d_update(F):
    #zmiana średnicy pod względem siły F
    '''
    result = 0
    if (F > F0_ox):
        if (F < F1_ox):
            result = z0_ox+(F-F0_ox)*(z1_ox-z0_ox)/(F1_ox-F0_ox)
        else:
            result = z1_ox
    else:
        result = z0_ox
    return result * dt_ox
    '''
    return (1-1/(1+np.exp(30*(F-0.3)))) * dt_ox


def update_graph(oxnow, oxresult, reg_reg_edges, reg_something_edges, in_edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])/l
            d += d_update(F)
            if d > dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        reg_reg_edges[i] = (n1, n2, d, l)


    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F=F_mult_ox*np.abs(oxnow[n1] - oxnow[n2])/l
            d += d_update(F)
            if d > dth:
                oxresult[n1] = 1
                oxresult[n2] = 1            
        reg_something_edges[i] = (n1, n2, d, l)

    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])/l
            d += d_update(F)
            if d > dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        in_edges[i] = (n1, n2, d, l)


    return reg_reg_edges, reg_something_edges, in_edges, oxresult

def update_oxresult(reg_reg_edges, reg_something_edges, in_edges, oxresult):
    for e in reg_reg_edges+reg_something_edges+in_edges:
        n1, n2, d, l = e
        if (oxresult[n1] == 1 and oxresult[n2] == 1):
            if d < dth:
                oxresult[n1] = 0
                oxresult[n2] = 0
        elif (oxresult[n1] == 1 or oxresult[n2] == 1):
            if d > dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        
                
    return oxresult