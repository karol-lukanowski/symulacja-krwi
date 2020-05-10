import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import numpy as np
from config import nkw, F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, D, k, dth, in_nodes, out_nodes



def solve_equation(matrix, oxresult):
    return sprlin.spsolve(matrix, oxresult)

def create_vector(G):
    oxresult = np.zeros(nkw)
    for node in in_nodes:
        oxresult[node] = 1
    for node in out_nodes:
        oxresult[node] = 0
    return oxresult

def update_matrix(oxresult, reg_reg_edges, reg_something_edges):

    data, row, col = [], [], []

    diag = np.ones(nkw) * k
    for n1, n2, d, l in reg_reg_edges:
        if d > dth and oxresult[n2] == 1:
            oxresult[n1] = 1
            diag[n1] = 1
        elif d > dth and oxresult[n1] == 1:
            oxresult[n2] = 1
            diag[n2] = 1
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
            
    for n1, n2, d, l in reg_something_edges:
        if d > dth and oxresult[n2] == 1:
            oxresult[n1] = 1
            diag[n1] = 1
        elif d > dth and oxresult[n1] == 1:
            oxresult[n2] = 1
            diag[n2] = 1
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

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw)), oxresult




def d_update(F):
    #zmiana średnicy pod względem siły F
    result = 0
    if (F > F0_ox):
        if (F < F1_ox):
            result = z0_ox+(F-F0_ox)*(z1_ox-z0_ox)/(F1_ox-F0_ox)
        else:
            result = z1_ox
    else:
        result = z0_ox
    return result * dt_ox

def update_graph(oxnow, oxresult, reg_reg_edges, reg_something_edges, in_edges):


    reg_reg_edges2, reg_something_edges2, in_edges2=[], [], []
    for n1, n2, d, l in reg_reg_edges:
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            F = F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])

            dnew=d+d_update(F)
        else:
            dnew = d
        reg_reg_edges2.append((n1, n2, dnew, l))

    for n1, n2, d, l in reg_something_edges:
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            F=F_mult_ox*np.abs(oxnow[n1] - oxnow[n2])
            dnew = d + d_update(F)
        else:
            dnew = d
        reg_something_edges2.append((n1, n2, dnew, l))
    for n1, n2, d, l in in_edges:
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            F = F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])
            dnew = d + d_update(F)
        else:
            dnew = d
        in_edges2.append((n1, n2, dnew, l))

    return reg_reg_edges2, reg_something_edges2, in_edges2

"""
def update_matrix(G, oxresult):

    matrix = np.zeros((n * n, n * n))

    for node in reg_nodes:
        matrix[node][node] = k
        for ind in G.neighbors(node):
            d = G[node][ind]["d"]
            l = G[node][ind]['length']
            if d > dth:
                matrix[node][node] = 1
                oxresult[node] = 1
                oxresult[ind] = 1
            else:
                matrix[node][ind] = D / l
                matrix[node][node] -= D / l

    for node in out_nodes:
        matrix[node][node] = 1

    for node in in_nodes:
        matrix[node][node] = 1

    return matrix, oxresult

def update_matrix(G, oxresult):

    data, row, col = [], [], []

    for node in reg_nodes:
        this_node = k
        for ind in G.neighbors(node):
            d = G[node][ind]["d"]
            l = G[node][ind]['length']
            if d > dth and oxresult[ind] == 1:
                this_node = 1
                oxresult[node] = 1
#                oxresult[ind] = 1
            else:
                data.append(D/l)
                row.append(node)
                col.append(ind)

                this_node -= D/l

        data.append(this_node)
        row.append(node)
        col.append(node)

    for node in in_nodes + out_nodes:
        data.append(1)
        row.append(node)
        col.append(node)

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw)), oxresult
"""
