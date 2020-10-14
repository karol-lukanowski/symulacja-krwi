import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import numpy as np
from collections import defaultdict
from build import nkw, F0, F1, z0, z1, F_mult, dt, c1, in_nodes, out_nodes, qin, presout




def solve_equation(matrix, presult):
    return sprlin.spsolve(matrix, presult)

def create_vector():
    global in_nodes, out_nodes
    presult = np.zeros(nkw)
    for node in in_nodes:
        presult[node] = -qin
    for node in out_nodes:
        presult[node] = presout
    return presult

def update_matrix(reg_reg_edges, reg_something_edges, in_edges):

    data, row, col = [], [], []

    diag = np.zeros(nkw)
    for n1, n2, d, l in reg_reg_edges:
        res = c1 * d ** 4 / l
        data.append(res)
        row.append(n1)
        col.append(n2)
        data.append(res)
        row.append(n2)
        col.append(n1)
        diag[n1] -= res
        diag[n2] -= res
    for n1, n2, d, l in reg_something_edges:
        res = c1 * d ** 4 / l
        data.append(res)
        row.append(n1)
        col.append(n2)
        diag[n1] -= res
    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
    for node in out_nodes:
        row.append(node)
        col.append(node)
        data.append(1)

    insert = defaultdict(float)
    for n1, n2, d, l in in_edges:
        insert[n2] += c1 * d ** 4 / l
    sum_insert = sum(insert.values())

    for node in in_nodes:
        data.append(-sum_insert)
        row.append(node)
        col.append(node)

        for ins_node, ins in insert.items():
            data.append(ins)
            row.append(node)
            col.append(ins_node)

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw))

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
#    return (1-1/(1+np.exp(10*(F-0.5))))*dt


def update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges):
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F)
        in_edges[i] = (n1, n2, d, l)

    return reg_reg_edges, reg_something_edges, in_edges

def update_network(G1,reg_reg_edges, reg_something_edges, pnow):
    Q_in = 0
    Q_out = 0
    for n1, n2, d, l in reg_reg_edges:
        G1[n1][n2]['d']= d
        q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        
        G1[n1][n2]['q'] = q

    for n1, n2, d, l in reg_something_edges:        
        G1[n1][n2]['d'] = d
        
        q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        G1[n1][n2]['q'] = q

        if n2 in in_nodes:
            Q_in += q
        if n2 in out_nodes:
            Q_out += q
    
    print('Q_in =', Q_in, 'Q_out =', Q_out)
    return G1