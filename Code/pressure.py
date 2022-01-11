from collections import defaultdict
import numpy as np
import scipy.sparse as spr

from utils import mu_d, d_update

from config import simInputData


def create_pressureFlowEqRHS(sid:simInputData, in_nodes, out_nodes):
    pressureFlowEqRHS = np.zeros(sid.networkSizeSquared)
    for node in in_nodes:
        pressureFlowEqRHS[node] = -sid.flow
    for node in out_nodes:
        pressureFlowEqRHS[node] = sid.outputPressure
    return pressureFlowEqRHS


def create_pressureFlowEqMatrix(sid:simInputData, edges, in_nodes, out_nodes):

    data, row, col = [], [], []
    insert = defaultdict(float)
    diag = np.zeros(sid.networkSizeSquared)
    for n1, n2, d, l, t in edges:
        res = sid.flowConstC1 / mu_d(d) * d ** 4 / l
        if t == 0:
            data.append(res)
            row.append(n1)
            col.append(n2)
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n1] -= res
            diag[n2] -= res
        elif t == 1:
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            insert[n1] += res
        elif t == 2:
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

    sum_insert = sum(insert.values())

    for node in in_nodes:
        data.append(-sum_insert)
        row.append(node)
        col.append(node)

        for ins_node, ins in insert.items():
            data.append(ins)
            row.append(node)
            col.append(ins_node)

    return spr.csr_matrix((data, (row, col)), shape=(sid.networkSizeSquared, sid.networkSizeSquared))

"""
def update_matrix_new(edges, networkSizeSquared):
    data, row, col = [], [], []
    diag = np.zeros(networkSizeSquared)
    for (n1, n2, d, l, edgeType) in edges:
        res = (np.pi * d ** 4) / (128 * mu_d(d) * l)
        if edgeType == 0:
            data.extend([res, res])
            row.extend([n1, n2])
            col.extend([n2, n1])
            diag[n1] -= res; diag[n2] -= res
"""



def update_graph(sid:simInputData, edges, pnow):
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        F = sid.forceMultiplier / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F, sid.F_p)
        if d < sid.minimalDiameter:
            d = sid.minimalDiameter
        elif d > sid.maximalDiameter:
             d = sid.maximalDiameter
        edges[i] = (n1, n2, d, l, t)

    return edges

def update_graph_gradp(sid:simInputData, edges, pnow):
    pin = np.max(pnow)
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        F = sid.cp * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F, sid.F_gradp)
        if d < sid.minimalDiameter:
            d = sid.minimalDiameter
        elif d > sid.maximalDiameter:
             d = sid.maximalDiameter
        edges[i] = (n1, n2, d, l, t)

    return edges

def update_graph_gradp_tips(sid:simInputData, edges, pnow, oxresult):
    oxresult2 = oxresult.copy()
    pin = np.max(pnow)
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = sid.cp * np.abs(pnow[n1] - pnow[n2]) / l
            d += d_update(F, sid.F_gradp)
        if d < sid.minimalDiameter:
            d = sid.minimalDiameter
        elif d > sid.maximalDiameter:
             d = sid.maximalDiameter
        edges[i] = (n1, n2, d, l, t)
    
    return edges


def update_network(G1, sid:simInputData, edges, pnow):
    Q_in = 0
    Q_out = 0
    
    for n1, n2, d, l, t in edges:
        G1[n1][n2]['d']= d
        q = sid.flowConstC1 / mu_d(d) * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        G1[n1][n2]['q'] = q

        if t == 1:    
            Q_in += q
        if t == 2:
            Q_out += q
    
    print(f'\tFlow of blood INTO the network: {Q_in:.4f}\n\tFlow of blood OUT of the network: {Q_out:.4f}')
   #print('Q_in =', Q_in, 'Q_out =', Q_out)
    return G1
