from collections import defaultdict
import numpy as np
import scipy.sparse as spr

from utils import d_update

from config import simInputData


def create_vector(sid:simInputData, in_nodes, out_nodes):
    presult = np.zeros(sid.nsq)
    for node in in_nodes:
        presult[node] = sid.presin
    for node in out_nodes:
        presult[node] = sid.presout
    return presult

def update_matrix(sid:simInputData, edges, in_nodes, out_nodes):

    data, row, col = [], [], []
    insert = defaultdict(float)
    diag = np.zeros(sid.nsq)
    for n1, n2, d, l, t in edges:
        res = sid.c1 / sid.mu * d ** 4 / l
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

    for node in in_nodes + out_nodes:
        row.append(node)
        col.append(node)
        data.append(1)


    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))

def update_graph(sid:simInputData, edges, pnow):
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        F = sid.F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        d += d_update(F, sid.F_p)
        if d < sid.dmin:
            d = sid.dmin
        elif d > sid.dmax:
             d = sid.dmax
        edges[i] = (n1, n2, d, l, t)

    return edges



def update_network(G1, sid:simInputData, edges, pnow):
    Q_in = 0
    Q_out = 0
    
    for n1, n2, d, l, t in edges:
        G1[n1][n2]['d']= d
        q = sid.c1 / sid.mu * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        G1[n1][n2]['q'] = q

        if t == 1:    
            Q_in += q
        if t == 2:
            Q_out += q
    
    print('Q_in =', Q_in, 'Q_out =', Q_out)
    return G1