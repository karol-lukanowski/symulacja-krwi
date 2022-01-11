import numpy as np
import scipy.sparse as spr

from config import simInputData
from utils import d_update


def update_matrix_upstream(sid:simInputData, vnow, pnow, edges, in_nodes):
    #sprawdzić zmiany ks
    data, row, col = [], [], []
    sresult = np.zeros(sid.networkSizeSquared)
    diag = np.ones(sid.networkSizeSquared) * sid.ks
    for n1, n2, d, l, t in edges:
        if t == 0:
            if pnow[n1] > pnow[n2]:
                res = sid.v / l
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n2] -= res
            elif pnow[n2] > pnow[n1]:
                res = sid.v / l
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n1] -= res
            sresult[n1] += sid.R * vnow[n2]
            sresult[n2] += sid.R * vnow[n1]
        elif t == 1:
            res = sid.v / l
            diag[n1] -= res

    row2 = np.array(row)
    col2 = np.array(col)    
    removetab = []
    for node, datum in enumerate(diag):
        
        #usuwamy z przepływu sygnału węzły, do których sygnał tylko wpływa, a nie wypływa (inaczej powstają osobliwosci) - poza węzłami wejsciowymi, w nich wartosci zerujemy w mainie
        if datum == sid.ks:
            tab = np.where(row2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)

            tab = np.where(col2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)
            sresult[node] = 0
        
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
            
    removetab = np.flip(np.sort(removetab))
    
    for i in removetab:
        row.pop(i)
        col.pop(i)
        data.pop(i)
            
    return spr.csr_matrix((data, (row, col)), shape=(sid.networkSizeSquared, sid.networkSizeSquared)), sresult

def update_matrix_downstream(sid, vnow, pnow, edges, out_nodes):
    #sprawdzić zmiany ks
    data, row, col = [], [], []
    sresult = np.zeros(sid.networkSizeSquared)
    diag = np.ones(sid.networkSizeSquared) * sid.ks
    for n1, n2, d, l, t in edges:
        if t == 0:
            if pnow[n1] < pnow[n2]:
                res = sid.v / l
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n2] -= res
            elif pnow[n2] < pnow[n1]:
                res = sid.v / l
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n1] -= res
            sresult[n1] += sid.R * vnow[n2]
            sresult[n2] += sid.R * vnow[n1]
        elif t == 1:
            res = sid.v / l
            diag[n1] -= res

    row2 = np.array(row)
    col2 = np.array(col)    
    removetab = []
    for node, datum in enumerate(diag):
        
        #usuwamy z przepływu sygnału węzły, do których sygnał tylko wpływa, a nie wypływa (inaczej powstają osobliwosci) - poza węzłami wejsciowymi, w nich wartosci zerujemy w mainie
        if datum == sid.ks:
            tab = np.where(row2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)

            tab = np.where(col2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)
            sresult[node] = 0
        
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
            
    removetab = np.flip(np.sort(removetab))
    
    for i in removetab:
        row.pop(i)
        col.pop(i)
        data.pop(i)
            
    return spr.csr_matrix((data, (row, col)), shape=(sid.networkSizeSquared, sid.networkSizeSquared)), sresult


    

def update_graph_upstream(sid:simInputData, snow, pnow, oxresult, edges):
    
    #dodać sigmoidy zamiast liniowych wzrostów
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] > pnow[n2]:
                F = sid.cs * snow[n1]
            else:
                F = sid.cs * snow[n2]
            d += d_update(F, sid.F_s)
            # if d > sid.dmax:
            #     d = sid.dmax
        edges[i] = (n1, n2, d, l, t)

    return edges

def update_graph_downstream(sid:simInputData, snow, pnow, oxresult, edges):
    
    #dodać sigmoidy zamiast liniowych wzrostów
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] < pnow[n2]:
                F = sid.cs * snow[n1]
            else:
                F = sid.cs * snow[n2]
            d += d_update(F, sid.F_s)
            # if d > sid.dmax:
            #     d = sid.dmax
        edges[i] = (n1, n2, d, l, t)
    
    return edges