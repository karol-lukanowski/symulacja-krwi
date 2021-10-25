import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import numpy as np

from build import nkw, F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, Dv, dth, in_nodes, out_nodes, in_nodes_ox
from config import kb, mu

def create_vector():
    bloodoxresult = np.zeros(nkw)
    for node in in_nodes_ox:
        bloodoxresult[node] = 1
    return bloodoxresult

def solve_equation(matrix, bloodoxresult):
    return sprlin.spsolve(matrix, bloodoxresult)

def update_matrix(pnow, oxresult, bloodoxresult, reg_reg_edges, reg_something_edges, other_edges):
    data, row, col = [], [], []
    diag = np.ones(nkw) * kb
    for n1, n2, d, l in reg_reg_edges+reg_something_edges+other_edges:
        if bloodoxresult[n1] == 1:
            diag[n1] = 1
            if bloodoxresult[n2] == 1:
                diag[n2] = 1
            elif oxresult[n2] == 1:
                if pnow[n2] < pnow[n1]:
                    res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * mu * l)
                    data.append(res)
                    row.append(n2)
                    col.append(n1)
                    diag[n2] -= res
        elif bloodoxresult[n2] == 1:
            diag[n2] = 1
            if oxresult[n1] == 1:
                if pnow[n1] < pnow[n2]:
                    res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * mu * l)
                    data.append(res)
                    row.append(n1)
                    col.append(n2)
                    diag[n1] -= res
        elif oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] < pnow[n2]:
                res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * mu * l)
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n1] -= res
            elif pnow[n2] < pnow[n1]:
                res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * mu * l)
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res

    
    row2 = np.array(row)
    col2 = np.array(col)    
    removetab = []
    for node, datum in enumerate(diag):
        
        #usuwamy z przepływu sygnału węzły, do których sygnał tylko wpływa, a nie wypływa (inaczej powstają osobliwosci) - poza węzłami wejsciowymi, w nich wartosci zerujemy w mainie
        if datum == kb:
            tab = np.where(row2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)

            tab = np.where(col2 == node)[0]
            for i in np.flip(tab):
                removetab.append(i)

        
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
            
    removetab = np.flip(np.sort(removetab))

    for i in removetab:
        row.pop(i)
        col.pop(i)
        data.pop(i)
            
    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw))