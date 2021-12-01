import numpy as np
import scipy.sparse as spr

from config import simInputData

def create_vector(sid:simInputData, in_nodes_ox):
    bloodoxresult = np.zeros(sid.nsq)
    for node in in_nodes_ox:
        bloodoxresult[node] = 1
    return bloodoxresult


def update_matrix(sid:simInputData, pnow, oxresult, bloodoxresult, edges):
    data, row, col = [], [], []
    diag = np.ones(sid.nsq) * sid.kb
    for n1, n2, d, l, t in edges:
        if bloodoxresult[n1] == 1:
            diag[n1] = 1
            if bloodoxresult[n2] == 1:
                diag[n2] = 1
            elif oxresult[n2] == 1:
                if pnow[n2] < pnow[n1]:
                    res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l)
                    data.append(res)
                    row.append(n2)
                    col.append(n1)
                    diag[n2] -= res
        elif bloodoxresult[n2] == 1:
            diag[n2] = 1
            if oxresult[n1] == 1:
                if pnow[n1] < pnow[n2]:
                    res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                    data.append(res)
                    row.append(n1)
                    col.append(n2)
                    diag[n1] -= res
        elif oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] < pnow[n2]:
                res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n1] -= res
            elif pnow[n2] < pnow[n1]:
                res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res

    
    row2 = np.array(row)
    col2 = np.array(col)    
    removetab = []
    for node, datum in enumerate(diag):
        
        #usuwamy z przepływu sygnału węzły, do których sygnał tylko wpływa, a nie wypływa (inaczej powstają osobliwosci) - poza węzłami wejsciowymi, w nich wartosci zerujemy w mainie
        if datum == sid.kb:
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
            
    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))