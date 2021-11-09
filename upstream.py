import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import numpy as np

from build import nkw, F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, Dv, dth, in_nodes
from config import ks, v, R, cs



def solve_equation(matrix, sresult):
    return sprlin.spsolve(matrix, sresult)



def update_matrix(vnow, pnow, oxresult, reg_reg_edges, reg_something_edges, other_edges):
    #sprawdzić zmiany ks
    data, row, col = [], [], []
    sresult = np.zeros(nkw)
    diag = np.ones(nkw) * ks
    for n1, n2, d, l in reg_reg_edges+reg_something_edges+other_edges:
        if oxresult[n1] == 1:
            if oxresult[n2] == 1:
                if pnow[n1] > pnow[n2]:
                    res = v
                    data.append(res)
                    row.append(n1)
                    col.append(n2)
                    diag[n2] -= res
                elif pnow[n2] > pnow[n1]:
                    res = v
                    data.append(res)
                    row.append(n2)
                    col.append(n1)
                    diag[n1] -= res
            elif oxresult[n2] == 0:
                sresult[n1] += R * vnow[n2]
        elif oxresult[n2] == 1:
            sresult[n2] += R * vnow[n1]

    
    row2 = np.array(row)
    col2 = np.array(col)    
    removetab = []
    for node, datum in enumerate(diag):
        
        #usuwamy z przepływu sygnału węzły, do których sygnał tylko wpływa, a nie wypływa (inaczej powstają osobliwosci) - poza węzłami wejsciowymi, w nich wartosci zerujemy w mainie
        if datum == ks and node not in in_nodes:
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
            

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw)), sresult




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
    return (z0_ox-1/(1+np.exp(F1_ox*(F-F0_ox)))) * dt_ox
    

def update_graph(snow, pnow, reg_reg_edges, reg_something_edges, in_edges):
    
    #dodać sigmoidy zamiast liniowych wzrostów
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        if pnow[n1] > pnow[n2]:
            d += cs * snow[n1]
        else:
            d += cs * snow[n2]
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        if pnow[n1] > pnow[n2]:
            d += cs * snow[n1]
        else:
            d += cs * snow[n2]
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        if pnow[n1] > pnow[n2]:
            d += cs * snow[n1]
        else:
            d += cs * snow[n2]
        in_edges[i] = (n1, n2, d, l)

    return reg_reg_edges, reg_something_edges, in_edges

