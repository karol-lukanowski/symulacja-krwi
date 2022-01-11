import numpy as np
import scipy.sparse as spr

from utils import d_update

from config import simInputData



def create_vector(sid:simInputData, oxresult, oxnow = []):
    if not sid.includeOxygenEvolution: # if oxygen is off, make oxnow such that vresult = 1
        oxnow = sid.R_c * np.ones(sid.networkSizeSquared)
    vresult = 1 / (1 + np.exp(sid.R_a * (oxnow - sid.R_c)))
    vresult = np.where(oxresult == 1, 0, vresult)
    vresult = -vresult
    return vresult

def update_matrix(sid:simInputData, vresult, edges):

    data, row, col = [], [], []

    diag = np.zeros(sid.networkSizeSquared)
    for n1, n2, d, l, t in edges:
        if vresult[n1] == 0:
            if vresult[n2] == 0:
                diag[n1] = 1
                diag[n2] = 1
            else:
                diag[n1] = 1
                res = sid.VEGFDiffusionConst / l
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res
        elif vresult[n2] == 0:
            diag[n2] = 1
            res = sid.VEGFDiffusionConst / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
        else:
            res = sid.VEGFDiffusionConst / l
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
        else:
            row.append(node)
            col.append(node)
            data.append(1)

    return spr.csr_matrix((data, (row, col)), shape=(sid.networkSizeSquared, sid.networkSizeSquared))
    

def update_graph(sid:simInputData, vnow, oxresult, edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = sid.F_mult_ox * np.abs(vnow[n1] - vnow[n2])/l
            d += d_update(F, sid.F_ox)
            if d > sid.thresholdDiameter:
                oxresult[n1] = 1
                oxresult[n2] = 1
        if d < sid.minimalDiameter:
            d = sid.minimalDiameter
        elif d > sid.maximalDiameter:
             d = sid.maximalDiameter
        edges[i] = (n1, n2, d, l, t)

    return edges


def update_blood(sid:simInputData, oxresult, edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if d > sid.thresholdDiameter:
            if oxresult2[n1] == 1:
                oxresult[n2] = 1
            elif oxresult2[n2] == 1:
                oxresult[n1] = 1
    return oxresult