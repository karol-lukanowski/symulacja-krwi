import numpy as np
import scipy.sparse as spr

from utils import solve_equation

from config import simInputData


def create_vector(sid:simInputData, inc_matrix, oxygen, blood_vessels, in_nodes, out_nodes):
    """ Creates vector result for pressure calculation.
    
    For inlet and outlet nodes elements of the vector correspond explicitly
    to the pressure in nodes, for regular nodes elements of the vector equal
    0 correspond to flow continuity.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        nsq - number of nodes in the network squared

    in_nodes : list
        list of inlet nodes
    
    Returns
    -------
    scipy sparse vector
        result vector for pressure calculation
    """
    inlet = np.zeros(sid.nsq)
    for node in in_nodes + out_nodes:
        inlet[node] = 1
    blood_nodes = 1 * ((np.abs(inc_matrix.transpose()) @ blood_vessels + inlet) != 0)
    vegf_b = (1 - blood_nodes) * sid.vegf_const / (1 + np.exp(-sid.vegf_a * (oxygen - sid.vegf_b)))    
    return vegf_b

def find_vegf(sid, inc_matrix, lens, blood_vessels, vegf_b, in_nodes, out_nodes):

    # we solve equation nabla^2 c_v = c1 sigmoid(c_0) in tissue
    # we set zero concentration of c_v on nodes that are part of blood vessels network
    # so we create full matrices and afterwards zero the rows for blood nodes and give them one on diagonal
    # first we create matrix for diffusion in tissue
    vegf_diff_matrix = inc_matrix.transpose() @ spr.diags((1 - blood_vessels) / lens) @ inc_matrix
    inlet = np.zeros(sid.nsq)
    for node in in_nodes + out_nodes:
        inlet[node] = 1
    blood_nodes = 1 * ((np.abs(inc_matrix.transpose()) @ blood_vessels + inlet) != 0)
    vegf_matrix = spr.diags(1 - blood_nodes) * vegf_diff_matrix + spr.diags(blood_nodes)
    
    vegf = solve_equation(vegf_matrix, vegf_b)
    return vegf






def update_matrix(sid:simInputData, vresult, edges):

    data, row, col = [], [], []

    diag = np.zeros(sid.nsq)
    for n1, n2, d, l, t in edges:
        if vresult[n1] == 0:
            if vresult[n2] == 0:
                diag[n1] = 1
                diag[n2] = 1
            else:
                diag[n1] = 1
                res = sid.Dv / l
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res
        elif vresult[n2] == 0:
            diag[n2] = 1
            res = sid.Dv / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
        else:
            res = sid.Dv / l
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

    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))
    

def update_graph(sid:simInputData, vnow, oxresult, edges):
    d_vegf = np.zeros(len(edges))
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if (oxresult[n1] != 0 or oxresult[n2] != 0):
            F = np.abs(vnow[n1] - vnow[n2]) / l
            d_vegf[i] = d_update(F, sid.F_ox)
    return d_vegf


def update_blood(sid:simInputData, oxresult, edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if d > sid.dth:
            if oxresult2[n1] != 0 and oxresult[n2] == 0:
                oxresult[n2] = oxresult2[n1]
            elif oxresult2[n2] != 0 and oxresult[n1] == 0:
                oxresult[n1] = oxresult2[n2]
    return oxresult