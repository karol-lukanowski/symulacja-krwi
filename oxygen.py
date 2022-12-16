import numpy as np
import scipy.sparse as spr

from utils import solve_equation

from config import simInputData


def create_vector(sid:simInputData, in_nodes):
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
    data, row, col = [], [], []
    for node in in_nodes:
        data.append(1)
        row.append(node)
        col.append(0)
    return spr.csc_matrix((data, (row, col)), shape=(sid.nsq, 1))


def find_oxygen(sid, inc_matrix, diams, lens, blood_vessels, flow, in_nodes, out_nodes, oxygen_b):

    # we solve equation d nabla^2 c_o - c1 q / d nabla c_o = c2 c_o where c1 = 0 in tissue (and c2 is different)
    # we set constant concentration of c_o on input as the boundary condition
    # so we create full matrices and afterwards zero the rows for inlet nodes and give them one on diagonal
    # first we create matrix for diffusion, same everywhere
    ox_diff_matrix = inc_matrix.transpose() @ spr.diags((blood_vessels *  diams ** 2 + 1 - blood_vessels) / lens) @ inc_matrix
    # we need a matrix with only downstream elements (and only in blood vessels, for advection terms)
    ox_adv_matrix = np.abs(inc_matrix.transpose() @ (spr.diags(sid.adv_const * flow * blood_vessels) @ inc_matrix > 0))
    ox_adv_matrix.setdiag(0)
    # we create the diagonal for advection terms (equal -inflow for given node), including ones for input nodes
    #diag = -np.abs(inc_matrix.transpose()) @ np.abs(sid.adv_const * flow * blood_vessels) / 2 # find diagonal coefficients (inlet flow for each node)
    diag = -np.array(ox_adv_matrix.sum(axis = 1).flatten())[0]
    for node in in_nodes:
        diag[node] = 1 # set diagonal for input nodes to 1
    for node in out_nodes:
        if diag[node] != 0: # fix for nodes which are connected only to other out_nodes - without it we get a singular matrix (whole row of zeros)
            diag[node] *= 2 # multiply diagonal for output nodes (they have no outlet, so inlet flow is equal to whole flow)
        else:
            diag[node] = 1
    # reaction term for the diagonal
    reaction = sid.rec_bv_const * np.array(np.abs(inc_matrix.transpose()) @ (blood_vessels * diams)) + sid.rec_t_const * np.array(np.abs(inc_matrix.transpose()) @ (1 - blood_vessels))
    inlet = np.zeros(sid.nsq)
    for node in in_nodes:
        inlet[node] = 1
    oxygen_matrix = spr.diags(1 - inlet) * (ox_diff_matrix + ox_adv_matrix + spr.diags(reaction)) + spr.diags(diag)
    
    oxygen = solve_equation(oxygen_matrix, oxygen_b)
    return oxygen





def update_matrix(sid:simInputData, oxresult, bloodoxresult, pnow, edges):
    data, row, col = [], [], []

    diag = np.ones(sid.nsq)
    for n1, n2, d, l, t in edges:
        if bloodoxresult[n1] == 1:
            diag[n1] = 1
        elif oxresult[n1] != 0:
            diag[n1] = np.pi * d * sid.k_ox_b
        else:
            diag[n1] = sid.k_ox_t
        if bloodoxresult[n2] == 1:
            diag[n2] = 1
        elif oxresult[n2] != 0:
            diag[n2] = np.pi * d * sid.k_ox_b
        else:
            diag[n2] = sid.k_ox_t
    for n1, n2, d, l, t in edges:
        if bloodoxresult[n1] == 1:
            if oxresult[n2] == 1:
                if pnow[n2] < pnow[n1]:
                    res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l) #
                    data.append(res)
                    row.append(n2)
                    col.append(n1)
                    diag[n2] -= res
                res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            else:
                res = -sid.D_ox_t / l
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        elif bloodoxresult[n2] == 1:
            if oxresult[n1] == 1:
                if pnow[n1] < pnow[n2]:
                    res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                    data.append(res)
                    row.append(n1)
                    col.append(n2)
                    diag[n1] -= res
                res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            else:
                res = -sid.D_ox_t / l
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
            res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        elif oxresult[n1] == 2 and oxresult[n2] == 2:
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
            res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        else:
            res = -sid.D_ox_t / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res

    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)

    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))


def update_graph(sid:simInputData, oxnow, oxresult, edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = sid.F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])/l
            d += d_update(F, sid.F_ox)
            if d > sid.dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        edges[i] = (n1, n2, d, l)


    return edges, oxresult

def update_oxresult(sid:simInputData, edges, oxresult):
    for e in edges:
        n1, n2, d, l, t = e
        if (oxresult[n1] == 1 and oxresult[n2] == 1):
            if d < sid.dth:
                oxresult[n1] = 0
                oxresult[n2] = 0
        elif (oxresult[n1] == 1 or oxresult[n2] == 1):
            if d > sid.dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        
    return oxresult