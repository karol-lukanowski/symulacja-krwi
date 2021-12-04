import networkx as nx
import numpy as np
import os
import scipy.sparse.linalg as sprlin
from delaunay import find_node


def PosGauss(mean, sigma):
    x = np.random.normal(mean, sigma)
    return(x if x>=0 else PosGauss(mean,sigma))
    
def mu_d(d):
    #return mu * (1 + 1 / (15 * d + 0.1))
    #return mu
    d = 100 * d + 0.2
    return ((6*np.exp(-0.085 * d)+2.2-2.44*np.exp(-0.06*d**0.645))*(d/(d-1.1)) ** 2 + 1) * (d/(d-1.1)) ** 2 / 1000

def solve_equation(matrix, result):
    return sprlin.spsolve(matrix, result)

def d_update(F, t):
    if (F > t.F0):
        if (F < t.F1):
            result = t.z0+(F-t.F0)*(t.z1-t.z0)/(t.F1-t.F0)
        else:
            result = t.z1
    else:
        result = t.z0
    return result

def make_dir(sid):
        if not os.path.isdir(sid.dirname):
            os.makedirs(sid.dirname)

        i = 0
        dirname2 = sid.dirname
        while (sid.dirname == dirname2):
            if not os.path.isdir(sid.dirname + "/" + str(i)):
                sid.dirname = sid.dirname + "/" + str(i)
            else:
                i += 1
        os.makedirs(sid.dirname)

class fParams():
    def __init__(self, params):
        self.F0 = params['F0']
        self.F1 = params['F1']
        self.z0 = params['z0']
        self.z1 = params['z1']


def find_ladder(sid, G):
    in1 = find_node(G, [40, 40])
    in2 = find_node(G, [25, 65])
    out1 = find_node(G, [60, 60])
    out2 = find_node(G, [35, 75])
    p1 = nx.shortest_path(G, in1, in2)
    p2 = nx.shortest_path(G, out1, out2)
    for i in range(len(p1) - 1):
        G[p1[i]][p1[i+1]]['d'] = sid.dth
    for i in range(len(p2) - 1):
        G[p2[i]][p2[i+1]]['d'] = sid.dth
    in_nodes = [in1]
    out_nodes = [out1]
    in_nodes_ox = [in1]
    out_nodes_ox = [out1]
    oxresult = np.zeros(sid.nsq)
    for node in p1 + p2:
        oxresult[node] = 1
    return G, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, oxresult

class simAnalysisData:
    pressure = []
    oxygen = []
    vegf = []
    signal = []

def collect_data(sad:simAnalysisData, sid, in_nodes, pnow, vnow, oxnow):
    data = [pnow[in_nodes[0]], np.average(vnow), np.average(oxnow)]
    sad.pressure.append(data[0])
    sad.vegf.append(data[1])
    sad.oxygen.append(data[2])
    f = open(sid.dirname+'/params.txt', 'a')
    f.write(f"\n")
    np.savetxt(f, [[data[0], data[1], data[2]]])
    f.close()