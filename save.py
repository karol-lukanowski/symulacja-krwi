# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
import json
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
from collections import defaultdict
from config import (
        n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters,  # stałe 
        F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox,                                      # stałe dla oxygen
        in_nodes, out_nodes, reg_nodes, other_nodes,                                       # tablice
        in_nodes_ox, out_nodes_ox,                                                         # tablice dla oxygen
        in_edges,                                                                          # tablica tupli
        G)                                                                                 # sieć
def save_all(name, reg_reg_edges, reg_something_edges, other_edges, oxresult):
    consts = [n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
              F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, iters]
    pos = nx.get_node_attributes(G,'pos')
    new_pos = {}
    for x in pos:
        new_pos[x] = tuple(pos[x])
    #print(new_pos)
    lists = [in_nodes, out_nodes, reg_nodes, other_nodes, in_nodes_ox, out_nodes_ox, oxresult.tolist()]
    tuple_lists = [in_edges, reg_reg_edges, reg_something_edges, other_edges]
    All = [consts, lists, tuple_lists]
    S = json.dumps(All)
    f = open(name+".json","w")
    f.write(S)
    f.close()
    
    S = json.dumps(new_pos)
    f = open(name+"2.json","w")
    f.write(S)
    f.close()    

#    nx.write_gexf(G, name+".gexf")
#    nx.write_gml(G, name+'.gml', stringizer=True)
#    nx.write_adjlist(G, name+".adjlist")

def load(name):
    def array2tuple(tab):
        new_tab = []
        for x in tab:
            x = tuple(x)
            new_tab.append(x)
        return new_tab
    with open(name+".json") as json_file:
        All= json.load(json_file)
    with open(name+"2.json") as json_file:
        pos= json.load(json_file)
    
    consts = All[0]
    lists = All[1]
    tuple_lists = All[2]

#    G = nx.read_gexf(name+".gexf", node_type = int)
#    G = nx.read_gml(name+'.gml', destringizer=True)
#    G = nx.read_adjlist(name+".adjlist")

    (n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
     F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters) = tuple(consts)

    in_nodes, out_nodes, reg_nodes, other_nodes, in_nodes_ox, out_nodes_ox, oxresult = tuple(lists)
    oxresult = np.array(oxresult)
    in_edges, reg_reg_edges, reg_something_edges, other_edges = tuple(tuple_lists)
    
    in_edges = array2tuple(in_edges)
    reg_reg_edges = array2tuple(reg_reg_edges)
    reg_something_edges = array2tuple(reg_something_edges)
    other_edges = array2tuple(other_edges)

    def reproduct():
#        presult = create_vector()
#        pmatrix = update_matrix(reg_reg_edges, reg_something_edges, in_edges)
#        pnow = solve_equation(pmatrix, presult)
        G1 = nx.Graph()
        for node in pos:
            G1.add_node(node, pos = pos[node])
        for n1, n2, d, L in reg_reg_edges:
            G1.add_edge(n1,n2, d=d,q=0, length=L)
        for n1, n2, d, L in reg_something_edges:
            G1.add_edge(n1,n2, d=d,q=0, length=L)
        
        '''    
        for n1, n2, d, L in reg_reg_edges:
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / L
            G.add_edge(n1,n2, d=d,q=q, length=L)
        for n1, n2, d, L in reg_something_edges:        
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / L
            G.add_edge(n1,n2, d=d,q=q, length=L)
        '''
        return G1
    
    G1 = reproduct()
    return  (G1, n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
             F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters,
             in_nodes, out_nodes, reg_nodes, other_nodes, in_nodes_ox, out_nodes_ox, oxresult,
             in_edges, reg_reg_edges, reg_something_edges, other_edges)
    
    
