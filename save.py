# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
import json

                                                                            
def save_all(name, reg_reg_edges, reg_something_edges, other_edges, oxresult,
             n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters,
             F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox,
             in_nodes, out_nodes, reg_nodes, other_nodes,in_nodes_ox, out_nodes_ox,in_edges,
             G, boundary_nodes_out, boundary_nodes_in, boundary_edges):

    consts = [n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
              F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, iters]
    pos = nx.get_node_attributes(G,'pos')
    
    new_pos = {}
    for x in pos:
        new_pos[x] = tuple(pos[x])

    lists = [in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult.tolist()]
    tuple_lists = [in_edges, reg_reg_edges, reg_something_edges, other_edges, boundary_edges]
    All = [consts, lists, tuple_lists]
    S = json.dumps(All)
    f = open(name+".json","w")
    f.write(S)
    f.close()
    
    S = json.dumps(new_pos)
    f = open(name+"2.json","w")
    f.write(S)
    f.close()    


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


    (n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
     F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters) = tuple(consts)

    in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult = tuple(lists)
    oxresult = np.array(oxresult)
    in_edges, reg_reg_edges, reg_something_edges, other_edges, boundary_edges = tuple(tuple_lists)
    
    in_edges = array2tuple(in_edges)
    reg_reg_edges = array2tuple(reg_reg_edges)
    reg_something_edges = array2tuple(reg_something_edges)
    other_edges = array2tuple(other_edges)
    boundary_edges = array2tuple(boundary_edges)

    def reproduct():
        G1 = nx.Graph()
        new_pos = {}
        for key, value in pos.items():
            new_pos[int(key)] =  value
        
        for node in new_pos:
            G1.add_node(node, pos = new_pos[node])
        for n1, n2, d, L in reg_reg_edges:
            G1.add_edge(n1,n2, d=d,q=0, length=L)
        for n1, n2, d, L in reg_something_edges:
            G1.add_edge(n1,n2, d=d,q=0, length=L)
        return G1
    
    G1 = reproduct()
    return  (G1, boundary_edges, n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
             F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters,
             in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult,
             in_edges, reg_reg_edges, reg_something_edges, other_edges)
    
def save_config(name, n, geo, nettype, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                        F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, length_wiggle_param, noise):
    consts = [n, geo, nettype, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                        F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, length_wiggle_param, noise]
    constnames = ['n', 'geo', 'nettype', 'F0', 'F1', 'z0', 'z1', 'F_mult', 'dt', 'c1', 'c2', 'l', 'mu', 'qin', 'presout', 'D', 'Dv', 'k', 'dth',
                        'F0_ox', 'F1_ox', 'z0_ox', 'z1_ox', 'F_mult_ox', 'dt_ox', 'length_wiggle_param', 'noise']
    file1 = open(name+".txt","w")
    for i in range(len(consts)):
        file1.write(constnames[i])
        file1.write(" ")
        file1.write(str(consts[i]))
        file1.write("\n")
    file1.close()    
