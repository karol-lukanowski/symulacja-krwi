import numpy as np
import os

import triangular_net as Tr
import delaunay as De
from geometry import set_geometry, equidistant_geometry, create_edgelist
import save as Sv
from config import load, qdrawconst, ddrawconst, iters, save_every, geo, nettype, save_name, n









if load == 0:
    from config import (F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                        F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, length_wiggle_param, diameter_wiggle_param)
    nkw = n ** 2
    if nettype == 'de':   
        G, boundary_edges = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
        delfix = True
    elif nettype == 'tr':
        G, boundary_edges = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)
        delfix = False
    G2 = G.copy()
    in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in, in_edges = set_geometry(n, G, geo=geo, R=n//2.5, R_s=n//20, **{'del' : delfix})
    reg_reg_edges, reg_something_edges, other_edges = create_edgelist(G, in_nodes, out_nodes, reg_nodes, boundary_nodes_out, boundary_nodes_in)
    #in_nodes, out_nodes = equidistant_geometry(G, n, R = n//2.5, xrange = n, yrange = n, how_many = 200)
    #in_nodes, out_nodes, reg_nodes, other_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)
    
    
    
    in_nodes_ox = in_nodes
    out_nodes_ox = out_nodes
    
    dirname = nettype + geo + "n" + str(n) + "lw" + str(length_wiggle_param) + "dw" + str(diameter_wiggle_param) + "dt" + str(dt) + "dtox" + str(dt_ox)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    old_iters = 0
    
    def create_vector():
        oxresult = np.zeros(nkw)
        for node in in_nodes_ox:
            oxresult[node] = 1
    #    for node in out_nodes_ox:
    #        oxresult[node] = 1
        return oxresult
    oxresult = create_vector()
    Sv.save_all(dirname+'/'+'template', reg_reg_edges, reg_something_edges, other_edges, oxresult,
        n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
        F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
        in_nodes_ox, out_nodes_ox,in_edges,G2, boundary_nodes_out, boundary_nodes_in, boundary_edges)

elif load==1:
    from config import nettype, geo, length_wiggle_param, diameter_wiggle_param, dt, dt_ox, load_name
    dirname = nettype + geo + "n" + str(n) + "lw" + str(length_wiggle_param) + "dw" + str(diameter_wiggle_param) + "dt" + str(dt) + "dtox" + str(dt_ox)
#    dirname = "templatki"
    (G, boundary_edges, n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
     F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters,
     in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in,
     in_nodes_ox, out_nodes_ox, oxresult, in_edges, reg_reg_edges, reg_something_edges, other_edges) = Sv.load(dirname+"/"+load_name)
    nkw = n ** 2

else:
    from config import (F0, F1, z0, z1, F_mult, dt, c1, c2, mu, qin, presout, D, Dv, k, dth,
                    F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, load_name)
    (G, boundary_edges, n, _, _, _, _, _, _, _, _, l, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in,
     in_nodes_ox, out_nodes_ox, oxresult, in_edges, reg_reg_edges, reg_something_edges, other_edges) = Sv.load('templatki/'+load_name)

    dirname = nettype + geo + "n" + str(n) + "w" + load_name + "dt" + str(dt) + "dtox" + str(dt_ox)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    old_iters = 0
    nkw = n**2
    in_nodes_ox = in_nodes
    out_nodes_ox = out_nodes

    def create_vector():
        oxresult = np.zeros(nkw)
        for node in in_nodes_ox:
            oxresult[node] = 1
    #    for node in out_nodes_ox:
    #        oxresult[node] = 1
        return oxresult
    oxresult = create_vector()
    