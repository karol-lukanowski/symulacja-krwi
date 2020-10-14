import numpy as np
import os
import draw_net as Dr
import pressure as Pr
import oxygen as Ox
import vegf as Ve
import save as Sv
from config import save_every
from build import (G, in_nodes, out_nodes, reg_nodes, other_nodes, iters, old_iters,
                   in_edges, boundary_nodes_out, boundary_nodes_in, dirname,
                   n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                   F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes_ox, out_nodes_ox, oxresult,
                   boundary_edges, reg_reg_edges, reg_something_edges, other_edges)



oxcopy = oxresult.copy()
regcopy = reg_reg_edges.copy()
sthcopy = reg_something_edges.copy()
otcopy = other_edges.copy()
incopy = in_edges.copy()


X = 2
k = 0.1


for i in range(X):
    # wyzerowanie sieci
    oxresult = oxcopy.copy()
    reg_reg_edges = regcopy.copy()
    reg_something_edges = sthcopy.copy()
    other_edges = otcopy.copy()
    in_edges = incopy.copy()
    dirname = 'spis_powszechny/'+"k"+str(k)
    
    if not os.path.isdir('spis_powszechny'):
        os.makedirs('spis_powszechny')
    
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    
    presult = Pr.create_vector()
    
    for i in range(iters):
        print(f'Iter {i + 1}/{iters}')
    
        pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
        oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges)
    
        pnow = Pr.solve_equation(pmatrix, presult)
        
        oxnow = Ox.solve_equation(oxmatrix, oxresult)
    
        vresult = Ve.create_vector(oxnow, oxresult)
        vmatrix = Ve.update_matrix(vresult, reg_reg_edges, reg_something_edges, other_edges)
        vnow = Ve.solve_equation(vmatrix, vresult)
        for node in other_nodes:
            vnow[node] = 0
    
        if i%save_every == 0:
            G = Pr.update_network(G, reg_reg_edges, reg_something_edges, pnow)
            
    #        Dr.drawhist(name = f'{i//save_every:04d}.png',  oxnow = oxnow, oxresult = oxresult, vnow = vnow, dirname = dirname)
#            Dr.drawd(name = f'd{(i+old_iters)//save_every:04d}.png',  oxdraw = [], dirname = dirname)
#            Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png',  oxdraw = [], dirname = dirname)
#            Dr.drawq(name=f'veq{i+old_iters // save_every:04d}.png',  oxdraw=vnow/np.max(vnow), dirname = dirname)
#            Dr.drawq(name=f'oxq{i+old_iters // save_every:04d}.png',  oxdraw=oxresult, dirname = dirname)
            Dr.drawq(name=f'veq{(i+old_iters) // save_every:04d}.png', oxdraw=vnow / np.max(vnow)+oxresult, dirname = dirname)
            Dr.drawblood(name=f'q_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='q', dirname = dirname)
 #           Dr.drawblood(name=f'd_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='d', dirname = dirname)
            
    
    
        reg_reg_edges, reg_something_edges, in_edges = Pr.update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)
        reg_reg_edges, reg_something_edges, in_edges, oxresult=Ve.update_graph(vnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)
    #    oxresult = Ox.update_oxresult(reg_reg_edges, reg_something_edges, in_edges, oxresult)      #update oxresult gdy vegf jest wylaczony
    
    
     
    
    Sv.save_all(dirname+'/save', reg_reg_edges, reg_something_edges, other_edges, oxresult,
                n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
                F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
                in_nodes_ox, out_nodes_ox,in_edges,G, boundary_nodes_out, boundary_nodes_in, boundary_edges)
    
    k = 10 * k