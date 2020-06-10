import numpy as np
import os
import draw_net as Dr
import pressure as Pr
import oxygen as Ox
import vegf as Ve
import save as Sv
from config import nettype, geo
from build import (G, in_nodes, out_nodes, reg_nodes, other_nodes, iters, old_iters,
                   in_edges, save_every, boundary_nodes_out, boundary_nodes_in, dirname, save_name,
                   n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                   F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes_ox, out_nodes_ox, oxresult,
                   boundary_edges, reg_reg_edges, reg_something_edges, other_edges)

Initial_everything = (G, in_nodes, out_nodes, reg_nodes, other_nodes,
                   in_edges, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult,
                   boundary_edges, reg_reg_edges, reg_something_edges, other_edges)

X = 12
'''
dt = 0.9
dt_ox = 0.45


D = 1 # współczynnik dyfuzji
Dv = 0.5 # współczynnik dyfuzji VEGF

k = 0.1 # stała reakcji
dth = 3 # graniczna grubosć

'''


for i in range(X):
    # wyzerowanie sieci
    Initial_everything = (G, in_nodes, out_nodes, reg_nodes, other_nodes,
                       in_edges, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult,
                       boundary_edges, reg_reg_edges, reg_something_edges, other_edges)
# Zmienne dt_ox
    '''
    if i < 3:
        dt_ox *= 2 
        dirname = 'spis_powszechny/'+nettype + geo + "n" + str(n) + "dt" + str(dt) + "dtox" + str(dt_ox)
# Zmienne Dv, dt_ox już ustalone
    elif i < 6:
        dt_ox = 3
        Dv *= 2
        dirname = 'spis_powszechny/'+nettype + geo + "n" + str(n) + "D" + str(D)+ "Dv" + str(Dv)

# Zmienne k, Dv już ustalone

    elif i < 9:
        Dv = 1
        k += 1
        dirname = 'spis_powszechny/'+nettype + geo + "n" + str(n) + "k" + str(k)

# Zmienne dth, k już ustalone

    elif i < 12:
        k = 0.1
        dth += 5
        dirname = 'spis_powszechny/'+nettype + geo + "n" + str(n) + "dth" + str(dth)
    
    if not os.path.isdir('spis_powszechny'):
        os.makedirs('spis_powszechny')
    
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    '''
    presult = Pr.create_vector()
    
    for i in range(iters):
        print(f'Iter {i + 1}/{iters}')
    
        pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
        oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges, dt,dt_ox,D,Dv,k,dth)
    
        pnow = Pr.solve_equation(pmatrix, presult)
        print(pnow)
        oxnow = Ox.solve_equation(oxmatrix, oxresult)
    
        vresult = Ve.create_vector(oxnow, oxresult)
        vmatrix = Ve.update_matrix(vresult, reg_reg_edges, reg_something_edges, other_edges, dt,dt_ox,D,Dv,k,dth)
        vnow = Ve.solve_equation(vmatrix, vresult)
        for node in other_nodes:
            vnow[node] = 0
    
        if i%save_every == 0:
            G = Pr.update_network(G, reg_reg_edges, reg_something_edges, pnow)
            
    #        Dr.drawhist(dirname,name = f'{i//save_every:04d}.png',  oxnow = oxnow, oxresult = oxresult, vnow = vnow)
            Dr.drawd(dirname,name = f'd{(i+old_iters)//save_every:04d}.png',  oxdraw = [])
            Dr.drawq(dirname,name = f'q{(i+old_iters)//save_every:04d}.png',  oxdraw = [])
            Dr.drawq(dirname,name=f'veq{i // save_every:04d}.png',  oxdraw=vnow/np.max(vnow))
            Dr.drawq(dirname,name=f'oxq{i // save_every:04d}.png',  oxdraw=oxresult)
            Dr.drawq(dirname,name=f'veq{(i+old_iters) // save_every:04d}.png', oxdraw=vnow / np.max(vnow)+oxresult)
            Dr.drawblood(dirname,name=f'q_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='q')
            Dr.drawblood(dirname,name=f'd_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='d')
            
    
    
        reg_reg_edges, reg_something_edges, in_edges = Pr.update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges, dt)
        reg_reg_edges, reg_something_edges, in_edges, oxresult=Ve.update_graph(vnow, oxresult, reg_reg_edges, reg_something_edges, in_edges, dt,dt_ox,D,Dv,k,dth)
    #    oxresult = Ox.update_oxresult(reg_reg_edges, reg_something_edges, in_edges, oxresult)      #update oxresult gdy vegf jest wylaczony
    
    
     
    
    Sv.save_all(dirname+'/'+save_name, reg_reg_edges, reg_something_edges, other_edges, oxresult,
                n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
                F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
                in_nodes_ox, out_nodes_ox,in_edges,G, boundary_nodes_out, boundary_nodes_in, boundary_edges)
