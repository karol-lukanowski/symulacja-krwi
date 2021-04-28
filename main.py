import numpy as np
import draw_net as Dr
import pressure as Pr
import oxygen as Ox
import vegf as Ve
import save as Sv
import pruning as Prun
import analysis as An
from build import (G, in_nodes, out_nodes, reg_nodes, other_nodes, iters, old_iters,
                   in_edges, boundary_nodes_out, boundary_nodes_in, dirname,
                   n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                   F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes_ox, out_nodes_ox, oxresult,
                   boundary_edges, reg_reg_edges, reg_something_edges, other_edges)
from config import save_every


presult = Pr.create_vector()


prestab = []
oxtab = []
vegftab = []

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
        vnow2 = vnow.copy()
        for node in out_nodes:
            vnow2[node] = 0
        G = Pr.update_network(G, reg_reg_edges, reg_something_edges, pnow)
        
        Dr.drawhist(name = f'{(i+old_iters) // save_every:04d}.png', oxnow = oxnow, oxresult = oxresult, vnow = vnow, oxdraw = oxresult)
#        Dr.drawd(name = f'd{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=oxnow)
#        Dr.drawq(name=f'veq{(i+old_iters) // save_every:04d}.png', oxdraw=vnow2 / np.max(vnow2)+oxresult)
        Dr.drawblood(name=f'q_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='q')
        Dr.drawblood(name=f'd_blood{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, data='d')
            


    reg_reg_edges, reg_something_edges, in_edges = Pr.update_graph(G, pnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)    
    reg_reg_edges, reg_something_edges, in_edges, oxresult = Ve.update_graph(vnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)
    prestab.append(np.max(pnow))
    oxtab.append(np.average(oxnow ** 2))
    vegftab.append(np.average(vnow ** 2))

Dr.plot_params(prestab, oxtab, vegftab)


Sv.save_all(dirname+'/save', reg_reg_edges, reg_something_edges, other_edges, oxresult,
            n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
            F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
            in_nodes_ox, out_nodes_ox,in_edges,G, boundary_nodes_out, boundary_nodes_in, boundary_edges)

An.getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname)


Prun.pruning(G, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult)



