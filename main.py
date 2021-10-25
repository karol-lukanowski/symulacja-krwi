import numpy as np
import draw_net as Dr
import pressure as Pr
import oxygen as Ox
import vegf as Ve
import save as Sv
import pruning as Prun
import upstream as Up
import analysis as An
import blood_oxygen as Bo
from build import (G, in_nodes, out_nodes, reg_nodes, other_nodes, iters, old_iters,
                   in_edges, boundary_nodes_out, boundary_nodes_in, dirname,
                   n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
                   F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes_ox, out_nodes_ox, oxresult,
                   boundary_edges, reg_reg_edges, reg_something_edges, other_edges)
from config import save_every


presult = Pr.create_vector()


#prestab = []
#oxtab = []
#vegftab = []

#bloodoxresult = Bo.create_vector()
#oxresult = Ve.update_blood(oxresult, reg_reg_edges, reg_something_edges, in_edges)

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    #oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges)
    

    pnow = Pr.solve_equation(pmatrix, presult)
    

    #bloodoxmatrix = Bo.update_matrix(pnow, oxresult, bloodoxresult, reg_reg_edges, reg_something_edges, other_edges)
    #bloodoxnow = Bo.solve_equation(bloodoxmatrix, bloodoxresult)
    #oxnow = Ox.solve_equation(oxmatrix, bloodoxnow)
    #oxnow = Ox.solve_equation(oxmatrix, oxresult)

    #vresult = Ve.create_vector(G, n, oxnow, oxresult)
    vresult = Ve.create_vector(G, n, oxresult)
    vmatrix = Ve.update_matrix(vresult, reg_reg_edges, reg_something_edges, other_edges)
    vnow = Ve.solve_equation(vmatrix, vresult)

    for node in other_nodes:
        vnow[node] = 0
    
    smatrix, sresult = Up.update_matrix_upstream(vnow, pnow, oxresult,reg_reg_edges, reg_something_edges, other_edges)
    snow_upstream = Up.solve_equation(smatrix, sresult)
    smatrix, sresult = Up.update_matrix_downstream(vnow, pnow, oxresult, reg_reg_edges, reg_something_edges, other_edges)
    snow_downstream = Up.solve_equation(smatrix, sresult)

    snow = snow_upstream+snow_downstream


    
    if i%save_every == 0:
        vnow2 = vnow.copy()
        for node in out_nodes:
            vnow2[node] = 0

        #for node in in_nodes:
        #    snow[node] = 0

        G = Pr.update_network(G, reg_reg_edges, reg_something_edges, pnow)

        gradp = np.zeros(n**2)
        pin = np.max(pnow)
        for n1, n2, d, l in reg_reg_edges+reg_something_edges+in_edges:
            if np.abs((pnow[n1] - pnow[n2]) / (pin * l)) > gradp[n1]:
                gradp[n1] = np.abs((pnow[n1] - pnow[n2]) / (pin * l))
            if np.abs((pnow[n1] - pnow[n2]) / (pin * l)) > gradp[n2]:
                gradp[n2] = np.abs((pnow[n1] - pnow[n2]) / (pin * l))
            
        Dr.drawhist(f'{(i+old_iters) // save_every:04d}.png', oxresult = oxresult, vnow = gradp, oxdraw = oxresult)
#        Dr.drawd(name = f'd{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=snow)



        Dr.drawblood(name=f'q_bloodp2{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, oxdraw = gradp, data='q')
        Dr.drawblood(name=f'd_bloodp2{(i+old_iters) // save_every:04d}.png', oxresult=oxresult, oxdraw = gradp, data='d')

            
    if i%100 == 0:
        Sv.save_all(dirname+'/save'+f'{i+old_iters}', reg_reg_edges, reg_something_edges, other_edges, oxresult,
             n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
             F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
             in_nodes_ox, out_nodes_ox,in_edges,G, boundary_nodes_out, boundary_nodes_in, boundary_edges)    

    reg_reg_edges, reg_something_edges, in_edges = Pr.update_graph(pnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)    
    reg_reg_edges, reg_something_edges, in_edges, oxresult = Ve.update_graph(vnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)
    reg_reg_edges, reg_something_edges, in_edges = Up.update_graph_upstream(snow_upstream, pnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)
    reg_reg_edges, reg_something_edges, in_edges = Up.update_graph_downstream(snow_downstream, pnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)
    reg_reg_edges, reg_something_edges, in_edges = Pr.update_graph_gradp(pnow, reg_reg_edges, reg_something_edges, in_edges)

    
    #prestab.append(np.max(pnow))
    #oxtab.append(np.average(oxnow ** 2))
    #vegftab.append(np.average(vnow ** 2))

#Dr.plot_params(prestab, oxtab, vegftab)


# Sv.save_all(dirname+'/save', reg_reg_edges, reg_something_edges, other_edges, oxresult,
#             n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth, iters+old_iters,
#             F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, in_nodes, out_nodes, reg_nodes, other_nodes,
#             in_nodes_ox, out_nodes_ox,in_edges,G, boundary_nodes_out, boundary_nodes_in, boundary_edges)

#An.getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname)


#DiG = An.getDiGraph(G, pnow, oxresult, in_nodes)
#DiG = An.strahlerOrder(DiG)
#An.plotStrahlerGraph(DiG, 'deown101deown_drabina/26')

#Prun.pruning(G, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult)



