#import analysis as An
import blood_oxygen as Bo
import decrease as Dc
import draw_net as Dr
import new_oxygen as Ox
#import bloodox_everywhere as Ox
#import oxygen as Ox
import pressure as Pr
#import pruning as Prun
import save as Sv
import upstream as Up
import vegf as Ve

from build import build
from utils import solve_equation, simAnalysisData, collect_data

from config import simInputData

from delaunay import find_node


sid = simInputData()
sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = build(sid)
sad = simAnalysisData()
out_node = find_node(G, sid.out_nodes_own[0])

oxnowtab = []

presult = Pr.create_vector(sid, in_nodes, out_nodes)
if sid.oxygen:
    #bloodoxresult_a = Bo.create_vector_a(sid, in_nodes_ox)
    #bloodoxresult_v = Bo.create_vector_a(sid, out_nodes_ox)
    bloodoxresult = Ox.create_vector2(sid, in_nodes_ox)

iters = sid.old_iters + sid.iters
for i in range(sid.old_iters, iters):
    print(f'Iter {i + 1}/{iters}')

    pmatrix = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    pnow = solve_equation(pmatrix, presult)
    
    
    if sid.oxygen:
        # bloodoxmatrix = Bo.update_matrix_a(sid, pnow, oxresult, bloodoxresult_a, edges)
        # bloodoxnow_a = solve_equation(bloodoxmatrix, bloodoxresult_a)
        # bloodoxmatrix = Bo.update_matrix_v(sid, pnow, oxresult, bloodoxresult_v, edges)
        # bloodoxnow_v = solve_equation(bloodoxmatrix, bloodoxresult_v)
        # bloodoxnow = bloodoxnow_a + bloodoxnow_v
        # oxmatrix = Ox.update_matrix(sid, oxresult, edges)
        # oxnow = solve_equation(oxmatrix, bloodoxnow_a)
        oxmatrix = Ox.update_matrix(sid, oxresult, bloodoxresult, pnow, edges)
        oxnow = solve_equation(oxmatrix, bloodoxresult)
        vresult = Ve.create_vector(sid, oxresult, oxnow)
    else:
        vresult = Ve.create_vector(sid, oxresult)

    #vresult = Ve.create_vector(sid, oxnow, oxresult)
    vmatrix = Ve.update_matrix(sid, vresult, edges)
    vnow = solve_equation(vmatrix, vresult)

    # for node in other_nodes:
    #     vnow[node] = 0

    if sid.signal:
        smatrix, sresult = Up.update_matrix_upstream(sid, vnow, pnow, edges, in_nodes)
        snow_upstream = solve_equation(smatrix, sresult)
        smatrix, sresult = Up.update_matrix_downstream(sid, vnow, pnow, edges, out_nodes)
        snow_downstream = solve_equation(smatrix, sresult)
        snow = snow_upstream+snow_downstream



    if i % sid.plot_every == 0:
        # print (oxnow[out_node])
        # oxnowtab.append(oxnow[out_node])
        # vnow2 = vnow.copy()
        # for node in out_nodes:
        #     vnow2[node] = 0

        G = Pr.update_network(G, sid, edges, pnow)

        #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, oxnow, vnow, snow_upstream, snow, name=f'histogram{sid.old_iters // sid.plot_every:04d}.png')
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, pnow, pnow, pnow, pnow, name=f'histogram{sid.old_iters // sid.plot_every:04d}.png')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, name=f'd{sid.old_iters // sid.save_every:04d}.png', data='d')
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=snow)
#        Dr.drawblood(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_blood{sid.old_iters // sid.save_every:04d}.png', oxresult=oxresult, oxdraw = vnow, data='q')
        #Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=oxresult, oxdraw = oxnow, data='q')
#        Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=oxresult, oxdraw = oxnow, data='d')
    
    if sid.shear_d:
        edges = Pr.update_graph(sid, edges, pnow)    
    if sid.vegf_d:
        edges = Ve.update_graph(sid, vnow, oxresult, edges)
    if sid.signal_d:
        edges = Up.update_graph_upstream(sid, snow_upstream, pnow, oxresult, edges)
        edges = Up.update_graph_downstream(sid, snow_downstream, pnow, oxresult, edges)
    if sid.gradp_d:
        edges = Pr.update_graph_gradp(sid, edges, pnow)
        #edges = Pr.update_graph_gradp_tips(sid, edges, pnow, oxresult)
    if sid.decrease_d:
        edges = Dc.update_graph(sid, edges)
    oxresult = Ve.update_blood(sid, oxresult, edges)

    if sid.data_collection:
        collect_data(sad, sid, in_nodes, pnow, vnow, oxnow)


    if i % sid.save_every == 0 and i != 0:
        Sv.save(f'/save{sid.old_iters}.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)

    sid.old_iters += 1

Sv.save('/save.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
if sid.data_collection:
    Dr.plot_params(sid)

import numpy as np
np.savetxt(sid.dirname+'/oxout.txt', oxnowtab)

#An.getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname)


#DiG = An.getDiGraph(G, pnow, oxresult, in_nodes)
#DiG = An.strahlerOrder(DiG)
#An.plotStrahlerGraph(DiG, 'deown101deown_drabina/26')

#Prun.pruning(G, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult)