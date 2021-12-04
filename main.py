import analysis2 as An
import blood_oxygen as Bo
import decrease as Dc
import draw_net as Dr
import oxygen as Ox
import pressure2 as Pr
#import pruning as Prun
import save as Sv
import upstream as Up
import vegf as Ve

from build import build
from utils import solve_equation

from config import simInputData



sid = simInputData()
sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = build(sid)
sad = An.simAnalysisData()

presult = Pr.create_vector(sid, in_nodes, out_nodes)
if sid.oxygen:
    bloodoxresult = Bo.create_vector(sid, in_nodes_ox)

iters = sid.old_iters + sid.iters
for i in range(sid.old_iters, iters):
    print(f'Iter {i + 1}/{iters}')

    pmatrix = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    pnow = solve_equation(pmatrix, presult)
    
    
    if sid.oxygen:
        bloodoxmatrix = Bo.update_matrix(sid, pnow, oxresult, bloodoxresult, edges)
        bloodoxnow = solve_equation(bloodoxmatrix, bloodoxresult)
        oxmatrix = Ox.update_matrix(sid, oxresult, edges)
        oxnow = solve_equation(oxmatrix, bloodoxnow)
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



    if i%sid.save_every == 0:
        vnow2 = vnow.copy()
        for node in out_nodes:
            vnow2[node] = 0

        G = Pr.update_network(G, sid, edges, pnow)

        # gradp = np.zeros(n**2)
        # pin = np.max(pnow)
        # for n1, n2, d, l in reg_reg_edges+reg_something_edges+in_edges:
        #     if np.abs((pnow[n1] - pnow[n2]) / (pin * l)) > gradp[n1]:
        #         gradp[n1] = np.abs((pnow[n1] - pnow[n2]) / (pin * l))
        #     if np.abs((pnow[n1] - pnow[n2]) / (pin * l)) > gradp[n2]:
        #         gradp[n2] = np.abs((pnow[n1] - pnow[n2]) / (pin * l))
            
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, vnow, snow_upstream, snow, name=f'histogram{sid.old_iters // sid.save_every:04d}.png')
        Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, name=f'd{sid.old_iters // sid.save_every:04d}.png', data='d')
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=snow)
#        Dr.drawblood(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_blood{sid.old_iters // sid.save_every:04d}.png', oxresult=oxresult, oxdraw = vnow, data='q')
#        Dr.drawblood(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_blood{sid.old_iters // sid.save_every:04d}.png', oxresult=oxresult, oxdraw = snow, data='d')
    
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
        An.collect_data(sad, sid, in_nodes, pnow, vnow, snow)


    sid.old_iters += 1

Sv.save('/save.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
#An.plot_data(sid)

#An.getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname)


#DiG = An.getDiGraph(G, pnow, oxresult, in_nodes)
#DiG = An.strahlerOrder(DiG)
#An.plotStrahlerGraph(DiG, 'deown101deown_drabina/26')

#Prun.pruning(G, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult)