import data as Da
import draw_net as Dr
import growth as Gr
import incidence as In
import pressure as Pr
import save as Sv

from build import build
from utils import initialize_iterators, update_iterators

from config import simInputData
import numpy as np

sid = simInputData()
sid, G, in_nodes, out_nodes, boundary_edges = build(sid)
pressure_b = Pr.create_vector(sid, in_nodes)
inc_matrix, mid_matrix, bound_matrix, in_matrix, diams, lens, in_edges, \
    out_edges, edge_list, boundary_edge_list = In.create_matrices(sid, G, in_nodes, out_nodes, boundary_edges)

iters, tmax, i, t, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')

    pressure, flow = Pr.find_flow(sid, diams, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_nodes)
    if i % sid.plot_every == 0:
        Da.check_flow(flow, in_edges, out_edges)
        G = Pr.update_network(G, edge_list, diams, flow)
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edge_list, diams, flow, name=f'd_{sid.old_iters:.2f}_{sid.old_t:.2f}.png', draw = 'd')
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edge_list, diams, flow, name=f'q_{sid.old_iters:.2f}_{sid.old_t:.2f}.png', draw = 'q')
        
    # diams, dt = update_diameters(sid, flow, cb, diams, lens, inc_matrix, out_edges)
    diams, dt = Gr.update_diameters(sid, flow, diams)

    i, t = update_iterators(sid, i, t, dt)
    Da.collect_data(sid, pressure)

if i != 1:
    G = Pr.update_network(G, edge_list, diams, flow)
    Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edge_list, diams, flow, name=f'd_{sid.old_iters:.2f}_{sid.old_t:.2f}.png', draw = 'd')
    Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edge_list, diams, flow, name=f'q_{sid.old_iters:.2f}_{sid.old_t:.2f}.png', draw = 'q')
    Da.check_flow(flow, in_edges, out_edges)
    Sv.save('/save.dill', sid, G, in_nodes, out_nodes, boundary_edges)
    Da.plot_data(sid)