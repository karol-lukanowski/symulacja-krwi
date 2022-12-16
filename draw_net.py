import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec

from config import simInputData


def uniform_hist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, diams, flow, blood_vessels, oxygen, name, draw = 'd'):
    if draw == 'd':
        qs = diams * (1 - boundary_edges)
    else:
        qs = flow * (1 - boundary_edges)
        
    shear = np.abs(flow) / diams ** 3

    if sid.normalize:
        qs = 5 * qs / np.max(qs)

    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])

    
    # x_tr, y_tr = [], []
    # for node in triangles_pos:
    #     x_tr.append(node[0])
    #     y_tr.append(node[1])

    cols = 5
    plt.figure(figsize=(sid.figsize * 1.25, sid.figsize))
    spec = gridspec.GridSpec(ncols=cols, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=cols))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    #plt.scatter(x_tr, y_tr, s=1*(vols < 9), facecolors='red', edgecolors='black')
    nx.draw_networkx_edges(G, pos, edge_color = blood_vessels, width=sid.ddrawconst * np.array(qs))
    #nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.ddrawconst * np.array(qs), edge_color=colors)
    nx.draw_networkx_nodes(G, pos, node_size = 1, node_color = oxygen, cmap='Blues')
    plt.axis('equal')
    
    plt.subplot(spec[cols]).set_title('Diameter')
    plt.hist(diams, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 1]).set_title('Flow')
    plt.hist(np.abs(flow), bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 2]).set_title('Shear')
    plt.hist(shear, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 3]).set_title('Oxygen')
    plt.hist(oxygen, bins=50)
    plt.yscale("log")

    # plt.subplot(spec[cols + 4]).set_title('vola')
    # plt.hist(vols, bins=50)
    # plt.yscale("log")
    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()
