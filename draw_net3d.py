import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from config import simInputData

def drawvessels(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, name, oxresult, oxdraw, data='q'):
    """
    rysowanie krwi, data to q albo d
    """
    plt.figure(figsize=(40, 40), dpi = 200)
    plt.axis('off')
    pos = nx.get_node_attributes(G, 'pos')

    qmax = max([edge[2] for edge in G.edges(data=data)])

    edges = []
    qs = []
    ds = []
    for edge in G.edges(data=data):
            x, y, q = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
                edges.append((x, y))
                qs.append(q)
                ds.append(q)

    #draw only those between oxygen nodes:
    # for i, edge in enumerate(edges):
    #     n1, n2 = edge
    #     if (oxresult[n1] != 1 or oxresult[n2] != 1):
    #         qs[i] = 0
    #         ds[i] = 0


    
    vessels = np.zeros(2 * sid.nsq)

    def find_veins(G, node, oxresult):
        vessels[node] = 2
        for neigh in G.neighbors(node):
            if oxresult[neigh] == 1 and vessels[neigh] == 0 and neigh != node - sid.nsq:
                find_veins(G, neigh, oxresult)

    def find_arts(G, node, oxresult):
        vessels[node] = 1
        for neigh in G.neighbors(node):
            if oxresult[neigh] == 1 and vessels[neigh] == 0 and neigh != node + sid.nsq:
                find_arts(G, neigh, oxresult)

    for node in in_nodes:
        find_arts(G, node, oxresult)
    for node in out_nodes:
        find_veins(G, node, oxresult)

    colors = []
    for edge in edges:
        n1, n2 = edge
        if vessels[n1] == 1 and vessels[n2] == 1:
            colors.append('r')
        elif vessels[n1] == 2 and vessels[n2] == 2:
            colors.append('b')
        else:
            colors.append('k')

    # colors = []
    # for edge in edges:
    #     n1, n2 = edge
    #     if oxresult[n1] == 1 and oxresult[n2] == 1:
    #         if n1 < sid.nsq and n2 < sid.nsq:
    #             colors.append('r')
    #         elif n1 > sid.nsq and n2 > sid.nsq:
    #             colors.append('b')
    #     else:
    #         colors.append('k')

    # for i, edge in enumerate(G.edges(data=data)):
    #     if edge[2] < G[edge[0]][edge[1]]['q']:
    #         qs[i] = 0
        

    if data == 'q':
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.qdrawconst * np.array(qs) / qmax, edge_color=colors)
    else:
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.ddrawconst * np.array(ds) / qmax, edge_color=colors)
    #nx.draw_networkx_nodes(G, pos, node_size = 15, node_color = oxdraw, cmap='coolwarm')

    #### IN_NODES i OUT_NODES ####
    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')

    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    
    plt.axis('equal')
    plt.xticks([], [])
    plt.yticks([], [])
    plt.savefig(sid.dirname + "/" + name)
    plt.close()