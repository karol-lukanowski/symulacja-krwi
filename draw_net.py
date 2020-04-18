#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# normalizacja rysowania (maksymalna grubość krawędzi)
qdrawconst = 5
ddrawconst = 3


def drawq(G, n, name, qdrawconst=qdrawconst, normalize=True, in_nodes=[], out_nodes=[]):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')

    if (normalize == False):
        qmax = 1
        qdrawconst = 1
    else:
        qmax = max([edge[2] for edge in G.edges(data='q')])

    edges = []
    qs = []
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x, y, q = edge
            edges.append((x, y))
            qs.append(q)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=qdrawconst * np.array(qs) / qmax)

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
    plt.savefig(name)
    plt.close()


def drawd(G, n, name, ddrawconst=ddrawconst, normalize=True, in_nodes=[], out_nodes=[]):
    """
    rysowanie srednic
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')

    if (normalize == False):
        dmax = 1
        ddrawconst = 1
    else:
        dmax = max([edge[2] for edge in G.edges(data='d')])

    edges = []
    ds = []
    for edge in G.edges(data='d'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x, y, d = edge
            edges.append((x, y))
            ds.append(d)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=ddrawconst * np.array(ds) / dmax)

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
    plt.savefig(name)
    plt.close()
