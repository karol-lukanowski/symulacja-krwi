#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import networkx as nx

# normalizacja rysowania (maksymalna grubość krawędzi)
qdrawconst = 5
ddrawconst = 3

def drawq(G, n, name, qdrawconst = qdrawconst, normalize = True):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    #nx.draw(G, pos, node_size=0, with_labels=False)
    
    qmax = 0
    for edge in G.edges(data='q'):
        if (edge[2] > qmax):
            qmax=edge[2]

    if (normalize == False):
        qmax = 1
        qdrawconst = 1
    Ed = []
    Qu = []
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x,y,q = edge
            Ed.append((x,y))
            Qu.append(q)
    nx.draw_networkx_edges(G, pos, edgelist=Ed, width=qdrawconst * Qu / qmax)
    
    plt.axis('equal')
    plt.savefig(name)
    plt.close()

def drawd(G, n, name, ddrawconst = ddrawconst, normalize = True):
    """
    rysowanie srednic
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    #nx.draw(G, pos, node_size=0)
    
    dmax = 0
    for edge in G.edges(data='d'):
        if (edge[2] > dmax):
            dmax=edge[2]

    if (normalize == False):
        dmax = 1
        ddrawconst = 1
    Ed = []
    Dia = []
    for edge in G.edges(data='d'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x,y,d = edge
            Ed.append((x,y))
            Dia.append(d)
    nx.draw_networkx_edges(G, pos, edgelist=Ed, width=qdrawconst * Dia / dmax)

    plt.axis('equal')
    plt.savefig(name)
    plt.close()
