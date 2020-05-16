#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np



from matplotlib import gridspec
from config import G, n, qdrawconst, ddrawconst, in_nodes, out_nodes, F_mult_ox, F_mult, c2
import pressure as Pr
import oxygen as Ox

# normalizacja rysowania (maksymalna grubość krawędzi)



def drawq(name, normalize=True, oxresult = []):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(20, 20))
    pos = nx.get_node_attributes(G, 'pos')

    if (normalize == False):
        qmax = 1
        drawconst = 1
    else:
        qmax = max([edge[2] for edge in G.edges(data='q')])
        drawconst = qdrawconst

    edges = []
    qs = []
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x, y, q = edge
            edges.append((x, y))
            qs.append(q)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=drawconst * np.array(qs) / qmax)

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

    #### OXYGEN NODES ####
    oxresult = np.abs(np.array(oxresult))
    oxresult=oxresult.astype(int)
    oxresult2 = oxresult-0.5
    oxresult = list(oxresult)

#    nx.draw_networkx_nodes(G, pos, node_size=25 * oxresult2, node_color=oxresult, cmap = 'bwr')

    plt.axis('equal')
    plt.savefig(name)
    plt.close()    
    



def drawd(name, normalize=True, oxresult = []):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(20, 20))
    pos = nx.get_node_attributes(G, 'pos')

    if (normalize == False):
        dmax = 1
        drawconst = 1
    else:
        dmax = max([edge[2] for edge in G.edges(data='d')])
        drawconst = ddrawconst

    edges = []
    ds = []
    for edge in G.edges(data='d'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x, y, d = edge
            edges.append((x, y))
            ds.append(d)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=drawconst * np.array(ds) / dmax)

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

    #### OXYGEN NODES ####
    oxresult = np.abs(np.array(oxresult))
    oxresult=oxresult.astype(int)
    oxresult2 = oxresult-0.5
    oxresult = list(oxresult)

#    nx.draw_networkx_nodes(G, pos, node_size=25 * oxresult2, node_color=oxresult, cmap = 'bwr')

    plt.axis('equal')
    plt.savefig(name)
    plt.close()    
    




def drawhist(name, oxnow=[], oxresult=[]):
    """
    rysowanie przepływów
    """
    dhist=[[], [], [], [], [], [], []]
    qhist=[[], [], [], [], [], [], []]
    shearhist=[[], [], [], [], [], [], []]
    Fhist=[[], [], [], [], [], [], []]
    shearhistox=[[], [], [], [], [], [], []]

    colors=[]
    dmax=1
    shearmax=0.001
    Fmax=0
    shearmaxox=0.001

    
    qmax = max([edge[2] for edge in G.edges(data='q')])
    for n1, n2 in G.edges():
        q=G[n1][n2]['q']
        d=G[n1][n2]['d']
        F=F_mult*c2*q/d**3
        shear=Pr.d_update(F)
        if d>dmax:
            dmax=d
        if shear>shearmax:
            shearmax=shear
        if F>Fmax:
            Fmax=F
        F_ox=F_mult_ox*np.abs(oxnow[n1] - oxnow[n2])
        shearox=Ox.d_update(F_ox)
        if shearox>shearmaxox:
            shearmaxox=shearox
        qhist[int(6*q/qmax)].append(q)
        dhist[int(6*q/qmax)].append(d)
        Fhist[int(6*q/qmax)].append(F)
        shearhist[int(6*q/qmax)].append(shear)
        shearhistox[int(6*q/qmax)].append(shearox)
    colortab = []

    
    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',]
    for edge in G.edges(data="q"):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            colors.append(color[int(6*edge[2]/qmax)])
        
    pos = nx.get_node_attributes(G, 'pos')




    edges = []
    qs = []
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            x, y, q = edge
            edges.append((x, y))
            qs.append(q)

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    
    oxresult2=oxresult.copy()        
    oxresult2 = np.abs(np.array(oxresult2))
    oxresult2 = oxresult2.astype(int)
    oxresult3 = oxresult2-0.5
    oxresult2 = list(oxresult2)
    
    
    plt.figure(figsize=(20, 20))
    plt.suptitle('Flow for n = '+str(n))
    spec = gridspec.GridSpec(ncols=4, nrows=2, height_ratios=[5, 1])
    plt.subplot(spec.new_subplotspec((0, 0), colspan=4))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=qdrawconst * np.array(qs) / qmax)
    nx.draw_networkx_nodes(G, pos, node_size=25 * oxresult3, node_color=colortab, cmap = 'bwr')
    plt.axis('equal')
    plt.subplot(spec[4]).set_title('Diameter')
    cindex=0
    plt.xlim((1,1.1*dmax))
    for hist in dhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    plt.subplot(spec[5]).set_title('Flow')
    cindex=0
    plt.xlim((0, 1.1*qmax))
    for hist in qhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    plt.subplot(spec[6]).set_title('Oxygen')
    cindex=0
#    plt.xlim((0, 1.1*shearmaxox))
#    for hist in shearhistox:
#        if len(hist)>1:
#            plt.hist(hist, bins=50, color=color[cindex])
#        cindex+=1
    plt.hist(oxnow, bins=50)
    plt.yscale("log")        
    plt.subplot(spec[7]).set_title('Shear')
    cindex=0
    plt.xlim((0,1.1*shearmax))
    for hist in shearhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log") 
    plt.savefig(name)
    plt.close()
        
    