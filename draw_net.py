#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec
from config import qdrawconst, ddrawconst, boundary_edges, dirname, Load, load_name, n, in_nodes, out_nodes, F_mult_ox, F_mult, c2, dt, dt_ox, dth
import pressure as Pr
import oxygen as Ox
import save as Sv
import vegf as Ve
'''
if Load == True:
    (G, n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
     F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, iters,
     in_nodes, out_nodes, reg_nodes, other_nodes, in_nodes_ox, out_nodes_ox, oxresult,
     in_edges, reg_reg_edges, reg_something_edges, other_edges)  = Sv.load(dirname+'/'+load_name)
    nkw = n**2
else:
    from config import G, n, in_nodes, out_nodes, F_mult_ox, F_mult, c2, dt, dt_ox, dth
'''
# normalizacja rysowania (maksymalna grubość krawędzi)



def drawq(G,name, normalize=True, oxdraw=[]):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(20, 20))
    pos = nx.get_node_attributes(G, 'pos')
    new_pos = {}
    for x in pos:
        new_pos[x] = tuple(pos[x])
    pos = new_pos
    #print(G.edges())
#    print(nx.get_edge_attributes(G, 'q'))
    if (normalize == False):
        qmax = 1
        drawconst = 1
    else:
        qmax = max([edge[2] for edge in G.edges(data='q')])
        drawconst = qdrawconst

    edges = []
    qs = []
    for edge in G.edges(data='q'):
        x, y, q = edge
        
        if (x, y) not in boundary_edges and (y, x) not in boundary_edges:            
            edges.append((x, y))
            qs.append(q)
    print(pos)
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
    nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')

    #nx.draw(G, pos, with_labels=True)
    #plt.show()

    plt.axis('equal')
    plt.savefig(dirname + "/" + name)
    plt.close()




def drawd(G,name, normalize=True, oxdraw = []):
    """
    rysowanie srednic
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
            x, y, d = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
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
    nx.draw_networkx_nodes(G, pos, node_size=25 * oxdraw, node_color=oxdraw, cmap='Reds')


    plt.axis('equal')
    plt.savefig(dirname + "/" + name)
    plt.close()    
    




def drawhist(G,name, oxnow=[], oxresult=[], vnow = [], oxdraw = []):
    """
    rysowanie histogramów
    """
    dhist=[[], [], [], [], [], [], []]
    qhist=[[], [], [], [], [], [], []]
    shearhist=[[], [], [], [], [], [], []]
    shearhistox=[[], [], [], [], [], [], []]
    colors=[]
    dmax=1
    shearmax=0.001
    shearmaxox=0.001
    qmax = max([edge[2] for edge in G.edges(data='q')])
    for n1, n2 in G.edges():
        q=G[n1][n2]['q']
        d=G[n1][n2]['d']
        F=F_mult*c2*q/d**3 * dt
        shear=Pr.d_update(F)
        F_ox=F_mult_ox*np.abs(vnow[n1] - vnow[n2]) * dt_ox
        shearox=Ve.d_update(F_ox)
        if d>dmax:
            dmax=d
        if shear>shearmax:
            shearmax=shear
        if shearox>shearmaxox:
            shearmaxox=shearox
        qhist[int(6*q/qmax)].append(q)
        dhist[int(6*q/qmax)].append(d)
        shearhist[int(6*q/qmax)].append(shear)
        shearhistox[int(6*q/qmax)].append(shearox)
    
    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',]
    for edge in G.edges(data="q"):
        x, y, q = edge
        if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
            colors.append(color[int(6*edge[2]/qmax)])       
    
    pos = nx.get_node_attributes(G, 'pos')
    
    edges = []
    qs = []
    for edge in G.edges(data='q'):
            x, y, q = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
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
        

 
    plt.figure(figsize=(20, 20))
    plt.suptitle('Flow for n = '+str(n))
    spec = gridspec.GridSpec(ncols=4, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=4))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=qdrawconst * np.array(qs) / qmax)
    nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')   
    plt.axis('equal')
    
    plt.subplot(spec[4]).set_title('Diameter')
    cindex=0
    plt.xlim((1,1.1*dmax))
    for hist in dhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.axvline(dth, color='k', linestyle='dashed', linewidth=1)
    plt.yscale("log")
    
    plt.subplot(spec[5]).set_title('Flow')
    cindex=0
    plt.xlim((0, 1.1*qmax))
    for hist in qhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[6]).set_title('Shear growth')
    cindex=0
    plt.xlim((0, 1.1*shearmax))
    for hist in shearhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")        
    
    plt.subplot(spec[7]).set_title('Oxygen growth')
    cindex=0
    plt.xlim((0,1.1*shearmaxox))
    for hist in shearhistox:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.savefig(dirname + "/" + name)
    plt.close()


def drawblood(G,name, oxresult, data='q'):
    """
    rysowanie krwi, data to q albo d
    """
    plt.figure(figsize=(20, 20))
    pos = nx.get_node_attributes(G, 'pos')

    qmax = max([edge[2] for edge in G.edges(data=data)])

    edges = []
    qs = []
    for edge in G.edges(data=data):
            x, y, q = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
                edges.append((x, y))
                qs.append(q)

    #draw only those between oxygen nodes:
    for i, edge in enumerate(edges):
        n1, n2 = edge
        if (oxresult[n1] != 1 or oxresult[n2] != 1):
            qs[i] = 0

    nx.draw_networkx_edges(G, pos, edgelist=edges, width=qdrawconst * np.array(qs) / qmax, edge_color='r')

    """
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
    """
    plt.axis('equal')
    plt.xticks([], [])
    plt.yticks([], [])
    plt.savefig(dirname + "/" + name)
    plt.close()
