#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec
import pressure as Pr
import vegf as Ve
from build import (G, in_nodes, out_nodes, F_mult_ox, F_mult,
                   c2, dt, dt_ox, dth, boundary_edges, dirname, n)
from config import qdrawconst, ddrawconst
from utils import mu_d

# normalizacja rysowania (maksymalna grubość krawędzi)



def drawq(name, normalize=True, oxdraw=[]):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(20, 20))
    plt.axis('off')
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
        x, y, q = edge
        
        if (x, y) not in boundary_edges and (y, x) not in boundary_edges:            
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
    nx.draw_networkx_nodes(G, pos, node_size = 15, node_color = oxdraw, cmap='plasma', vmax = 10000)
    
    #nx.draw(G, pos, with_labels=True)
    #plt.show()

    plt.axis('equal')
    plt.savefig(dirname + "/" + name)
    plt.close()




def drawd(name, normalize=True, oxdraw = []):
    """
    rysowanie srednic
    """
    plt.figure(figsize=(20, 20))
    plt.axis('off')
    pos = nx.get_node_attributes(G, 'pos')

    if (normalize == False):
        dmax = 1
        drawconst = 1
    else:
        dmax = max([edge[2] for edge in G.edges(data='d')])
        drawconst = ddrawconst

    color = ['white', 'lightgrey', 'darkgrey', 'grey', 'dimgrey', 'black', 'black','black','black','black','black']
    colors = []

    edges = []
    ds = []
    for edge in G.edges(data='d'):
            x, y, d = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
                edges.append((x, y))
                ds.append(d)
                colors.append(color[int(6*edge[2]/dmax)])
                
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color = colors, width=drawconst * np.array(ds) / dmax)

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
    #nx.draw_networkx_nodes(G, pos, node_size=25 * oxdraw, node_color=oxdraw, cmap='Reds')


    plt.axis('equal')
    plt.savefig(dirname + "/" + name)
    plt.close()    
    




def drawhist(name, snow, oxnow=[], oxresult=[], vnow = [], oxdraw = []):
    """
    rysowanie histogramów
    """
    dhist=[[], [], [], [], [], [], [], []]
    qhist=[[], [], [], [], [], [], [], []]
    shearhist=[[], [], [], [], [], [], [], []]
    dshearhist=[[], [], [], [], [], [], [], []]
    vegfgradhist=[[], [], [], [], [], [], [], []]
    dvegfhist=[[], [], [], [], [], [], [], []]
    colors=[]
    dmax=1
    shearmax=0.001
    dshearmax = 0.001
    vegfgradmax=0.001
    dvegfmax = 0.001
    qmax = max([edge[2] for edge in G.edges(data='q')])
    
    for n1, n2 in G.edges():
        q=G[n1][n2]['q']
        d=G[n1][n2]['d']
        F=F_mult*c2 * mu_d(d) *q/d**3
        shear=F
        if shear>shearmax:
            shearmax=shear       
        
    for n1, n2 in G.edges():
        q=G[n1][n2]['q']
        d=G[n1][n2]['d']
        F=F_mult*c2* mu_d(d) *q/d**3
        shear=F
        dshear = Pr.d_update(F)
        if (oxresult[n1] == 1 or oxresult[n2] == 1):    
            F_ox=F_mult_ox*np.abs(vnow[n1] - vnow[n2])
        else:
            F_ox = 0
        vegfgrad=F_ox
        dvegf = Ve.d_update(F_ox)
        if d>dmax:
            dmax=d
        if shear>shearmax:
            shearmax=shear
        if vegfgrad>vegfgradmax:
            vegfgradmax=vegfgrad
        if dshear > dshearmax:
            dshearmax = dshear
        if dvegf > dvegfmax:
            dvegfmax = dvegf
        
        qhist[int(6*q/qmax)].append(q)
        dhist[int(6*q/qmax)].append(d)
        shearhist[int(6*q/qmax)].append(shear)
        dshearhist[int(6*q/qmax)].append(dshear)
        vegfgradhist[int(6*q/qmax)].append(vegfgrad)
        dvegfhist[int(6*q/qmax)].append(dvegf)  
        """
        qhist[int(6*shear/shearmax)].append(q)
        dhist[int(6*shear/shearmax)].append(d)
        if (oxresult[n1] == 1 or oxresult [n2] == 1) and d > dth: 
            shearhist[int(6*shear/shearmax)].append(shear)
        else:
            shearhist[int(6*shear/shearmax)].append(0)
        dshearhist[int(6*shear/shearmax)].append(dshear)
        vegfgradhist[int(6*shear/shearmax)].append(vegfgrad)
        dvegfhist[int(6*shear/shearmax)].append(dvegf)     
        """
    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'k', 'k', 'k', 'k']
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
        

 
    plt.figure(figsize=(25, 20))
    plt.suptitle('Flow for n = '+str(n))
    spec = gridspec.GridSpec(ncols=6, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=6))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=qdrawconst * np.array(qs) / qmax)
    #nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')   
    plt.axis('equal')
    
    plt.subplot(spec[6]).set_title('Diameter')
    cindex=0
    plt.xlim((0,1.1*dmax))
    for hist in dhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.axvline(dth, color='k', linestyle='dashed', linewidth=1)
    plt.yscale("log")
    
    plt.subplot(spec[7]).set_title('Flow')
    cindex=0
    plt.xlim((0, 1.1*qmax))
    for hist in qhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[8]).set_title('Shear')
    cindex=0
    plt.xlim((0, 1.1*shearmax))
    for hist in shearhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")        
    
    plt.subplot(spec[9]).set_title('Shear growth')
    cindex=0
    plt.xlim((-0.01,1.1*dshearmax))
    for hist in dshearhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[10]).set_title('VEGF gradient')
    cindex=0
    plt.xlim((0,1.1*vegfgradmax))
    for hist in vegfgradhist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[11]).set_title('VEGF growth')
    #cindex=0
    #plt.xlim((0,1.1*dvegfmax))
    #for hist in dvegfhist:
    #    if len(hist)>1:
    #        plt.hist(hist, bins=50, color=color[cindex])
    #    cindex+=1
    plt.yscale("log")
    plt.hist(snow, bins=50, color=color[cindex])
    
    plt.savefig(dirname + "/" + name)
    plt.close()


def drawblood(name, oxresult, oxdraw, data='q'):
    """
    rysowanie krwi, data to q albo d
    """
    plt.figure(figsize=(40, 40), dpi = 200)
    plt.axis('off')
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
    nx.draw_networkx_nodes(G, pos, node_size = 30, node_color = oxdraw, cmap='Blues')
    
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

def draw_initial_d():
    ds = []
    for edge in G.edges(data='d'):
            x, y, d = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
                ds.append(d)
    plt.hist(ds, bins = 100, density = True)
    plt.xlabel('diameter')
    plt.ylabel('count (normalized)')
    plt.title('Initial diameters distribution')
    
def plot_params(prestab, oxtab, vegftab):
    plt.figure(figsize=(15, 5))
    plt.suptitle('Parameters')
    spec = gridspec.GridSpec(ncols=3, nrows=1)
    
    plt.subplot(spec[0]).set_title('Pressure')
    x = np.linspace(0, len(prestab), len(prestab))
    plt.plot(x, prestab)
    plt.xlabel('iterations')
    plt.ylabel(r'$\Delta p$')
    
    plt.subplot(spec[1]).set_title('Oxygen')
    x = np.linspace(0, len(prestab), len(prestab))
    plt.plot(x, oxtab)
    plt.xlabel('iterations')
    plt.ylabel(r'average $c_{ox}^2$')
    
    plt.subplot(spec[2]).set_title('VEGF')
    x = np.linspace(0, len(prestab), len(prestab))
    plt.plot(x, vegftab)
    plt.xlabel('iterations')
    plt.ylabel(r'average $c_{vegf}^2$')
    
    plt.savefig(dirname +'/parameters.png')
    plt.close()
    