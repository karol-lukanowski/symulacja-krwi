import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec

import pressure as Pr
import vegf as Ve

from utils import mu_d, d_update

from config import simInputData


# normalizacja rysowania (maksymalna grubość krawędzi)

def draw(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, oxresult, name, data='q'):
    """
    rysowanie krwi, data to q albo d
    """
    plt.figure(figsize=(40, 40), dpi = 200)
    plt.axis('off')
    pos = nx.get_node_attributes(G, 'pos')

    qmax = max([edge[2] for edge in G.edges(data=data)])

    edges = []
    qs = []
    colors = []
    for edge in G.edges(data=data):
            x, y, q = edge
            if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
                edges.append((x, y))
                qs.append(q)
                if oxresult[x] == 1 and oxresult[y] == 1:
                    colors.append('r')
                else:
                    colors.append('k')


    nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.qdrawconst * np.array(qs) / qmax, edge_color = colors)
    #nx.draw_networkx_nodes(G, pos, node_size = 30, node_color = oxdraw, cmap='Blues')
    
    
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


def drawblood(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, name, oxresult, oxdraw, data='q'):
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

    nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.qdrawconst * np.array(qs) / qmax, edge_color='r')
    nx.draw_networkx_nodes(G, pos, node_size = 30, node_color = oxdraw, cmap='Blues')
    
    
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


def drawhist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, vnow, snow_upstream, snow_downstream, name):
    d_hist = [[], [], [], [], [], [], [], []]
    q_hist = [[], [], [], [], [], [], [], []]
    shear_hist = [[], [], [], [], [], [], [], []]
    vegf_hist = [[], [], [], [], [], [], [], []]
    gradp_hist = [[], [], [], [], [], [], [], []]
    upstream_hist = [[], [], [], [], [], [], [], []]
    downstream_hist = [[], [], [], [], [], [], [], []]
    colors = []
    d_max= 1
    shear_max=0.001
    vegf_max = 0.001
    gradp_max = 0.001
    upstream_max = 0.001
    downstream_max = 0.001
    q_max = max([edge[2] for edge in G.edges(data='q')])
    pin = np.max(pnow)

    for n1, n2 in G.edges():
        q = G[n1][n2]['q']
        d = G[n1][n2]['d']
        l = G[n1][n2]['length']

        shear = sid.F_mult * sid.c2 * mu_d(d) * q / d ** 3

        if (oxresult[n1] == 1 or oxresult[n2] == 1):    
            vegf = sid.F_mult_ox * np.abs(vnow[n1] - vnow[n2])
            gradp = sid.cp * np.abs(pnow[n1] - pnow[n2]) / (pin * l)
        else:
            vegf = 0
            gradp = 0

        if oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] > pnow[n2]:
                upstream = sid.cs * snow_upstream[n1]
                downstream = sid.cs * snow_downstream[n2]
            else:
                upstream = sid.cs * snow_upstream[n2]
                downstream = sid.cs * snow_downstream[n1]
        else:
            upstream = 0
            downstream = 0
        
        if d > d_max:
            d_max = d
        if shear > shear_max:
            shear_max = shear
        if vegf > vegf_max:
            vegf_max = vegf
        if gradp > gradp_max:
            gradp_max = gradp
        if upstream > upstream_max:
            upstream_max = upstream
        if downstream > downstream_max:
            downstream_max = downstream

        q_hist[int(6 * q / q_max)].append(q)
        d_hist[int(6 * q / q_max)].append(d)
        shear_hist[int(6 * q / q_max)].append(shear)
        vegf_hist[int(6 * q / q_max)].append(vegf)
        gradp_hist[int(6 * q / q_max)].append(gradp)
        upstream_hist[int(6 * q / q_max)].append(upstream)
        downstream_hist[int(6 * q / q_max)].append(downstream)

    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'k', 'k', 'k', 'k']

    edges = []
    qs = []
    for edge in G.edges(data="q"):
        x, y, q = edge
        if (x, y) not in boundary_edges and (y, x) not in boundary_edges:
            colors.append(color[int(6 * edge[2] / q_max)])
            edges.append((x, y))
            qs.append(q)   

    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
        

 
    plt.figure(figsize=(25, 20))
    spec = gridspec.GridSpec(ncols=7, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=7))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=sid.qdrawconst * np.array(qs) / q_max)
    #nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')   
    plt.axis('equal')
    
    plt.subplot(spec[7]).set_title('Diameter')
    cindex=0
    plt.xlim(0,1.1*d_max)
    for hist in d_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.axvline(sid.dth, color='k', linestyle='dashed', linewidth=1)
    plt.yscale("log")
    
    plt.subplot(spec[8]).set_title('Flow')
    cindex=0
    plt.xlim(0, 1.1*q_max)
    for hist in q_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[9]).set_title('Shear')
    cindex=0
    plt.xlim(0, 1.1*shear_max)
    for hist in shear_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")        
    
    plt.subplot(spec[10]).set_title('VEGF')
    cindex=0
    plt.xlim(0,1.1*vegf_max)
    for hist in vegf_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[11]).set_title('Pressure gradient')
    cindex=0
    plt.xlim(0,1.1*gradp_max)
    for hist in gradp_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    
    plt.subplot(spec[12]).set_title('Upstream')
    cindex=0
    plt.xlim(0,1.1*upstream_max)
    for hist in upstream_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    plt.hist(vnow, bins=50, color=color[cindex])

    plt.subplot(spec[13]).set_title('Downstream')
    cindex=0
    plt.xlim(0,1.1*downstream_max)
    for hist in downstream_hist:
        if len(hist)>1:
            plt.hist(hist, bins=50, color=color[cindex])
        cindex+=1
    plt.yscale("log")
    plt.hist(vnow, bins=50, color=color[cindex])
    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()


def uniform_hist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, vnow, snow_upstream, snow_downstream, name):
    d_hist = []
    q_hist = []
    shear_hist = []
    vegf_hist = []
    gradp_hist = []
    upstream_hist = []
    downstream_hist = []
    d_shear_hist = []
    d_vegf_hist = []
    d_gradp_hist = []


    for n1, n2 in G.edges():
        q = G[n1][n2]['q']
        d = G[n1][n2]['d']
        l = G[n1][n2]['length']

        shear = sid.F_mult * sid.c2 * mu_d(d) * q / d ** 3
        d_shear = d_update(shear, sid.F_p)

        if (oxresult[n1] == 1 or oxresult[n2] == 1):    
            vegf = sid.F_mult_ox * np.abs(vnow[n1] - vnow[n2])
            gradp = sid.cp * np.abs(pnow[n1] - pnow[n2]) / l
        else:
            vegf = 0
            gradp = 0

        d_vegf = d_update(vegf, sid.F_ox)
        d_gradp = d_update(gradp, sid.F_gradp)

        if oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] > pnow[n2]:
                upstream = sid.cs * snow_upstream[n1]
                downstream = sid.cs * snow_downstream[n2]
            else:
                upstream = sid.cs * snow_upstream[n2]
                downstream = sid.cs * snow_downstream[n1]
        else:
            upstream = 0
            downstream = 0
        
        d_upstream = d_update(upstream, sid.F_s)
        d_downstream = d_update(downstream, sid.F_s)


        q_hist.append(q)
        d_hist.append(d)
        shear_hist.append(shear)
        vegf_hist.append(vegf)
        gradp_hist.append(gradp)
        upstream_hist.append(d_upstream)
        downstream_hist.append(d_downstream)
        d_shear_hist.append(d_shear)
        d_vegf_hist.append(d_vegf)
        d_gradp_hist.append(d_gradp)

    edges = []
    qs = []
    colors = []
    for n1, n2 in G.edges():
        q = G[n1][n2]['q']
        d = G[n1][n2]['d']
        l = G[n1][n2]['length']
        if (n1, n2) not in boundary_edges and (n2, n1) not in boundary_edges:
            edges.append((n1, n2))
            qs.append(q)
            F = sid.F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
            colors.append(d_update(F, sid.F_p))

    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
        
    q_max = max([edge[2] for edge in G.edges(data='q')])
 
    plt.figure(figsize=(25, 20))
    spec = gridspec.GridSpec(ncols=5, nrows=3, height_ratios=[5, 1, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=5))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color = colors, width=sid.qdrawconst * np.array(qs) / q_max, edge_cmap = plt.cm.plasma)
    #nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')   
    plt.axis('equal')
    
    plt.subplot(spec[5]).set_title('Diameter')
    plt.hist(d_hist, bins=50)
    plt.axvline(sid.dth, color='k', linestyle='dashed', linewidth=1)
    plt.yscale("log")

    plt.subplot(spec[6]).set_title('Flow')
    plt.hist(q_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[7]).set_title('Shear')
    plt.hist(shear_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[8]).set_title('VEGF')
    plt.hist(vegf_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[9]).set_title('Grad p')
    plt.hist(gradp_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[10]).set_title('Upstream growth')
    plt.hist(upstream_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[11]).set_title('Signal')
    plt.hist(snow_upstream, bins=50)
    plt.yscale("log")

    plt.subplot(spec[12]).set_title('Shear growth')
    plt.hist(d_shear_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[13]).set_title('VEGF growth')
    plt.hist(d_vegf_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[14]).set_title('Grad p growth')
    plt.hist(d_gradp_hist, bins=50)
    plt.yscale("log")    

    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()