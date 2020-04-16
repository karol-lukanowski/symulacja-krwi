#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import triangular_net as Tr
import draw_net as Dr

#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 101  # rozmiar siatki
iters = 205  # liczba iteracji

length_wiggle_param = 1
diameter_wiggle_param = 0

SPARSE = 0 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

def solve_equation_for_pressure(matrix, presult):
    """
    Zamieniamy macierz w równaniu na formę macierzy rzadkiej
    w celu usprawnienia obliczeń
    """
    if (SPARSE == 0): pnow = sprlin.spsolve(spr.csc_matrix(matrix), presult)
    elif (SPARSE == 1): pnow = sprlin.spsolve(matrix, presult)

    # jeżeli jest geometria donutowa to robimy fix
    try: pnow = donut_fix(inner_circle, boundary_nodes, pnow)
    except NameError:
        pass

    return pnow

def update_graph(G, pnow):
    def d_update(F):
        #zmiana średnicy pod względem siły F
        if F < 0.00002:
            return 0
        else:
            if F < 0.002:
                return 1000 * F
            else:
                return 0.2
    lengths = np.fromiter(nx.get_edge_attributes(G, 'length').values(), dtype=float)

    # ustawienie nowych ciśnień w węzłach
    pnow_dict = {node: p for node, p in enumerate(pnow)}
    nx.set_node_attributes(G, name='p', values=pnow_dict)

    # obliczenie nowych przepływów w krawędziach
    ds = np.fromiter(nx.get_edge_attributes(G, 'd').values(), dtype=float)
    dps = np.array([(pnow[n1] - pnow[n2]) for (n1, n2) in G.edges()])
    new_qs = c1 * ds ** 4 * dps / lengths

    # obliczenie nowych średnic krawędzi
    fs = c2 * np.abs(new_qs) / ds ** 3
    new_ds = [d + d_update(F) for d, F in zip(ds, fs)]

    # ustawienie nowych atrybutów krawędzi
    new_values = {edge: {'q': new_qs[i], 'd': new_ds[i]} for i, edge in enumerate(G.edges())}
    nx.set_edge_attributes(G, new_values)

    return G

def r_squared(node):
    x0, y0 = G.nodes[n*n//2]["pos"]
    x, y = G.nodes[node]['pos']
    r_sqr = (x - x0) ** 2 + (y - y0) ** 2
    return r_sqr

#=================  IMPORT GEOMETRII  ==================================================================================

G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#from rectangular_geom import *

from cylindrical_geom import *

#from donut_geom import *
#inner_circle, boundary_nodes, boundary_edges = donut_setup(G)

#=================  PROGRAM WŁAŚCIWY  ==================================================================================

ds = []
lens = []
for (n1, n2) in G.edges():
    ds.append(G[n1][n2]['d'])
    lens.append(G[n1][n2]['length'])



presult = create_pressure_flow_vector(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)

matrix = update_matrix(G, matrix, SPARSE)
pnow = solve_equation_for_pressure(matrix, presult)
G = update_graph(G, pnow)

q_dict = nx.get_edge_attributes(G, 'q')

for key in q_dict:
    if q_dict[key] != 0:
        q_dict[key] = 1

nx.set_edge_attributes(G, q_dict, 'q')

pos = nx.get_node_attributes(G, 'pos')

edges = []
widths = []
for edge in G.edges(data='q'):
    if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
        edges.append(edge)
        widths.append(edge[2]/2)
id_center = n**2 // 2

### boundary nodes calculation and sorting
boundary_nodes = []
for (n1, n2) in G.edges():
    r1, r2 = r_squared(n1), r_squared(n2)
    if r1 > r2:
        n1, n2 = n2, n1
        r1, r2 = r2, r1

    n_b=n2

    if r2 >= R**2 and r1 <= R**2:
        x, y = G.nodes[n_b]['pos'][0] - G.nodes[n**2 // 2]['pos'][0], G.nodes[n_b]['pos'][1] - G.nodes[n**2 // 2]['pos'][1]

        if (x>=0 and y>=0):
            fi = np.arctan(y/x)
        elif (x<0 and y>=0):
            fi = np.pi/2 + np.arctan(-x / y)
        elif (x<0 and y<0):
            fi = np.pi + np.arctan(y/x)
        else:
            fi = (3/2)*np.pi + np.arctan(x / -y)
        boundary_nodes.append([n_b, fi])
boundary_nodes.sort(key=lambda node: node[1])

boundary_nodes, fis = [node[0] for node in boundary_nodes], [node[1] for node in boundary_nodes]
boundary_nodes = list(dict.fromkeys(boundary_nodes))
fis = list(dict.fromkeys(fis))

###


ls = []


for target_num in range(len(boundary_nodes)):
    l = nx.shortest_path_length(G, source=id_center, target=boundary_nodes[target_num], weight='length')
    ls.append(l)
    path = nx.shortest_path(G, source=id_center, target=boundary_nodes[target_num], weight='length')

    print(target_num, l)

    fig = plt.figure(figsize=(7, 7))
    fig.suptitle(f'Siatka {n}x{n} węzłów, szum o sile {length_wiggle_param}')

    nx.draw_networkx_edges(G, pos, edgelist=edges, width=0.5*np.array(widths))
    nx.draw_networkx_nodes(G, pos, nodelist=[id_center], node_size=3, node_color='r')
    nx.draw_networkx_nodes(G, pos, nodelist=boundary_nodes, node_size=3, node_color='g')
    nx.draw_networkx_nodes(G, pos, nodelist=path, node_size=3, node_color='b')
    plt.axis('off')


    #ax = plt.subplot(111, projection='polar')
    #c = plt.plot(fis, 20 * np.ones_like(fis) + test)
    plt.axis('equal')

    plt.gcf().add_subplot(121, projection='polar')
    ax = plt.gca()

    r_small = 41
    r_big = 45

    ax.plot(fis[:len(ls)], ls, marker='.', markersize=0.8, linewidth=0.3)
    ax.plot(fis, r_big*np.ones_like(fis), linestyle='--', linewidth=0.5)
    ax.plot(fis, r_small*np.ones_like(fis), linestyle='--', linewidth=0.5)

    plt.text(0, r_small, str(r_small), fontsize=7, horizontalalignment='center', verticalalignment='bottom')
    plt.text(0, r_big, str(r_big), fontsize=7, horizontalalignment='center', verticalalignment='bottom')

    box = ax.get_position()
    box.x1 += 0.4215



    move_ax = 0.25
    box.x0 -= move_ax
    box.y0 -= move_ax
    box.x1 += move_ax
    box.y1 += move_ax


    ax.set_position(box)




    ax.patch.set_alpha(0)
    ax.axis('off')




    #plt.show()
    plt.savefig(f'{target_num:04d}.png')
    plt.close()

print(min(ls), max(ls))
