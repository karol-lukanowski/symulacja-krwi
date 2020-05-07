#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import triangular_net as Tr
import draw_net as Dr
import delaunay as De
from collections import defaultdict
from geometry import set_geometry
import time

#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 101 # rozmiar siatki
nkw=n**2

length_wiggle_param = 1
diameter_wiggle_param = 3

SPARSE = 1 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

D = 0.138 # współczynnik dyfuzji
k = 0.137 # stała reakcji
dth = 10 # graniczna grubosć

F0=0.2
F1=1
z0=0
z1=1
F_mult = 10000
dt = 0.8

iters = 301  # liczba iteracji
save_every = 30

#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

def solve_equation_for_pressure(matrix, presult):
    """
    Zamieniamy macierz w równaniu na formę macierzy rzadkiej
    w celu usprawnienia obliczeń
    """
    if (SPARSE == 0): pnow = sprlin.spsolve(spr.csc_matrix(matrix), presult)
    elif (SPARSE == 1): pnow = sprlin.spsolve(matrix, presult)

    return pnow

def update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges):
    def d_update(F):
        #zmiana średnicy pod względem siły F
        result = 0
        if (F > F0):
            if (F < F1):
                result = z0+(F-F0)*(z1-z0)/(F1-F0)
            else:
                result = z1
        else:
            result = z0
        return result * dt


    reg_reg_edges2, reg_something_edges2, in_edges2=[], [], []
    for n1, n2, d, l in reg_reg_edges:
        F=F_mult*c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew=d+d_update(F)
        reg_reg_edges2.append((n1, n2, dnew, l))

    for n1, n2, d, l in reg_something_edges:
        F=F_mult*c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew=d+d_update(F)
        reg_something_edges2.append((n1, n2, dnew, l))
    for n1, n2, d, l in in_edges:
        F=F_mult*c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew=d+d_update(F)
        in_edges2.append((n1, n2, dnew, l))

    return reg_reg_edges2, reg_something_edges2, in_edges2

def create_pressure_flow_vector(G, n, qin=1, presout=0):
    global in_nodes, out_nodes
    presult = np.zeros(n * n)
    for node in in_nodes:
        presult[node] = -qin
    for node in out_nodes:
        presult[node] = presout
    return presult

def update_matrix():
    global in_nodes, out_nodes, reg_nodes, in_edges
    data, row, col = [], [], []

    diag = np.zeros(n * n)
    for n1, n2, d, l in reg_reg_edges:
        res = c1 * d ** 4 / l
        data.append(res)
        row.append(n1)
        col.append(n2)
        data.append(res)
        row.append(n2)
        col.append(n1)
        diag[n1] -= res
        diag[n2] -= res
    for n1, n2, d, l in reg_something_edges:
        res = c1 * d ** 4 / l
        data.append(res)
        row.append(n1)
        col.append(n2)
        diag[n1] -= res
    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
    for node in out_nodes:
        row.append(node)
        col.append(node)
        data.append(1)

    insert = defaultdict(float)
    for n1, n2, d, l in in_edges:
        insert[n2] += c1 * d ** 4 / l
    sum_insert = sum(insert.values())

    for node in in_nodes:
        data.append(-sum_insert)
        row.append(node)
        col.append(node)

        for ins_node, ins in insert.items():
            data.append(ins)
            row.append(node)
            col.append(ins_node)

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw))

def equidistant_geometry(R, xrange, yrange, how_many):
    id_center = De.find_center_node(G, n, xrange=xrange, yrange=yrange)

    def r_squared(node):
        # x0, y0 = G.nodes[n*n//2]["pos"]
        x0, y0 = G.nodes[id_center]['pos']
        x, y = G.nodes[node]['pos']
        r_sqr = (x - x0) ** 2 + (y - y0) ** 2
        return r_sqr

    boundary_nodes = []
    for (n1, n2) in G.edges():
        r1, r2 = r_squared(n1), r_squared(n2)
        if r1 > r2:
            n1, n2 = n2, n1
            r1, r2 = r2, r1

        n_b = n2

        if r2 >= R ** 2 and r1 <= R ** 2:
            # x, y = G.nodes[n_b]['pos'][0] - G.nodes[n**2 // 2]['pos'][0], G.nodes[n_b]['pos'][1] - G.nodes[n**2 // 2]['pos'][1]
            x, y = G.nodes[n_b]['pos'][0] - G.nodes[id_center]['pos'][0], G.nodes[n_b]['pos'][1] - \
                   G.nodes[id_center]['pos'][1]

            if x == 0: x = 0.000001
            if y == 0: y = 0.000001

            if (x >= 0 and y >= 0):
                fi = np.arctan(y / x)
            elif (x < 0 and y >= 0):
                fi = np.pi / 2 + np.arctan(-x / y)
            elif (x < 0 and y < 0):
                fi = np.pi + np.arctan(y / x)
            else:
                fi = (3 / 2) * np.pi + np.arctan(x / -y)
            boundary_nodes.append([n_b, fi])
    boundary_nodes.sort(key=lambda node: node[1])

    boundary_nodes, fis = zip(*boundary_nodes)

    num_of_out_nodes = how_many
    out_indexes = np.round(np.linspace(0, len(boundary_nodes) - 1, num_of_out_nodes + 1)).astype(int)
    out_nodes = list(np.array(boundary_nodes)[out_indexes[:-1]])
    in_nodes = [id_center]

    return in_nodes, out_nodes


#=================  FUNKCJE TLEN  ======================================================================================


def create_ox_vector(G, n):
    global in_nodes_ox, out_nodes_ox
    oxresult = np.zeros(nkw)
    for node in in_nodes_ox:
        oxresult[node] = 1
    for node in out_nodes_ox:
        oxresult[node] = 1
    return oxresult

def update_matrix_ox(G):

    matrix = np.zeros((n * n, n * n))

    for node in reg_nodes:
        matrix[node][node] = k
        for ind in G.neighbors(node):
            d = G[node][ind]["d"]
            l = G[node][ind]['length']
            if d > dth:
                matrix[node][node] = 1
                oxresult[node] = 1
                oxresult[ind] = 1
            else:
                matrix[node][ind] = D / l
                matrix[node][node] -= D / l

    for node in out_nodes:
        matrix[node][node] = 1

    for node in in_nodes:
        matrix[node][node] = 1

    return matrix

def update_sparse_matrix_ox(G):

    data, row, col = [], [], []

    for node in reg_nodes:
        this_node = k
        for ind in G.neighbors(node):
            d = G[node][ind]["d"]
            l = G[node][ind]['length']
            if d > dth:
                this_node = 1
                oxresult[node] = 1
                oxresult[ind] = 1
            else:
                data.append(D/l)
                row.append(node)
                col.append(ind)

                this_node -= D/l

        data.append(this_node)
        row.append(node)
        col.append(node)

    for node in in_nodes + out_nodes:
        data.append(1)
        row.append(node)
        col.append(node)

    return spr.csr_matrix((data, (row, col)), shape=(nkw, nkw))



#=================  GRAF I GEOMETRIA  ==================================================================================

#G = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#in_nodes, out_nodes = equidistant_geometry(R = n//2.5, xrange = n, yrange = n, how_many = 100)

#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='cylindrical', R=n//2.5)
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=n//2.5, R_s=n//20)
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='rect')
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)

in_nodes_ox, out_nodes_ox, reg_nodes_ox = in_nodes.copy(), out_nodes.copy(), reg_nodes.copy()

reg_reg_edges, reg_something_edges = [],  []
for n1, n2 in G.edges():
    d = G[n1][n2]['d']
    l = G[n1][n2]['length']
    if (n1 not in in_nodes and n1 not in out_nodes) and (n2 not in in_nodes and n2 not in out_nodes):
        reg_reg_edges.append((n1, n2, d, l))
    elif (n1 not in in_nodes and n1 not in out_nodes):
        reg_something_edges.append((n1, n2, d, l))
    elif (n2 not in in_nodes and n2 not in out_nodes):
        reg_something_edges.append((n2, n1, d, l ))

#=================  PROGRAM WŁAŚCIWY  ==================================================================================

presult = create_pressure_flow_vector(G, n, qin, presout)

oxresult = create_ox_vector(G, n)

#matrix_ox = np.zeros((n * n, n * n))


for i in range(iters):
    print(f'Iter {i + 1}/{iters}')
    matrix = update_matrix()

    matrix_ox = update_sparse_matrix_ox(G)

    pnow = solve_equation_for_pressure(matrix, presult)

    oxnow = sprlin.spsolve(matrix_ox, oxresult)


    if i%save_every == 0:
        Q_in = 0
        Q_out = 0
        for n1, n2, d, l in reg_reg_edges:
            G[n1][n2]['d']=d
            q=c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q']=q

        for n1, n2, d, l in reg_something_edges:        
            G[n1][n2]['d']=d
            q=c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q']=q
            
            if n2 in in_nodes:
                Q_in+=q
            if n2 in out_nodes:
                Q_out+=q
        
        print('Q_in =', Q_in, 'Q_out =', Q_out)

        #nx.draw_networkx_nodes(G, pos, node_size=100 * oxresult, node_color='b')
        Dr.drawq(G, n, f'{i//save_every:04d}.png', in_nodes=in_nodes, out_nodes=out_nodes, oxresult = oxresult)
    reg_reg_edges, reg_something_edges, in_edges=update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)
