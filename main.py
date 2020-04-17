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

n = 71 # rozmiar siatki
iters = 301  # liczba iteracji

length_wiggle_param = 1
diameter_wiggle_param = 3

SPARSE = 0 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły


#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

def create_matrix(G, SPARSE=0):
    n = int(len(G)**0.5)
    if (SPARSE == 0): return np.zeros((n * n, n * n))
    if (SPARSE == 1): return spr.csc_matrix(([], ([], [])), shape=(n * n, n * n))

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

def create_pressure_flow_vector2(G, n, qin=1, presout=0):
    global in_nodes, out_nodes
    presult = np.zeros(n * n)
    for node in in_nodes:
        presult[node] = -qin
    for node in out_nodes:
        presult[node] = presout
    return presult

def update_matrix2(G, matrix, SPARSE):
    global in_nodes, out_nodes, reg_nodes, in_edges
    n = int(len(G) ** 0.5)
    def update_nonsparse_matrix(G):
        """
        macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
        z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
        srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
        """
        #matrix = np.zeros((n * n, n * n))

        for node in reg_nodes:
            matrix[node][node] = 0
            for ind in G.neighbors(node):
                d = G[node][ind]["d"]
                l = G[node][ind]['length']
                matrix[node][ind] = c1 * d ** 4 / l
                matrix[node][node] -= c1 * d ** 4 / l

        for node in out_nodes:
            matrix[node][node] = 1

        insert = [0]*len(G)
        for n1, n2 in in_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            insert[n2] += c1*d**4 / l

        sum_insert = np.sum(insert)
        for node in in_nodes:
            matrix[node][:] = insert
            matrix[node][node] -= sum_insert

        return matrix
    def update_sparse_matrix(G):
        """
        macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
        z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
        srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
        """
        row = np.array([])  # numer wiersza w którym będzie wartoć
        col = np.array([])  # numer kolumny w której będzie wartoć
        data = np.array([])  # wartosć

        for node in G.nodes:
            if (node >= n and node < n * (n - 1)):
                row = np.append(row, node)
                col = np.append(col, node)
                data = np.append(data, 0)
                k = data.size - 1
                for ind in G.nodes[node]["neigh"]:
                    d = G[node][ind]["d"]
                    l = G[node][ind]['length']
                    row = np.append(row, node)
                    col = np.append(col, ind)
                    data = np.append(data, c1 * d ** 4 / l)
                    data[k] -= c1 * d ** 4 / l
            elif (node < n):
                row = np.append(row, node)
                col = np.append(col, node)
                data = np.append(data, 0)
                k = data.size - 1
                for ind in range(n):
                    d = G[ind][ind + n]["d"]
                    l = G[ind][ind + n]['length']
                    row = np.append(row, node)
                    col = np.append(col, ind + n)
                    data = np.append(data, c1 * d ** 4 / l)
                    data[k] -= c1 * d ** 4 / l
            else:
                row = np.append(row, node)
                col = np.append(col, node)
                data = np.append(data, 1)
        S = spr.csc_matrix((data, (row, col)), shape=(n * n, n * n))
        return S

    if SPARSE == 1: return update_sparse_matrix(G)
    elif SPARSE == 0: return update_nonsparse_matrix(G)

#=================  IMPORT GEOMETRII  ==================================================================================

G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#from rectangular_geom import *

#from donut_geom import *
#inner_circle, boundary_nodes, boundary_edges = donut_setup(G)

#from cylindrical_geom import *

from donut_geom import *

#=================  NODES  =============================================================================================
# Aby samemu ustawić, trzeba tutaj zadeklarować in_nodes i out_nodes
#in_nodes = [50]
#out_nodes = [5020]
in_nodes, out_nodes = default_nodes(n, G)
reg_nodes = [node for node in G.nodes() if (node not in in_nodes) and (node not in out_nodes)]
in_edges = []
for node in in_nodes:
    for neigh in G.neighbors(node):
        in_edges.append((node, neigh))

#=================  PROGRAM WŁAŚCIWY  ==================================================================================

presult = create_pressure_flow_vector2(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    matrix = update_matrix2(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)
    G = update_graph(G, pnow)

    if i%50 == 0:
        Dr.drawq(G, n, f'{i//50:04d}.png', in_nodes=in_nodes, out_nodes=out_nodes)

