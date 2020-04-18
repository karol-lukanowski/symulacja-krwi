#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import triangular_net as Tr
import draw_net as Dr
from geometry import set_geometry

#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 201 # rozmiar siatki
iters = 501  # liczba iteracji

length_wiggle_param = 1
diameter_wiggle_param = 3

SPARSE = 1 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły


#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

def create_matrix(G, SPARSE=0):
    def create_sparse_matrix(G):
        data, row, col = [], [], []
        for node in reg_nodes:
            data_for_this_node = 0
            for neigh in G.neighbors(node):
                d = G[node][neigh]['d']
                l = G[node][neigh]['length']
                row.append(node)
                col.append(neigh)
                data.append(c1 * d ** 4 / l)
                data_for_this_node -= c1 * d ** 4 / l
            row.append(node)
            col.append(node)
            data.append(data_for_this_node)

        for node in out_nodes:
            row.append(node)
            col.append(node)
            data.append(1)

        insert = [0] * len(G)
        for n1, n2 in in_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            insert[n2] += c1 * d ** 4 / l

        sum_insert = np.sum(insert)

        for node in in_nodes:
            for i in range(len(G)):
                if i == node or insert[i] == 0: continue
                row.append(node)
                col.append(i)
                data.append(insert[i])

            row.append(node)
            col.append(node)
            data.append(-sum_insert)

        # posortujmy teraz dane tak, aby były najpierw po row, potem po column
        to_sort = list(zip(data, row, col))
        to_sort = sorted(to_sort, key=lambda elem: elem[1] * n ** 2 + elem[2])
        data, row, col = zip(*to_sort)

        return spr.csr_matrix((data, (row, col)), shape=(n * n, n * n))
    if (SPARSE == 0): return np.zeros((n * n, n * n))
    if (SPARSE == 1): return create_sparse_matrix(G)

def solve_equation_for_pressure(matrix, presult):
    """
    Zamieniamy macierz w równaniu na formę macierzy rzadkiej
    w celu usprawnienia obliczeń
    """
    if (SPARSE == 0): pnow = sprlin.spsolve(spr.csc_matrix(matrix), presult)
    elif (SPARSE == 1): pnow = sprlin.spsolve(matrix, presult)

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

def create_pressure_flow_vector(G, n, qin=1, presout=0):
    global in_nodes, out_nodes
    presult = np.zeros(n * n)
    for node in in_nodes:
        presult[node] = -qin
    for node in out_nodes:
        presult[node] = presout
    return presult

def update_matrix(G, matrix, SPARSE):
    global in_nodes, out_nodes, reg_nodes, in_edges
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
        data = []
        pos = []
        for node in reg_nodes:
            data_for_this_node = 0
            for neigh in G.neighbors(node):
                d = G[node][neigh]['d']
                l = G[node][neigh]['length']
                data.append(c1 * d ** 4 / l)
                pos.append(node * n ** 2 + neigh)
                data_for_this_node -= c1 * d ** 4 / l
            data.append(data_for_this_node)
            pos.append(node * n ** 2 + node)

        for node in out_nodes:
            data.append(1)
            pos.append(node * n ** 2 + node)

        insert = [0] * len(G)
        for n1, n2 in in_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            insert[n2] += c1 * d ** 4 / l

        sum_insert = np.sum(insert)

        for node in in_nodes:
            for i in range(len(G)):
                if i == node or insert[i] == 0: continue
                data.append(insert[i])
                pos.append(node * n ** 2 + i)

            data.append(-sum_insert)
            pos.append(node * n ** 2 + node)

        to_sort = zip(data, pos)
        to_sort = sorted(to_sort, key=lambda elem: elem[1])
        data, pos = zip(*to_sort)

        global indices, indptr
        return spr.csr_matrix((data, indices, indptr), shape=(n * n, n * n))

    if SPARSE == 1: return update_sparse_matrix(G)
    elif SPARSE == 0: return update_nonsparse_matrix(G)

#=================  GRAF I GEOMETRIA  ==================================================================================

G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='cylindrical', R=80)

#=================  PROGRAM WŁAŚCIWY  ==================================================================================

presult = create_pressure_flow_vector(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)
if SPARSE: indices, indptr = matrix.indices, matrix.indptr

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    matrix = update_matrix(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)
    G = update_graph(G, pnow)

    print(pnow[in_nodes])

    if i%100 == 0:
        Dr.drawq(G, n, f'{i//50:04d}.png', in_nodes=in_nodes, out_nodes=out_nodes)
