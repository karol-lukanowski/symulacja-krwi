#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import triangular_net as Tr
import draw_net as Dr
from collections import defaultdict
from geometry import set_geometry

#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 51 # rozmiar siatki
iters = 50  # liczba iteracji

length_wiggle_param = 1
diameter_wiggle_param = 3

SPARSE = 1 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły


#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

def create_matrix(G, SPARSE=0):
    def create_sparse_matrix(G):
        data, row, col = [], [], []

        diag = np.zeros(n * n)
        for n1, n2 in reg_reg_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            res = c1 * d ** 4 / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n1] -= res
            diag[n2] -= res
        for n1, n2 in reg_something_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
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
        for n1, n2 in in_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            insert[n2] += c1 * d**4 / l

        sum_insert = sum(insert.values())

        for node in in_nodes:
            data.append(-sum_insert)
            row.append(node)
            col.append(node)

            for ins_node, ins in insert.items():
                data.append(ins)
                row.append(node)
                col.append(ins_node)


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

        diag = np.zeros(n * n)
        for n1, n2 in reg_reg_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            res = c1 * d ** 4 / l
            data.append(res)
            pos.append(n1 * n ** 2 + n2)
            data.append(res)
            pos.append(n2 * n ** 2 + n1)
            diag[n1] -= res
            diag[n2] -= res
        for n1, n2 in reg_something_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            res = c1 * d ** 4 / l
            data.append(res)
            pos.append(n1 * n ** 2 + n2)
            diag[n1] -= res

        for node, datum in enumerate(diag):
            if datum != 0:
                pos.append(node * n ** 2 + node)
                data.append(datum)

        for node in out_nodes:
            data.append(1)
            pos.append(node * n ** 2 + node)

        insert = defaultdict(float)
        for n1, n2 in in_edges:
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            insert[n2] += c1 * d ** 4 / l

        sum_insert = sum(insert.values())

        for node in in_nodes:
            data.append(-sum_insert)
            pos.append(node * n ** 2 + node)

            for ins_node, ins in insert.items():
                data.append(ins)
                pos.append(node * n ** 2 + ins_node)

        to_sort = zip(data, pos)
        to_sort = sorted(to_sort, key=lambda elem: elem[1])
        data, pos = zip(*to_sort)

        global indices, indptr
        return spr.csr_matrix((data, indices, indptr), shape=(n * n, n * n))

    if SPARSE == 1: return update_sparse_matrix(G)
    elif SPARSE == 0: return update_nonsparse_matrix(G)

#=================  GRAF I GEOMETRIA  ==================================================================================

G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=20, R_s=5)
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=[15, 25, 30], out_nodes=[950, 903])


reg_reg_edges, reg_something_edges = [],  []
for n1, n2 in G.edges():
    if (n1 not in in_nodes and n1 not in out_nodes) and (n2 not in in_nodes and n2 not in out_nodes):
        reg_reg_edges.append((n1, n2))
    elif (n1 not in in_nodes and n1 not in out_nodes):
        reg_something_edges.append((n1, n2))
    elif (n2 not in in_nodes and n2 not in out_nodes):
        reg_something_edges.append((n2, n1))



#=================  PROGRAM WŁAŚCIWY  ==================================================================================

presult = create_pressure_flow_vector(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)
if SPARSE: indices, indptr = matrix.indices, matrix.indptr

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    matrix = update_matrix(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)
    G = update_graph(G, pnow)

    #if i%50 == 0:
    #    Dr.drawq(G, n, f'{i//50:04d}.png', in_nodes=in_nodes, out_nodes=out_nodes)
