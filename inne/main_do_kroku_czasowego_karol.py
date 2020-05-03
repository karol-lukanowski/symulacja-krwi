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


#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 151 # rozmiar siatki
nkw=n**2

length_wiggle_param = 0
diameter_wiggle_param = 3

SPARSE = 1 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

F0=0.2
F1=1
z0=0
z1=1
F_mult = 10000
dt = 0.1

iters = int(50/dt)+1  # liczba iteracji


#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

def create_matrix(G, SPARSE=0):
    def create_sparse_matrix(G):
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
        to_sort = sorted(to_sort, key=lambda elem: elem[1] * nkw + elem[2])
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
        for n1, n2, d, l in in_edges:
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


        diag = np.zeros(nkw)
        for n1, n2, d, l in reg_reg_edges:
            res = c1 * d ** 4 / l
            data.append((res, n1 * nkw + n2))
            data.append((res, n2 * nkw + n1))
            diag[n1] -= res
            diag[n2] -= res

        insert = defaultdict(float)

        for n1, n2, d, l in reg_something_edges:
            res = c1 * d ** 4 / l
            data.append((res, n1 * nkw + n2))
            diag[n1] -= res
            



        for node, datum in enumerate(diag):
            if datum != 0:
                data.append((datum, node * nkw + node))

        for node in out_nodes:
            data.append((1, node * nkw + node))


        for n1, n2, d, l in in_edges:
            insert[n2] += c1 * d ** 4 / l

        sum_insert = sum(insert.values())

        for node in in_nodes:
            data.append((-sum_insert, node * nkw + node))

            for ins_node, ins in insert.items():
                data.append((ins, node * nkw + ins_node))                

        data.sort(key=lambda elem: elem[1])

        global indices, indptr
        return spr.csr_matrix(([datum[0] for datum in data], indices, indptr), shape=(nkw, nkw))


    if SPARSE == 1: return update_sparse_matrix(G)
    elif SPARSE == 0: return update_nonsparse_matrix(G)

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


#=================  GRAF I GEOMETRIA  ==================================================================================

#G = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#nx.write_edgelist(G,'test.edgelist',data=['d'])
example_graph = nx.read_edgelist('test.edgelist', nodetype=int, data=(('d', float),))
# wczytaj zapisane srednice
for n1, n2 in G.edges():
    G[n1][n2]['d'] = example_graph[n1][n2]['d']


in_nodes, out_nodes = equidistant_geometry(R = 50, xrange = n*np.sqrt(3)/2, yrange = n, how_many = 12)

#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='cylindrical', R=n//2.5)
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=n//2.5, R_s=n//20)
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='rect')
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)


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
matrix = create_matrix(G, SPARSE=SPARSE)
if SPARSE: indices, indptr = matrix.indices, matrix.indptr

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')
    matrix = update_matrix(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)
    if i%((iters-1)//5) == 0:
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
       
        Dr.drawq(G, n, f'd{i//((iters-1)//5):04d}.png', in_nodes=in_nodes, out_nodes=out_nodes)
    reg_reg_edges, reg_something_edges, in_edges=update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)

for n1, n2 in G.edges():
    print(n1, n2, G[n1][n2]['q'])