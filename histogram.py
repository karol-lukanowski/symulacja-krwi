#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin
import triangular_net as Tr
import delaunay as De
import draw_net as Dr
from collections import defaultdict
from geometry import set_geometry
from matplotlib import gridspec

#=================  WARUNKI POCZĄTKOWE i STAŁE  ========================================================================

n = 101 # rozmiar siatki
nkw=n**2
iters = 501  # liczba iteracji

length_wiggle_param = 0.3
diameter_wiggle_param = 0.1

SPARSE = 1 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

F0=0.02
F1=0.1
z0=0
z1=0.1

center_node=0

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


def d_update(F):
    #zmiana średnicy pod względem siły F
    if (F > F0):
        if (F < F1):
            return z0 + (F - F0) * (z1 - z0) / (F1 - F0)
        else:
            return z1
    else:
        return z0

def update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges):
    reg_reg_edges2, reg_something_edges2, in_edges2=[], [], []
    for n1, n2, d, l in reg_reg_edges:
        F = 100 * c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew = d + d_update(F)
        reg_reg_edges2.append((n1, n2, dnew, l))


    for n1, n2, d, l in reg_something_edges:
        F = 100 * c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew = d + d_update(F)
        reg_something_edges2.append((n1, n2, dnew, l))
    for n1, n2, d, l in in_edges:
        F = 100 * c1 * c2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        dnew = d + d_update(F)
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

        data = sorted(data, key=lambda elem: elem[1])

        global indices, indptr
        return spr.csr_matrix(([datum[0] for datum in data], indices, indptr), shape=(nkw, nkw))


    if SPARSE == 1: return update_sparse_matrix(G)
    elif SPARSE == 0: return update_nonsparse_matrix(G)



#=================  GRAF I GEOMETRIA  ==================================================================================
#G, center_node=De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
G = Tr.Build_triangular_net(n, length_wiggle_param=length_wiggle_param, diameter_wiggle_param=diameter_wiggle_param)

#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=20, R_s=5)
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='rect', R=40, R_s=5, center_node=center_node)


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
color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',]
presult = create_pressure_flow_vector(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)
if SPARSE: indices, indptr = matrix.indices, matrix.indptr

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')
    matrix = update_matrix(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)
    
    if i%1 == 0:
        dhist=[[], [], [], [], [], [], []]
        qhist=[[], [], [], [], [], [], []]
        shearhist=[[], [], [], [], [], [], []]
        Fhist=[[], [], [], [], [], [], []]
        colors=[]
        dmax=1
        shearmax=0
        Fmax=0
        
        for n1, n2, d, l in reg_reg_edges:
            G[n1][n2]['d']=d
            q=c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q']=q

        for n1, n2, d, l in reg_something_edges:        
            G[n1][n2]['d']=d
            q=c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q']=q
        
        qmax = max([edge[2] for edge in G.edges(data='q')])
        for n1, n2 in G.edges():
            q=G[n1][n2]['q']
            d=G[n1][n2]['d']
            F=100*c2*q/d**3
            shear=d_update(F)
            if d>dmax:
                dmax=d
            if shear>shearmax:
                shearmax=shear
            if F>Fmax:
                Fmax=F
            qhist[int(6*q/qmax)].append(q)
            dhist[int(6*q/qmax)].append(d)
            Fhist[int(6*q/qmax)].append(F)
            shearhist[int(6*q/qmax)].append(shear)
            
        for edge in G.edges(data="q"):
            if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
                colors.append(color[int(6*edge[2]/qmax)])

        plt.figure(figsize=(20, 20))
        plt.suptitle('Flow for n = '+str(n)+", iters = " + str(i))
        spec = gridspec.GridSpec(ncols=4, nrows=2, height_ratios=[5, 1])
        plt.subplot(spec.new_subplotspec((0, 0), colspan=4))
        Dr.drawhist(G, n, f'{i//50:04d}.png', colors=colors, in_nodes=in_nodes, out_nodes=out_nodes)
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
        plt.subplot(spec[6]).set_title('Force')
        cindex=0
        plt.xlim((0,1.1*Fmax))
        for hist in Fhist:
            if len(hist)>1:
                plt.hist(hist, bins=50, color=color[cindex])
            cindex+=1
        plt.yscale("log")
        plt.subplot(spec[7]).set_title('Shear')
        cindex=0
        plt.xlim((0,1.1*shearmax))
        for hist in shearhist:
            if len(hist)>1:
                plt.hist(hist, bins=50, color=color[cindex])
            cindex+=1
        plt.yscale("log")
        plt.savefig(f'{i:04d}.png')
        plt.close()
        
    reg_reg_edges, reg_something_edges, in_edges=update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)

