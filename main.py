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
diameter_wiggle_param = 1

SPARSE = 0 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

G = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)

#=================  FUNKCJE WSPÓLNE DLA KAŻDEJ GEOMETRII  ==============================================================

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

#=================  IMPORT GEOMETRII  ==================================================================================

from rectangular_geom import *

#from cylindrical_geom import *

#from donut_geom import *
#inner_circle, boundary_edges, boundary_nodes = donut_setup(G)

#=================  PROGRAM WŁAŚCIWY  ==================================================================================

presult = create_pressure_flow_vector(G, n, qin, presout)
matrix = create_matrix(G, SPARSE=SPARSE)

for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    matrix = update_matrix(G, matrix, SPARSE)
    pnow = solve_equation_for_pressure(matrix, presult)

    # TYLKO DLA DONUT GEOM - DLA INNEJ ZAKOMENTUJ
    #pnow = donut_fix(inner_circle, boundary_nodes, pnow)

    G = update_graph(G, pnow)

    if (i%50 == 0):
        Dr.drawq(G, n, f"q{i}.png")
