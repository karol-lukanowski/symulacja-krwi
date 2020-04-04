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
iters = 350 # liczba iteracji

R= 40#n/2.5
R_s = 28

length_wiggle_param = 1
diameter_wiggle_param = 1

qdrawconst = 5
ddrawconst = 3

SPARSE = 0 # 1 = twórz macierz rzadką, 0 = twórz zwykłą macierz

qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły


#=================  ITERACJA  ==========================================================================================
"""
W każdym kroku czasowym: 
1.  obliczamy cisnienia na podstawie wzoru na przepływy między węzłami
    oraz warunków brzegowych tj. stałe cisnienie na brzegach siatki + zerowanie się przepływów w węźle
2.  Na podstawie nowych cisnień obliczamy przepływy w każdej krawędzi
3.  Z przepływów obliczamy siłę działającą na scianki kanałów
4.  Następnie rozszerzamy scianki na podstawie wyliczonej siły
"""

def r_squared(node):
    x0, y0 = G.nodes[n*n//2]["pos"]
    x, y = G.nodes[node]['pos']
    r_sqr = (x - x0) ** 2 + (y - y0) ** 2
    return r_sqr
def create_pressure_flow_vector(qin, presout):
    """
    Wektor przepływów/cisnień, będący wynikiem działania macierzy na wektor cisnień
    pierwsze n elementów jest takie same, co odpowiada stałemu wpływowi do sieci
    dalej mamy elementy równe 0, co odpowiada sumie wpływów do każdego węzła równej 0
    ostatnie n elementów jest 0, bo ostatnie równania odpowiadają trzymaniu cisnienia
    wyjsciowego stale równego 0 (inaczej w układzie szybko pojawiają się ujemne cisnienia)
    """
    id_center = int(n*n/2)
    x0,y0 = G.nodes[id_center]["pos"]

    presult = np.zeros(n * n)
    for i in range(n * n):
        pos = G.nodes[i]["pos"]
        r = np.sqrt((pos[0] - x0)**2 + (pos[1] - y0)**2)
            
        if r <= R_s:
            presult[i] = -qin
        elif r >= R:
            presult[i] = presout
    return presult

def create_matrix(G):
    if (SPARSE == 0): return np.zeros((n * n, n * n))
    if (SPARSE == 1): return spr.csc_matrix(([], ([], [])), shape=(n * n, n * n))
def update_matrix(G):
    def update_nonsparse_cylindrical_matrix(G):

        for node in G.nodes:
            r_sqr = r_squared(node)
            if r_sqr < R**2 and r_sqr > R_s**2:
                matrix[node][node] = 0
                for ind in G.nodes[node]["neigh"]:
                    d = G[node][ind]["d"]
                    l = G[node][ind]['length']
                    matrix[node][ind] = c1 * d ** 4 / l
                    matrix[node][node] -= c1 * d ** 4 / l
            else:
                matrix[node][node] = 1

        for node in boundary_nodes:
            for i in range(n*n):
                matrix[node][i] = 0

        for edge in boundary_edges:
            n1, n2 = edge
            d = G[n1][n2]['d']
            l = G[n1][n2]['length']
            for node in boundary_nodes:
                matrix[node][node] -= c1 * d**4 / l
                matrix[node][n2] += c1 * d**4 / l

        return matrix
    def update_nonsparse_matrix(G):
        """
        macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
        z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
        srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
        """
        #matrix = np.zeros((n * n, n * n))
        for node in G.nodess:
            if (node >= n and node < n * (n - 1)):
                matrix[node][node] = 0
                for ind in G.nodess[node]["neigh"]:
                    d = G[node][ind]["d"]
                    l = G[node][ind]['length']
                    matrix[node][ind] = c1 * d ** 4 / l
                    matrix[node][node] -= c1 * d ** 4 / l
            elif (node < n):
                matrix[node][node] = 0
                for ind in range(n):
                    d = G[ind][ind + n]["d"]
                    l = G[ind][ind + n]['length']
                    matrix[node][ind + n] = c1 * d ** 4 / l
                    matrix[node][node] -= c1 * d ** 4 / l
            else:
                matrix[node][node] = 1
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

        for node in G.nodess:
            if (node >= n and node < n * (n - 1)):
                row = np.append(row, node)
                col = np.append(col, node)
                data = np.append(data, 0)
                k = data.size - 1
                for ind in G.nodess[node]["neigh"]:
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
    elif SPARSE == 0: return update_nonsparse_cylindrical_matrix(G)

def solve_equation_for_pressure(matrix, presult):
    """
    Zamieniamy macierz w równaniu na formę macierzy rzadkiej
    w celu usprawnienia obliczeń
    """
    if (SPARSE == 0): pnow = sprlin.spsolve(spr.csc_matrix(matrix), presult)
    elif (SPARSE == 1): pnow = sprlin.spsolve(matrix, presult)
    return pnow
def update_pressure_in_nodes(G, pnow):
    for node in G.nodes:
        G.nodes[node]["p"] = pnow[node]
    return G

def find_edge_parameters(G, node1, node2):
    p1 = G.nodes[node1]["p"]
    p2 = G.nodes[node2]["p"]
    d = G[node1][node2]["d"]
    l = G[node1][node2]['length']
    return (p1, p2, d, l)
def find_flow(p1, p2, l, d):
    """
    Przepływ pomiędzy dwoma węzłami
    """
    dp = p1 - p2
    q = c1 * d ** 4 * dp / l
    return q
def find_force(q, d):
    """
    Siła z jaką ciecz działa na scianki kanału
    """
    f = c2 * np.abs(q) / d ** 3
    return f

def update_flow_in_edge(G, node1, node2, q):
    G[node1][node2]['q'] = q
    return G
def update_diameter_in_edge(G, node1, node2, F):
    """
    Sposób, w jaki siła działa na srednice jest największym problemem, wynik bardzo mocno od tego zależy
    """
    if (F > 0.00002):
        if (F < 0.0002):
            G[node1][node2]["d"] += 1000 * F
        else:
            G[node1][node2]["d"] += 0.2
    return G


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     PROGRAM WŁASCIWY     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

G = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)

#tablica ktora trzyma info czy jest node w wewnetrznym okregu czy nie
inner_circle = np.zeros(len(G))
for node in G.nodes():
    if r_squared(node) < R_s**2:
        inner_circle[node] = 1

#lista krawedzi przekraczajacych boundary, pierwszy node to niech bedzie ten w srodku
boundary_edges = []
for edge in G.edges():
    n1, n2 = edge
    if inner_circle[n1] + inner_circle[n2] == 1:
        #czyli jestesmy na krawedzi przekraczajacej obwod wewnetrznego okregu
        r1, r2 = r_squared(n1), r_squared(n2)
        if (r1 <= r2):
            boundary_edges.append((n1, n2))
        else:
            boundary_edges.append((n2, n1))

#i jeszcze zbiór brzegowych nodes:
boundary_nodes = set([edge[0] for edge in boundary_edges])
boundary_nodes=sorted(list(boundary_nodes))


#presult = create_pressure_flow_vector(qin,presout)

# nowy sposob na stworzenie presulta
presult = np.zeros(n*n)
presult[boundary_nodes] = -qin

matrix = create_matrix(G)
for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    matrix = update_matrix(G)
    pnow = solve_equation_for_pressure(matrix, presult)

    p_inner = pnow[boundary_nodes[0]]
    for j in range(len(inner_circle)):
        if inner_circle[j] == 1:
            pnow[j] = p_inner

    G = update_pressure_in_nodes(G, pnow)

    for node in G.nodes:
        for ind in G.nodes[node]["neigh"]:
            if (ind > node):
                p1, p2, d, l = find_edge_parameters(G, node, ind)
                q = find_flow(p1, p2, l, d)
                G = update_flow_in_edge(G, node, ind, q)
                F = find_force(q, d)
                G = update_diameter_in_edge(G, node, ind, F)

    """
    Q = 0
    for edge in boundary_edges:
        n1, n2 = edge
        Q += abs(G[n1][n2]['q'])

    print(Q)
    """


    if i % 50 == 0:
        Dr.drawq(G, n, f"q{i}.png", qdrawconst)
        # Dr.drawd(G, n, f"d{i}.png", ddrawconst)

