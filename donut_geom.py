#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin

mu = 0.0035  # współczynnik lepkosci
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

"""
W każdym kroku czasowym: 
1.  obliczamy cisnienia na podstawie wzoru na przepływy między węzłami
    oraz warunków brzegowych tj. stałe cisnienie na brzegach siatki + zerowanie się przepływów w węźle
2.  Na podstawie nowych cisnień obliczamy przepływy w każdej krawędzi
3.  Z przepływów obliczamy siłę działającą na scianki kanałów
4.  Następnie rozszerzamy scianki na podstawie wyliczonej siły
"""

R = 25
R_s = 5

def r_squared(G, node):
    n = int(len(G) ** 0.5)
    x0, y0 = G.nodes[n*n//2]["pos"]
    x, y = G.nodes[node]['pos']
    r_sqr = (x - x0) ** 2 + (y - y0) ** 2
    return r_sqr

#tworzenie tablicy ktora trzyma info czy jest node w wewnetrznym okregu czy nie
def get_inner_circle(G):
    inner_circle = np.zeros(len(G))
    for node in G.nodes():
        if r_squared(G, node) < R_s**2:
            inner_circle[node] = 1
    return inner_circle

def donut_setup(G):
    # lista krawedzi przekraczajacych boundary, pierwszy node to niech bedzie ten w srodku
    global inner_circle, boundary_nodes, boundary_edges
    inner_circle = get_inner_circle(G)
    boundary_edges = []
    for edge in G.edges():
        n1, n2 = edge
        if inner_circle[n1] + inner_circle[n2] == 1:
            # czyli jestesmy na krawedzi przekraczajacej obwod wewnetrznego okregu
            r1, r2 = r_squared(G, n1), r_squared(G, n2)
            if (r1 <= r2):
                boundary_edges.append((n1, n2))
            else:
                boundary_edges.append((n2, n1))
    # i jeszcze zbiór brzegowych nodes:
    boundary_nodes = set([edge[0] for edge in boundary_edges])
    boundary_nodes = sorted(list(boundary_nodes))

    return inner_circle, boundary_nodes, boundary_edges


def create_pressure_flow_vector(G, n, qin, presout):
    presult = np.zeros(n*n)
    presult[boundary_nodes] = -qin
    return presult

def create_matrix(G, SPARSE=0):
    n = int(len(G)**0.5)
    if (SPARSE == 0): return np.zeros((n * n, n * n))
    if (SPARSE == 1): return spr.csc_matrix(([], ([], [])), shape=(n * n, n * n))

def update_matrix(G, matrix, SPARSE=0):
    global boundary_edges, boundary_nodes
    n = int(len(G) ** 0.5)
    def update_nonsparse_cylindrical_matrix(G):

        for node in G.nodes:
            r_sqr = r_squared(G, node)
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

def donut_fix(inner_circle, boundary_nodes, pnow):
    p_inner = pnow[boundary_nodes[0]]
    for j in range(len(inner_circle)):
        if inner_circle[j] == 1:
            pnow[j] = p_inner
    return pnow