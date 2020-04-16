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

def create_pressure_flow_vector(G, n, qin=1, presout=0):
    """
    Wektor przepływów/cisnień, będący wynikiem działania macierzy na wektor cisnień
    pierwsze n elementów jest takie same, co odpowiada stałemu wpływowi do sieci
    dalej mamy elementy równe 0, co odpowiada sumie wpływów do każdego węzła równej 0
    ostatnie n elementów jest 0, bo ostatnie równania odpowiadają trzymaniu cisnienia
    wyjsciowego stale równego 0 (inaczej w układzie szybko pojawiają się ujemne cisnienia)
    """
    id_center = int(n * n / 2)
    x0, y0 = G.nodes[id_center]["pos"]

    presult = np.zeros(n * n)
    for i in range(n * n):
        pos = G.nodes[i]["pos"]
        r = np.sqrt((pos[0] - x0) ** 2 + (pos[1] - y0) ** 2)

        if (i == id_center):
            presult[i] = -qin
        elif r >= R:
            presult[i] = presout
    return presult

def create_matrix(G, SPARSE=0):
    n = int(len(G)**0.5)
    if (SPARSE == 0): return np.zeros((n * n, n * n))
    if (SPARSE == 1): return spr.csc_matrix(([], ([], [])), shape=(n * n, n * n))

def update_matrix(G, matrix, SPARSE=0):
    n = int(len(G) ** 0.5)
    def update_nonsparse_cylindrical_matrix(G):
        """
        macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
        z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
        srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
        """
        id_center = int(n*n/2)
        x0,y0 = G.nodes[id_center]["pos"]
        #matrix = np.zeros((n * n, n * n))
        for node in G.nodes:
            pos = G.nodes[node]["pos"]
            r = np.sqrt((pos[0] - x0)**2 + (pos[1] - y0)**2)
            if r < R:
                matrix[node][node] = 0
                for ind in G.nodes[node]["neigh"]:
                    d = G[node][ind]["d"]
                    l = G[node][ind]['length']
                    matrix[node][ind] = c1 * d ** 4 / l
                    matrix[node][node] -= c1 * d ** 4 / l

            else:
                matrix[node][node] = 1
        return matrix
    def update_nonsparse_matrix(G):
        """
        macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
        z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
        srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
        """
        #matrix = np.zeros((n * n, n * n))
        for node in G.nodes:
            if (node >= n and node < n * (n - 1)):
                matrix[node][node] = 0
                for ind in G.nodes[node]["neigh"]:
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
    elif SPARSE == 0: return update_nonsparse_cylindrical_matrix(G)
