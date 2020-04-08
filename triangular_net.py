#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

"""
Funkcja budująca siatkę trójkątną o wymiarach n x n węzłów, jako graf
"""
   
def Build_triangular_net(n, l = 1, length_wiggle_param = 0, diameter_wiggle_param = 0):
    G = nx.Graph()    
    lens = []
    # sąsiedztwo: pierwsza i ostatnia kolumna siatki ma tylko po jednym połączeniu,
    # żeby łatwiej trzymać stały wpływ i cisnienie; poza tym standardowe sąsiedztwo
    # dla trójkatnej siatki
    def neighbours(index):
        # czysto techniczne funkcje do warunków brzegowych
        def boundaryr(i):
            if (i < n - 1):
                return i + 1
            else:
                return 0
        def boundaryl(i):
            if (i > 0):
                return i - 1
            else:
                return n - 1

        j = index % n
        i = int((index - j) / n)
        neightab = []
        if (i > 1) and (i < n - 2):
            neightab.append(n * i + boundaryr(j))
            neightab.append(n * i + boundaryl(j))
            if (i % 2 == 0):
                if (i > 0):
                    neightab.append(n * (i - 1) + j)
                    neightab.append(n * (i - 1) + boundaryl(j))
                if (i < n - 1):
                    neightab.append(n * (i + 1) + j)
                    neightab.append(n * (i + 1) + boundaryl(j))
            else:
                if (i > 0):
                    neightab.append(n * (i - 1) + j)
                    neightab.append(n * (i - 1) + boundaryr(j))
                if (i < n - 1):
                    neightab.append(n * (i + 1) + j)
                    neightab.append(n * (i + 1) + boundaryr(j))
        elif (i == 0):
            neightab.append(n * (i + 1) + j)
        elif (i == 1):
            neightab.append(n * (i - 1) + j)
            neightab.append(n * i + boundaryr(j))
            neightab.append(n * i + boundaryl(j))
            neightab.append(n * (i + 1) + j)
            neightab.append(n * (i + 1) + boundaryr(j))
        elif (i == n - 2):
            neightab.append(n * (i + 1) + j)
            neightab.append(n * i + boundaryr(j))
            neightab.append(n * i + boundaryl(j))
            if (i % 2 == 0):
                neightab.append(n * (i - 1) + j)
                neightab.append(n * (i - 1) + boundaryl(j))
            else:
                neightab.append(n * (i - 1) + j)
                neightab.append(n * (i - 1) + boundaryr(j))
        elif (i == n - 1):
            neightab.append(n * (i - 1) + j)
        return neightab

    def add_nodes():   
        h = 3**0.5 / 2           # wysokosc trojkata równobocznego
        for i in range(n):
            for j in range(n):
                index=n*i+j
                if (i%2==0):
                    pos=(i*h,j)
                else:
                    pos=(i*h,j+0.5)
                if (i==0):
                    G.add_node(index,p=1,pos=pos,neigh=neighbours(index))
                else:
                    G.add_node(index,p=0,pos=pos,neigh=neighbours(index))
    def wiggle_nodes(length_wiggle_param):
        """
        Przesuwanie wezlow zgodnie z length_wiggle_param, oprócz węzłów brzegowych
        """
        for node in G.nodes:
            if (node >= n and node < n*(n-1) and (node%n != 0) and ((node+1)%n != 0)):
                r = np.random.ranf() * 0.5
                fi = np.random.ranf() * 2 * np.pi
                dx = length_wiggle_param * r * np.cos(fi)
                dy = length_wiggle_param * r * np.sin(fi)
                pos = G.nodes[node]['pos']
                new_x, new_y = pos[0] + dx, pos[1] + dy
                G.nodes[node]['pos'] = (new_x, new_y)
        import matplotlib.pyplot as plt

    def add_edges(diameter_wiggle_param, l):
        """
        Deklaracja połączeń: węzeł początkowy, węzeł końcowy, grubosć (z wylosowanym szumem), przepływ
        """
        for node in G.nodes:
            for ind in G.nodes[node]["neigh"]:
                if (ind>node):
                    d0 = np.random.rand() * diameter_wiggle_param + 1
                    G.add_edge(node,ind, d=d0,q=0, length=l)
    def find_edges_lengths():
        for edge in G.edges():
            node, neigh = edge
            pos1 = G.nodes[node]['pos']
            pos2 = G.nodes[neigh]['pos']
        
            len = np.linalg.norm(np.array(pos1)-np.array(pos2))
            if (len <= 1+length_wiggle_param): # czyli nie zmieniaj dlugosci miedzy brzegowymi wezlami, ona jest stale 1
                G[node][neigh]['length'] = len
         
    add_nodes()
    wiggle_nodes(length_wiggle_param)
    add_edges(diameter_wiggle_param, l)
    find_edges_lengths()
    return G
