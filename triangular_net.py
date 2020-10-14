#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import utils as Ut

"""
Funkcja budująca siatkę trójkątną o wymiarach n x n węzłów, jako graf
"""
   
def Build_triangular_net(n, l = 1, length_wiggle_param = 0, noise = ["uniform", 1, 1]):
    G = nx.Graph()    
    
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
        def in_hexagon(x, y):
            x, y = abs(x), abs(y)
            if (x > 1 / (2 * (3 ** 0.5))) and (y > -(3 ** 0.5) * x + 1):
                return 0
            return 1
        for node in G.nodes:
            # stare, wadliwe losowanie w kole
            """
            if (node >= n and node < n*(n-1) and (node%n != 0) and ((node+1)%n != 0)):
                r = np.random.ranf() * 0.5
                fi = np.random.ranf() * 2 * np.pi
                dx = length_wiggle_param * r * np.cos(fi)
                dy = length_wiggle_param * r * np.sin(fi)
                pos = G.nodes[node]['pos']
                new_x, new_y = pos[0] + dx, pos[1] + dy
                G.nodes[node]['pos'] = (new_x, new_y)
            """
            # nowe losowanie w szesciokacie
            if (node >= n and node < n * (n - 1) and (node % n != 0) and ((node + 1) % n != 0)):
                dx, dy = np.random.uniform(-1 / 3 ** 0.5, 1 / 3 ** 0.5), np.random.uniform(-0.5, 0.5)
                while not in_hexagon(dx, dy):
                    dx, dy = np.random.uniform(-1 / 3 ** 0.5, 1 / 3 ** 0.5), np.random.uniform(-0.5, 0.5)
                dx *= length_wiggle_param
                dy *= length_wiggle_param
                pos = G.nodes[node]['pos']
                G.nodes[node]['pos'] = (pos[0]+dx, pos[1]+dy)

    def add_edges(diameter_wiggle_param, l):
        """
        Deklaracja połączeń: węzeł początkowy, węzeł końcowy, grubosć (z wylosowanym szumem), przepływ
        """
        for node in G.nodes:
            for ind in G.nodes[node]["neigh"]:
                if (ind>node):
                    if noise[0] == "uniform":
                        d0 = np.random.rand() * noise[2] + noise[1]
                    elif noise[0] == "gaussian":
                        d0 = Ut.PosGauss(noise[1], noise[2])
                    elif noise[0] == "lognormal":
                        d0 = np.random.lognormal(noise[1], noise[2])
                G.add_edge(node,ind, d=d0,q=0, length=l)
    def find_edges_lengths():
        for edge in G.edges():
            node, neigh = edge
            pos1 = G.nodes[node]['pos']
            pos2 = G.nodes[neigh]['pos']
        
            len = np.linalg.norm(np.array(pos1)-np.array(pos2))
            if (len <= 1+length_wiggle_param): # czyli nie zmieniaj dlugosci miedzy brzegowymi wezlami, ona jest stale 1
                G[node][neigh]['length'] = len
    def find_boundary_edges():
        boundary_edges = [] 
        for edge in G.edges(data="q"):
            if ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
                boundary_edges.append((edge[0], edge[1]))        
        return boundary_edges
    
    add_nodes()
    wiggle_nodes(length_wiggle_param)
    add_edges(diameter_wiggle_param, l)
    find_edges_lengths()
    boundary_edges = find_boundary_edges()
    return G, boundary_edges
