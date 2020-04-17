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

R = 25

def default_nodes(n, G=[]):
    id_center = n*n // 2
    in_nodes = [id_center]
    x0, y0 = G.nodes[id_center]["pos"]
    out_nodes = []
    for node in G.nodes:
        pos = G.nodes[node]["pos"]
        r = np.sqrt((pos[0] - x0) ** 2 + (pos[1] - y0) ** 2)
        if r > R and r < R+1:
            out_nodes.append(node)
    return in_nodes, out_nodes