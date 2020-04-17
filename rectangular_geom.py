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

def default_nodes(n, G=[]):
    in_nodes = list(range(n))
    out_nodes = list(range(n*(n-1), n*n))
    return in_nodes, out_nodes