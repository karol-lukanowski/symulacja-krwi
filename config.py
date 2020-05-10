import numpy as np

import triangular_net as Tr
import delaunay as De
from geometry import set_geometry


n = 101 # rozmiar siatki
nkw = n ** 2
iters = 501  # liczba iteracji
save_every = 100


length_wiggle_param = 0.5
diameter_wiggle_param = 0.3


qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły

F0 = 0.2
F1 = 1
z0 = 0
z1 = 1
F_mult = 10000
dt = 0.3


D = 1 # współczynnik dyfuzji
k = 0.01 # stała reakcji
dth = 8 # graniczna grubosć

F0_ox = 0.1
F1_ox = 1
z0_ox = 0
z1_ox = 1
F_mult_ox = 10
dt_ox = 0.05


qdrawconst = 5
ddrawconst = 3

#G = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
G = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)

#in_nodes, out_nodes = equidistant_geometry(R = n//2.5, xrange = n, yrange = n, how_many = 100)

#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='cylindrical', R=n//2.5)
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=n//2.5, R_s=n//20)
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo = 'rect')
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)


