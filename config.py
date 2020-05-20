import numpy as np

import triangular_net as Tr
import delaunay as De
from geometry import set_geometry, equidistant_geometry

n = 351 # rozmiar siatki
nkw = n ** 2
iters = 301  # liczba iteracji
save_every = 30


length_wiggle_param = 1
diameter_wiggle_param = 3


qin = 50  # ilosć wpływającej krwi
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
dt = 0.8
#dt = 0.1

D = 1 # współczynnik dyfuzji
k = 0.1 # stała reakcji
dth = 10 # graniczna grubosć
#dth = 15
Dv = 1 # współczynnik dyfuzji VEGF

F0_ox = 0.1
F1_ox = 1
z0_ox = 0
z1_ox = 1
F_mult_ox = 0.01
dt_ox = 0.8


qdrawconst = 10
#qdrawconst = 5
ddrawconst = 10

G = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
#G = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)

#in_nodes, out_nodes = equidistant_geometry(G, n, R = n//2.5, xrange = n, yrange = n, how_many = 5)

#in_nodes = [nkw//2+n//8, nkw//2-1+n//8, nkw//2+1+n//8, nkw//2-n+n//8, nkw//2+n+n//8]
#out_nodes = [nkw//2+n//8 + n//4, nkw//4, nkw - nkw//4 - n//2, nkw//4 + n//6, nkw - nkw//4 - n//2 + n//6, nkw//4+n//2, nkw - nkw//4]


#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='cylindrical', R=n//2.5, **{'del': True})
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=0.75*n//2.5, R_s=n//20, **{'del': True})
in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='donut', R=(0.9*n)//2.5, R_s=n//20, **{'del': True})
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo = 'rect')
#in_nodes, out_nodes, reg_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)


in_nodes_ox = in_nodes
out_nodes_ox = out_nodes
