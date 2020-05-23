import numpy as np
import os
import triangular_net as Tr
import delaunay as De
from geometry import set_geometry, equidistant_geometry

Load = False
save_name = 'xd'
load_name = 'xd'

n = 103 # rozmiar siatki
nkw = n ** 2
iters = 91  # liczba iteracji
save_every = 30


length_wiggle_param = 0.1
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
F_mult = 100
dt = 1


D = 1 # współczynnik dyfuzji
k = 0.1 # stała reakcji
dth = 8 # graniczna grubosć
#dth = 15
Dv = 1 # współczynnik dyfuzji VEGF

F0_ox = 0.1
F1_ox = 1
z0_ox = 0
z1_ox = 1
F_mult_ox = 0.001
dt_ox = 3


qdrawconst = 5
ddrawconst = 3




#G, boundary_edges, nettype = De.Build_delaunay_net(n, diameter_wiggle_param=diameter_wiggle_param)
G, boundary_edges, nettype = Tr.Build_triangular_net(n, length_wiggle_param = length_wiggle_param, diameter_wiggle_param = diameter_wiggle_param)



#geo = "cylindrical"
geo = "donut"
#geo = "rect"
#geo = "own"

in_nodes, out_nodes, reg_nodes, other_nodes, in_edges = set_geometry(n, G, geo=geo, R=n//2.5, R_s=n//20, **{'del': False})



#in_nodes, out_nodes = equidistant_geometry(G, n, R = n//2.5, xrange = n, yrange = n, how_many = 200)
#in_nodes, out_nodes, reg_nodes, other_nodes, in_edges = set_geometry(n, G, geo='own', in_nodes=in_nodes, out_nodes=out_nodes)



in_nodes_ox = in_nodes
out_nodes_ox = out_nodes

dirname = nettype + geo + "n" + str(n) + "lw" + str(length_wiggle_param) + "dw" + str(diameter_wiggle_param) + "dt" + str(dt) + "dtox" + str(dt_ox)
if not os.path.isdir(dirname):
    os.makedirs(dirname)