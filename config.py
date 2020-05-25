import numpy as np


n = 103 # rozmiar siatki
nkw = n ** 2
iters = 61  # liczba iteracji

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
F_mult = 500
dt = 0.8



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


load = 1
save_name = 'xd'
load_name = 'xd'

#geo = "cylindrical"
geo = "donut"
#geo = "rect"
#geo = "own"

#nettype = "tr"
nettype = "de"
