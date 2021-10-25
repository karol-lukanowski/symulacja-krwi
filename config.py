import numpy as np



n = 101 # rozmiar siatki
nkw = n ** 2
iters = 601  # liczba iteracji
save_every = 30

#noise = ["uniform", 1, 3] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
noise = ["gaussian", 1, 0.3] #gaussowski rozkład srednic, mu, sigma
#noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

qin = 10  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
mu = 0.0035  # współczynnik lepkosci
l = 1  # początkowa długosć krawędzi
c1 = np.pi / 128  # stała przepływu
c2 = 64 / np.pi  # stała siły

F0 = 0.2
F1 = 1
z0 = 0

z1 = 1
F_mult = 1
dt = 0.05

kb = 10

D = 1 # współczynnik dyfuzji
k = 0.1 # stała reakcji
dth = 4 # graniczna grubosć
#dth = 15
Dv = 0.2 # współczynnik dyfuzji VEGF


F0_p = 0.15
F1_p = 0.2
z0_p = 0
z1_p = 5
cp = 1

F0_ox = 0.1
F1_ox = 1
z0_ox = 0

z1_ox = 1
F_mult_ox = 0.005
dt_ox = 0.4

ks = 0.0001
v = -1
R = 1
cs = 0.00015

pruning_iters = 0
pruning_type = "flow"
pruning_th = 0.01
pruning_step = 0.005

qdrawconst = 5
ddrawconst = 3



load = 2 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci, 2- wczytanie jednej z templatek
templatki_names = ['deslabiak','deavr','dehard','trslabiak','travr','trhard']
#load_name = templatki_names[0] #load = 1 - nr folderu, load = 2 - nazwa templatki
#load_name = 'derect101template'
#load_name = 'deown101own101/0'
load_name = 'own101'
#load_name = 'derect101derect101template/8'

#geo = "cylindrical"
#geo = "donut"
geo = "rect"
#geo = "own"

in_nodes_own, out_nodes_own = [n/2, n/2], [n/2 + n/10, n/2 + n/10]

nettype = "de"
#nettype = "tr"

length_wiggle_param = 1
