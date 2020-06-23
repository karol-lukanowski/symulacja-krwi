import numpy as np



n = 151 # rozmiar siatki
nkw = n ** 2
iters = 301  # liczba iteracji
save_every = 30

#noise = ["uniform", 1, 1] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
noise = ["gaussian", 1, 0.3] #gaussowski rozkład srednic, mu, sigma
#noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

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
dt = 0.9


D = 1 # współczynnik dyfuzji
k = 0.1 # stała reakcji
dth = 8 # graniczna grubosć
#dth = 15
Dv = 1 # współczynnik dyfuzji VEGF

F0_ox = 0.1
F1_ox = 1
z0_ox = 0
z1_ox = 1
F_mult_ox = 0.005
dt_ox = 3.6


qdrawconst = 5
ddrawconst = 3



load = 0 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci, 2- wczytanie jednej z templatek
templatki_names = ['deslabiak','deavr','dehard','trslabiak','travr','trhard']
load_name = templatki_names[0] #load = 1 - nr folderu, load = 2 - nazwa templatki

#geo = "cylindrical"
geo = "donut"
#geo = "rect"
#geo = "own"

nettype = "de"
#nettype = "tr"

length_wiggle_param = 1