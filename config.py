import numpy as np
from utils import fParams

class simInputData:
    n = 101 # rozmiar siatki
    iters = 31  # liczba iteracji
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

    F_mult = 1
    dt = 0.05

    kb = 10

    D = 1 # współczynnik dyfuzji
    k = 0.1 # stała reakcji
    dth = 4 # graniczna grubosć
    #dth = 15
    Dv = 0.2 # współczynnik dyfuzji VEGF

    F_p = fParams({'F0': 0.2, 'F1': 1, 'z0': 0, 'z1': 0.05})
    F_gradp = fParams({'F0': 0.15, 'F1': 0.2, 'z0': 0, 'z1': 5})
    F_ox = fParams({'F0': 0.1, 'F1': 1, 'z0': 0, 'z1': 0.4})

    cp = 1

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



    load = 0 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci, 2- wczytanie jednej z templatek
    #load_name = template #load = 1 - nr folderu, load = 2 - nazwa templatki
    #load_name = 'derect101template'
    load_name = 'rect101/0'
    #load_name = 'deown201/2'
    #load_name = 'derect101derect101template/8'

    #geo = "cylindrical"
    geo = "donut"
    #geo = "rect"
    #geo = "own"

    in_nodes_own, out_nodes_own = [n/2, n/2], [n/2 + n/10, n/2 + n/10]

    length_wiggle_param = 1
    
    nsq = n ** 2
    old_iters = 0
    dirname = geo + str(n)