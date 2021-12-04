import numpy as np
from utils import fParams

class simInputData:
    n = 101 # rozmiar siatki
    iters = 901  # liczba iteracji
    save_every = 30

    #noise = ["uniform", 1, 3] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
    noise = ["gaussian", 2, 0.2] #gaussowski rozkład srednic, mu, sigma
    #noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

    oxygen = False
    signal = True

    shear_d = True
    vegf_d = True
    signal_d = True
    gradp_d = False
    decrease_d = True

    data_collection = False

    qin = 10  # ilosć wpływającej krwi
    presin = 1 / 15
    presout = 0  # cisnienie na wyjsciu
    mu = 0.0035  # współczynnik lepkosci
    l = 1  # początkowa długosć krawędzi
    c1 = np.pi / 128  # stała przepływu
    c2 = 64 / np.pi  # stała siły

    F_mult = 1000
    dt = 0.05

    kb = 10

    D = 1 # współczynnik dyfuzji
    k = 0.1 # stała reakcji
    dth = 5 # graniczna grubosć
    dmin = 1
    dmax = 20
    #dth = 15
    Dv = 0.2 # współczynnik dyfuzji VEGF

    F_p = fParams({'F0': 0.2, 'F1': 5, 'z0': 0, 'z1': 0.02})
    F_gradp = fParams({'F0': 0.0001, 'F1': 0.001, 'z0': 0, 'z1': 0.05})
    F_ox = fParams({'F0': 0.1, 'F1': 1.5, 'z0': 0, 'z1': 0.4})
    F_s = fParams({'F0': 0.2, 'F1': 1, 'z0': 0, 'z1': 0.08})
    
    R_c = 0.3
    R_a = 15

    cp = 1

    F_mult_ox = 0.005
    dt_ox = 0.4

    ks = 0.1
    v = -1
    R = 1
    cs = 0.000015

    dec = 0.005

    pruning_iters = 0
    pruning_type = "flow"
    pruning_th = 0.01
    pruning_step = 0.005

    qdrawconst = 5
    ddrawconst = 3



    load = 0 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci, 2- wczytanie jednej z templatek
    #load_name = template #load = 1 - nr folderu, load = 2 - nazwa templatki
    #load_name = 'own101/1/template/0'
    load_name = 'own101/7'
    #load_name = 'deown201/2'
    #load_name = 'derect101derect101template/8'

    #geo = "cylindrical"
    #geo = "donut"
    #geo = "rect"
    geo = "own"

    periodic = 'all'

    in_nodes_own, out_nodes_own = [[45, 45]], [[55, 55]]
    #in_nodes_own, out_nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60]], [[50, 64], [50, 36], [36, 50], [64, 50]] #lista pozycji nodów in i out


    length_wiggle_param = 1
    
    nsq = n ** 2
    old_iters = 0
    dirname = geo + str(n)