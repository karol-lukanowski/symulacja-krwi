import numpy as np
from utils import fParams

class simInputData:
    n = 51 # rozmiar siatki
    iters = 1801  # liczba iteracji
    plot_every = 30
    save_every = 200

    noise = ["uniform", 1, 0.9] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
    #noise = ["gaussian", 2, 0.2] #gaussowski rozkład srednic, mu, sigma
    #noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

    arterioles_and_venules = True

    oxygen = False
    signal = False

    shear_d = True
    vegf_d = False
    signal_d = False
    gradp_d = False
    decrease_d = False

    data_collection = False

    qin = 10  # ilosć wpływającej krwi
    presout = 0  # cisnienie na wyjsciu
    mu = 0.0035  # współczynnik lepkosci
    l = 1  # początkowa długosć krawędzi
    c1 = np.pi / 128  # stała przepływu
    c2 = 64 / np.pi  # stała siły

    F_mult = 1000
    dt = 0.05

    kb = 0.00001

    D = 0.1 # współczynnik dyfuzji
    k = 0.001 # stała reakcji
    dth = 5 # graniczna grubosć
    dmin = 1
    dmax = 20
    #dth = 15
    Dv = 0.2 # współczynnik dyfuzji VEGF

    F_p = fParams({'F0': 0.2, 'F1': 2, 'z0': 0, 'z1': 0.1})
    F_gradp = fParams({'F0': 0.0001, 'F1': 0.001, 'z0': 0, 'z1': 0.01})
    F_ox = fParams({'F0': 0.1, 'F1': 1.5, 'z0': 0, 'z1': 0.4})
    F_s = fParams({'F0': 0.5, 'F1': 2, 'z0': 0, 'z1': 0.1})
    
    R_c = 0.3
    R_a = 15

    cp = 1

    F_mult_ox = 0.005
    dt_ox = 0.4

    ks = 0.05
    v = -1
    R = 1
    cs = 0.00005

    cox_in = 1
    cox_out = 0.7

    dec = 0.001

    pruning_iters = 0
    pruning_type = "flow"
    pruning_th = 0.01
    pruning_step = 0.005

    qdrawconst = 6
    ddrawconst = 3



    load = 0 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci, 2- wczytanie jednej z templatek
    #load_name = 'templates/ladder101' #load = 1 - nr folderu, load = 2 - nazwa templatki
    load_name = 'rect51/4'
    #load_name = 'own201/6'
    #load_name = 'deown201/2'
    #load_name = 'derect101derect101template/8'

    #geo = "cylindrical"
    #geo = "donut"
    geo = "rect"
    #geo = "top"
    #geo = "own"

    periodic = 'top' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    #in_nodes_own, out_nodes_own = np.array([[40, 40]]) / 100 * n, np.array([[60, 60]]) / 100 * n
    in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n

    length_wiggle_param = 1
    
    nsq = n ** 2
    old_iters = 0
    dirname = geo + str(n)