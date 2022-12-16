import numpy as np

class simInputData:
    n = 100 # rozmiar siatki
    iters = 1000 # liczba iteracji
    tmax = 1000
    plot_every = 100
    #save_every = 100

    # INCLUDE IN SIMULATION
    include_mu_d = False # don't use until you check formula for mu
    include_adaptive_dt = True
    
    include_shear_growth = True
    include_vegf_growth = True
    include_signal_growth = False
    include_shrink = False

    # SHEAR
    shear_impact = 1
    shear_a = 10
    shear_b = 1

    # OXYGEN
    adv_const = 1
    rec_bv_const = 0.02
    rec_t_const = 0.1
    
    # VEGF
    vegf_const = 1
    vegf_a = 10
    vegf_b = 0.3
    vegf_force_a = 10
    vegf_force_b = 1


    noise = ["gaussian", 1, 0.1] # fix truncnorm in Delaunay

    qin = 1 # przepływ na wejściu
    pout = 0  # jednostki? ciśnienie na wyjściu
    l = 1  # początkowa długosć krawędzi


    d_min = 0.1
    d_max = 1000
    d_th = 5

    # TIMESTEP
    dt = 0.01
    dt_growth_rate = 0.05 # maksymalny procent średnicy o jaki może urosnąć krawędź
    dt_max = 1

    # DRAWING
    figsize = 10
    normalize = True
    qdrawconst = 1
    ddrawconst = 1

    load = 0 # 0 - dane z config, 1 - wczytanie danych z ewoluowanej sieci (plik save), 2 - wczytanie templatki (plik template)
    load_name = 'rect100/G1.00Daeff1.00/7'
    vtk_name = 'network_100x100.vtk'

    #geo = "cylindrical"
    #geo = "donut"
    geo = "rect"
    #geo = "top"
    #geo = "own"

    periodic = 'top' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    #in_nodes_own, out_nodes_own = np.array([[25, 50]]) / 100 * n, np.array([[75, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    in_nodes_own, out_nodes_own = np.array([[20, 50]]) / 100 * n, np.array([[80, 50], [70, 25], [70, 75]]) / 100 * n

    nsq = n ** 2
    old_iters = 0
    old_t = 0
    dirname = geo + str(n) + '/'