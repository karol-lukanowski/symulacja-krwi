import numpy as np
from utils import fParams

class simInputData:
    """Class containing the parameters of a given simulation, transfered between files."""

    ####################################################################################################################
    # Loading a previously evolved network / loading a template / starting a new simulation with the config data below
    ####################################################################################################################

    loadMode = 'load_template'                         # 'load_network' / 'load_template' / 'new_simulation'
    loadDirectory = 'templates/ladder101'               # directory of existing network/template if loading

    ####################################################################################################################
    #######################################  Data for starting a new simulation  #######################################
    ####################################################################################################################
    # Simulation size settings
    ####################################################################################################################

    networkSize = 101                                     # network size in one dimension (networkSize**2 is the number of nodes)
    networkSizeSquared = networkSize ** 2
    iterations = 35                                    # number of iterations to be performed by the sim
    plot_every = 30                                     # make plots every _ iterations
    save_every = 200                                    # save the network every _ iterations

    ####################################################################################################################
    # Initial statistics of capillary diameters
    ####################################################################################################################

    diameterStats = ["gaussian", 2, 0.2]                # [gaussian distribution, average diameter, standard deviation]

    # other possible choices:
    #diameterStats = ["uniform", 1, 3]                  # [uniform distribution, srednica poczatkowa???, diamtere_wiggle_param???]
    #diameterStats = ["lognormal", 1, 0.3]              # [log-normal distribution, mu, sigma]

    ####################################################################################################################
    # Evolution mechanisms to include
    ####################################################################################################################
    includeOxygenEvolution = True
    includeSignalEvolution = True

    shear_d = True
    vegf_d = False
    signal_d = True
    gradp_d = False
    decrease_d = False

    data_collection = True

    ####################################################################################################################
    # Evolution parameters
    ####################################################################################################################

    ### general
    flow = 10                                           # constant flow of blood (amount flowing in, amount flowing out)
    outputPressure = 0                                  # constant value of pressure in output nodes
    viscosityCoeff = 0.0035                             # viscosity coefficient

    ### constants appearing in equations
    flowConstC1 = np.pi / 128  # stała przepływu
    forceConstC2 = 64 / np.pi  # stała siły
    forceMultiplier = 1000

    ### timestep in a single iteration
    dt = 0.05

    ### ???
    kb = 0.01

    ### oxygen equation parameters
    oxygenDiffusionConst = 1                            # oxygen diffusion constant
    oxygenReactionConst = 0.1                           # oxygen reaction constant
    thresholdDiameter = 5                               # threshold diameter; a capillary exceeding it becomes an arteriole/venule
    minimalDiameter = 1                                 # minimal possible diameter of a vessel
    maximalDiameter = 20                                # maximal possible diameter of a vessel

    ### VEGF equation parameters
    VEGFDiffusionConst = 0.2                            # VEGF diffusion constant

    F_p = fParams({'F0': 0.2, 'F1': 2, 'z0': 0, 'z1': 0.01})
    F_gradp = fParams({'F0': 0.0001, 'F1': 0.001, 'z0': 0, 'z1': 0.01})
    F_ox = fParams({'F0': 0.1, 'F1': 1.5, 'z0': 0, 'z1': 0.4})
    F_s = fParams({'F0': 0.2, 'F1': 1, 'z0': 0, 'z1': 0.1})
    
    R_c = 0.3
    R_a = 15

    cp = 1

    F_mult_ox = 0.005
    dt_ox = 0.4

    ks = 0.05
    v = -1
    R = 1
    cs = 0.000015

    dec = 0.002

    pruning_iters = 0
    pruning_type = "flow"
    pruning_th = 0.01
    pruning_step = 0.005

    qdrawconst = 5
    ddrawconst = 3



    #geo = "cylindrical"
    #geo = "donut"
    geo = "rect"
    #geo = "own"

    periodic = 'none'

    #in_nodes_own, out_nodes_own = [[45, 45]], [[55, 55]]
    in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * networkSize, \
                                  np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * networkSize #lista pozycji nodów in i out


    length_wiggle_param = 1
    
    old_iters = 0
    dirname = geo + str(networkSize)


