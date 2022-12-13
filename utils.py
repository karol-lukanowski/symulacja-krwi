import numpy as np
import os
import scipy.sparse.linalg as sprlin
import scipy.sparse as spr

def solve_equation(A, b):
    """ Solve matrix equation A * x = b.

    Parameters
    -------
    A : scipy sparse matrix
        matrix A from equation

    b : scipy sparse vector
        result b from equation

    Returns
    -------
    scipy sparse vector
        result x from equation    
    """
    return sprlin.spsolve(A, b)

def initialize_iterators(sid):
    """ Create iterators for simulation steps, time and other conditions.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        iters - max no. of iterations of new simulation
        tmax - max time of new simulation
        old_iters - no. of previous iterations (if loaded from saved file)
        old_t - time of previous simulation (if loaded from saved file)
    
    Returns
    -------
    iters : int 
        max no. of new iterations

    i : int
        iterator in range from old iterations to sum of old and new

    tmax : float
        max new time

    t : float
        time iterator in range from old time to sum of old and new

    breakthrough : bool
        parameter stating if the system was dissolved (if diameter of edge
        going to the output grew at least sid.d_break times)
    """
    iters = sid.old_iters + sid.iters
    tmax = sid.old_t + sid.tmax
    i = sid.old_iters
    t = sid.old_t
    breakthrough = False
    return iters, tmax, i, t, breakthrough

def update_iterators(sid, i, t, dt):
    """ Update iterators in simulation and in configuration class.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        old_iters - no. of previous iterations (if loaded from saved file)
        old_t - time of previous simulation (if loaded from saved file)

    i : int
        current iteration

    t : float
        current time

    dt : float
        current timestep

    Returns
    -------
    i : int
        current iteration

    t : float
        current time
    """
    i += 1
    sid.old_iters += 1 # update simulation iterations in configuration class 
    t += dt
    sid.old_t += dt # update simulation time in configuration class
    return i, t

def find_viscosity(sid, diams):
    """ Find viscosity of each edge, depending on diameter.
    
    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        mu0 - characteristic viscosity of an edge
        mu_d - use mu(d) dependence or not
    
    diams : numpy array
        current diameters of edges
        
    Returns
    -------
    numpy array
        viscosities of edges
    """
    if sid.include_mu_d:
        diams = 100 * diams + 0.2
        return ((6*np.exp(-0.085 * diams)+2.2-2.44*np.exp(-0.06*diams**0.645))*(diams/(diams-1.1)) ** 2 + 1) * (diams/(diams-1.1)) ** 2 / 1000 # check if makes sense
    else:
        return np.ones_like(diams)

def make_dir(sid):
    ''' Create directory for current simulation with geometry name, size and simulation index.
    
    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        dirname - directory name consisting of geometry and size
        
    Returns
    -------
    None   
    '''
    i = 0
    dirname2 = sid.dirname
    while (sid.dirname == dirname2):
        if not os.path.isdir(sid.dirname + "/" + str(i)):
            sid.dirname = sid.dirname + "/" + str(i)
        else:
            i += 1
    os.makedirs(sid.dirname)