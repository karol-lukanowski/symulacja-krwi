import numpy as np
import scipy.sparse as spr

from config import simInputData
from utils import find_viscosity

def sigmoid(x, alpha, beta):
    """ Return sigmoidal function s(x) = 1 / (1 + exp(-alpha * (x - beta)))
    
    Parameters
    -------
    x : numpy array
        arguments vector
    alpa, beta : float
        sigmoid parameters
    
    Returns
    -------
    numpy array
        results vector
    """
    return 1 / (1 + np.exp(-alpha * (x - beta)))

def update_diameters(sid, flow, diams):
    growth = np.zeros_like(diams)
    if sid.include_shear_growth:
        growth += sid.shear_impact * find_shear_growth(sid, flow, diams)
    if sid.include_vegf_growth:
        growth += find_vegf_growth()
    if sid.include_signal_growth:
        growth += find_signal_growth()
    if sid.include_shrink:
        growth += find_shrink()
        
    if sid.include_adaptive_dt:
        dt = sid.dt_growth_rate / np.max(growth / diams)
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    
    diams += growth * dt
    diams = (diams > sid.dmin) * diams + (diams <= sid.dmin) * sid.dmin
    return diams, dt
        
def find_shear_growth(sid, flow, diams):
    """_summary_

    Args:
        sid (_type_): _description_
        flow (_type_): _description_
        diams (_type_): _description_

    Returns:
        _type_: _description_
    """
    viscs = find_viscosity(sid, diams)
    force = viscs * np.abs(flow) / diams
    return sigmoid(force, sid.shear_a, sid.shear_b)