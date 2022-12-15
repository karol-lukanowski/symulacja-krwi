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

def update_diameters(sid, inc_matrix, flow, diams, blood_vessels, in_nodes, out_nodes):
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
    diams = (diams > sid.d_min) * diams + (diams <= sid.d_min) * sid.d_min
    inlet = np.zeros(sid.nsq)
    for node in in_nodes + out_nodes:
        inlet[node] = 1
    blood_vessels = ((diams > sid.d_th) * 1) * (np.abs(inc_matrix) @ (np.abs(inc_matrix.transpose()) @ blood_vessels + inlet) != 0) # MAYBE?
    return diams, blood_vessels, dt
        
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