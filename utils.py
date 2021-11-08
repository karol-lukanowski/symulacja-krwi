import numpy as np
import os
import scipy.sparse.linalg as sprlin


def PosGauss(mean, sigma):
    x = np.random.normal(mean, sigma)
    return(x if x>=0 else PosGauss(mean,sigma))
    
def mu_d(d):
    #return mu * (1 + 1 / (15 * d + 0.1))
    #return mu
    d = 100 * d + 0.2
    return ((6*np.exp(-0.085 * d)+2.2-2.44*np.exp(-0.06*d**0.645))*(d/(d-1.1)) ** 2 + 1) * (d/(d-1.1)) ** 2

def solve_equation(matrix, result):
    return sprlin.spsolve(matrix, result)

def d_update(F, t):
    if (F > t.F0):
        if (F < t.F1):
            result = t.z0+(F-t.F0)*(t.z1-t.z0)/(t.F1-t.F0)
        else:
            result = t.z1
    else:
        result = t.z0
    return result

def make_dir(sid):
        if not os.path.isdir(sid.dirname):
            os.makedirs(sid.dirname)

        i = 0
        dirname2 = sid.dirname
        while (sid.dirname == dirname2):
            if not os.path.isdir(sid.dirname + "/" + str(i)):
                sid.dirname = sid.dirname + "/" + str(i)
            else:
                i += 1
        os.makedirs(sid.dirname)

class fParams():
    def __init__(self, params):
        self.F0 = params['F0']
        self.F1 = params['F1']
        self.z0 = params['z0']
        self.z1 = params['z1']

