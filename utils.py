#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from config import mu


def PosGauss(mean, sigma):
    x = np.random.normal(mean, sigma)
    return(x if x>=0 else PosGauss(mean,sigma))
    
def mu_d(d):
    #return mu * (1 + 1 / (15 * d + 0.1))
    #return mu
    d = 100 * d+0.2
    return ((6*np.exp(-0.085 * d)+2.2-2.44*np.exp(-0.06*d**0.645))*(d/(d-1.1)) ** 2 + 1) * (d/(d-1.1)) ** 2
