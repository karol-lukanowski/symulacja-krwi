#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def PosGauss(mean, sigma):
    x = np.random.normal(mean, sigma)
    return(x if x>=0 else PosGauss(mean,sigma))
    
