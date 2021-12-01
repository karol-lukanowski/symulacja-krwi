import numpy as np

from config import simInputData

def update_graph(sid:simInputData, edges):
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e        
        d -= sid.dec
        if d < sid.dmin:
            d = sid.dmin
        elif d > sid.dmax:
             d = sid.dmax
        edges[i] = (n1, n2, d, l, t)

    return edges