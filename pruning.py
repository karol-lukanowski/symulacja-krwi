import pressure as Pr 
from build import F0, F1, z0, z1, F_mult, dt, c1
from config import pruning_iters, pruning_type, pruning_th, pruning_step, noise
import draw_net as Dr
import numpy as np
import utils as Ut

import matplotlib.pyplot as plt

def pruning(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult):
    th = pruning_th
    maxtab = []
    avrtab = []
    for i in range(pruning_iters):
        if pruning_type == "shear":
            G1, reg_reg_edges, reg_something_edges, in_edges, pnow = final_pruning_shear(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, th)
        elif pruning_type == "flow":
            G1, reg_reg_edges, reg_something_edges, in_edges, pnow = final_pruning_flow(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, th)
        Fmax = 0
        Favr = 0
        n = 0     
        if pruning_type == "shear":
            G1, reg_reg_edges, reg_something_edges, in_edges, pnow = final_pruning_shear(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, th)
        elif pruning_type == "flow":
            G1, reg_reg_edges, reg_something_edges, in_edges, pnow = final_pruning_flow(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, th)
        Dr.drawblood(name=f'prun_q_blood_shear{th:0.4f}.png', oxresult=oxresult, data='q')
        #Dr.drawblood(name=f'prun_d_blood_shear{th:0.4f}.png', oxresult=oxresult, data='d') 
        for i,e in enumerate(reg_reg_edges + reg_something_edges + in_edges):
            n1, n2, d, l = e
            F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
            if F > Fmax:
                Fmax = F
            Favr += F
            n += 1
        print ('Fmax = ', Fmax)
        print ('Favr = ', Favr / n)
        maxtab.append(Fmax)
        avrtab.append(Favr / n)
#        Dr.drawblood(name=f'prun_d_blood{th:0.4f}.png', oxresult=oxresult, data='d')
#        Dr.drawblood(name=f'prun_q_blood{th:0.4f}.png', oxresult=oxresult, data='q')
        th += pruning_step
    plt.plot(maxtab)
    plt.show()
    plt.plot(avrtab)
    plt.show()
                


def final_pruning_shear(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, Fth_prun):
    Fmax = 0
    for i,e in enumerate(reg_reg_edges + reg_something_edges + in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F > Fmax:
            Fmax = F
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 1
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 1
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        F = F_mult / 2 * d * np.abs(pnow[n1] - pnow[n2]) / l
        if F < Fmax * Fth_prun:
            d = 1
        in_edges[i] = (n1, n2, d, l)
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    pnow = Pr.solve_equation(pmatrix, presult)
    G1 = Pr.update_network(G1, reg_reg_edges, reg_something_edges, pnow)
    return G1, reg_reg_edges, reg_something_edges, in_edges, pnow

"""
def final_pruning_flow(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, qth_prun):
    qmax = max([edge[2] for edge in G1.edges(data='q')])
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
    
        q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        if q < qmax * qth_prun:               
            d = np.random.rand() * 3 + 1
        reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
    
        q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        if q < qmax * qth_prun:
            d = np.random.rand() * 3 + 1
        reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
    
        q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        if q < qmax * qth_prun:
            d = np.random.rand() * 3 + 1
        in_edges[i] = (n1, n2, d, l)
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    pnow = Pr.solve_equation(pmatrix, presult)
    G1 = Pr.update_network(G1, reg_reg_edges, reg_something_edges, pnow)
    return G1, reg_reg_edges, reg_something_edges, in_edges, pnow
"""
def final_pruning_flow(G1, reg_reg_edges, reg_something_edges, in_edges, pnow, presult, oxresult, qth_prun):
    qmax = max([edge[2] for edge in G1.edges(data='q')])
    for i,e in enumerate(reg_reg_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:               
                d = 1
            reg_reg_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(reg_something_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:
                d = 1
            reg_something_edges[i] = (n1, n2, d, l)
    for i,e in enumerate(in_edges):
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            if q < qmax * qth_prun:
                d = 1
            in_edges[i] = (n1, n2, d, l)
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    pnow = Pr.solve_equation(pmatrix, presult)
    G1 = Pr.update_network(G1, reg_reg_edges, reg_something_edges, pnow)
    return G1, reg_reg_edges, reg_something_edges, in_edges, pnow
