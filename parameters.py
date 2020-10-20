import numpy as np
import pressure as Pr
import oxygen as Ox
import vegf as Ve
import save as Sv
import draw_net as Dr
import matplotlib.pyplot as plt
import networkx as nx

import os
import re

i = 0
params = [x[1] for x in os.walk('./spis_powszechny')]
vegftab = []
oxtab = []
hittab = []
tab = []
for directory in params[0]:
    word = list(directory)
    hits = []
    hitangle = []
    hitfinal = 0
    if (word[0] == 'D' and word[1] == 'v'):
        Dv = float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", directory)[0])
        dirname = 'spis_powszechny/' + directory
        (G, boundary_edges, n, F0, F1, z0, z1, F_mult, dt, c1, c2, l, mu, qin, presout, D, Dv, k, dth,
         F0_ox, F1_ox, z0_ox, z1_ox, F_mult_ox, dt_ox, old_iters,
         in_nodes, out_nodes, reg_nodes, other_nodes, boundary_nodes_out, boundary_nodes_in,
         in_nodes_ox, out_nodes_ox, oxresult, in_edges, reg_reg_edges, reg_something_edges, other_edges) = Sv.load(dirname+"/save")
        nkw = n ** 2
        presult = Pr.create_vector()
        pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
        oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges, k, D) 
        pnow = Pr.solve_equation(pmatrix, presult)
        oxnow = Ox.solve_equation(oxmatrix, oxresult)
        vresult = Ve.create_vector(oxnow, oxresult)
        vmatrix = Ve.update_matrix(vresult, reg_reg_edges, reg_something_edges, other_edges, Dv)
        vnow = Ve.solve_equation(vmatrix, vresult)
        for node in other_nodes:
            vnow[node] = 0
        tab.append(Dv)
        oxtab.append(np.average(oxnow))
        vegftab.append(np.average(vnow))
#        Dr.drawq(name=f'veq-check-{i:04d}.png', oxdraw=vnow+oxresult, dirname = 'spis_powszechny')
        for node in out_nodes:
            if oxnow[node] == 1:
                hits.append(node)  
        
        pos = nx.get_node_attributes(G, 'pos')
        id_center = int(n * n / 2)
        for node in hits:
            hitangle.append(np.arctan2(pos[node][0]-pos[id_center][0], pos[node][1]-pos[id_center][1]) * 180 / np.pi)
        for j in range(len(hitangle)):
            if (np.abs(hitangle[j]-np.roll(hitangle, 1)[j])<10 and np.abs(hitangle[j]-np.roll(hitangle, -1)[j])<10):
                hitfinal += 1
        hittab.append(hitfinal)
        i += 1

plt.xlabel('Dv')
plt.xscale('log')
plt.ylabel('average VEGF')
plt.plot(tab, vegftab, '.')
plt.show()

plt.xlabel('Dv')
plt.xscale('log')
plt.ylabel('average oxygen')
plt.plot(tab, oxtab, '.')
plt.show()

plt.xlabel('Dv')
plt.xscale('log')
plt.ylabel('veins')
plt.plot(tab, hittab, '.')
plt.show()