import numpy as np

import draw_net as Dr
import pressure as Pr
import oxygen as Ox
from config import n, G, in_nodes, out_nodes, iters, in_edges, c1, save_every
import matplotlib.pyplot as plt

reg_reg_edges, reg_something_edges = [],  []
for n1, n2 in G.edges():
    d = G[n1][n2]['d']
    l = G[n1][n2]['length']
    if (n1 not in in_nodes and n1 not in out_nodes) and (n2 not in in_nodes and n2 not in out_nodes):
        reg_reg_edges.append((n1, n2, d, l))
    elif (n1 not in in_nodes and n1 not in out_nodes):
        reg_something_edges.append((n1, n2, d, l))
    elif (n2 not in in_nodes and n2 not in out_nodes):
        reg_something_edges.append((n2, n1, d, l ))

presult = Pr.create_vector(G)
oxresult = Ox.create_vector(G)



for i in range(iters):
    print(f'Iter {i + 1}/{iters}')
    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)
    oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges)
    pnow = Pr.solve_equation(pmatrix, presult)
    oxnow = Ox.solve_equation(oxmatrix, oxresult)        
    
    if i%save_every == 0:
        Q_in = 0
        Q_out = 0
        
        for n1, n2, d, l in reg_reg_edges:
            G[n1][n2]['d']= d
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q'] = q

        for n1, n2, d, l in reg_something_edges:        
            G[n1][n2]['d'] = d
            q = c1 * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
            G[n1][n2]['q'] = q
            
            if n2 in in_nodes:
                Q_in += q
            if n2 in out_nodes:
                Q_out += q
        
        print('Q_in =', Q_in, 'Q_out =', Q_out)
        
        
        Dr.drawhist(name = f'{i//save_every:04d}.png', oxnow = oxnow, oxresult = oxresult)
#        Dr.drawd(name = f'd{i//save_every:04d}.png', oxresult = oxresult)
#        Dr.drawq(name = f'q{i//save_every:04d}.png', oxresult = oxresult)
        a=0
        nowtab = []
        ntab = []
        for index in range(n**2):
            if index%n==0:
                nowtab.append(oxnow[index])
                ntab.append(a)
                a+=1
        plt.plot(ntab, nowtab)
        plt.show()


    reg_reg_edges, reg_something_edges, in_edges=Pr.update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)
    reg_reg_edges, reg_something_edges, in_edges, oxresult=Ox.update_graph(oxnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)


