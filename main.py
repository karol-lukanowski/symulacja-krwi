import numpy as np

import draw_net as Dr
import pressure as Pr
import oxygen as Ox
import vegf as Ve
from config import n, G, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, iters, in_edges, c1, save_every, dth
import matplotlib.pyplot as plt
import time

reg_reg_edges, reg_something_edges, other_edges = [],  [], []
for n1, n2 in G.edges():
    d = G[n1][n2]['d']
    l = G[n1][n2]['length']
    if (n1 not in in_nodes and n1 not in out_nodes) and (n2 not in in_nodes and n2 not in out_nodes):
        reg_reg_edges.append((n1, n2, d, l))
    elif (n1 not in in_nodes and n1 not in out_nodes):
        reg_something_edges.append((n1, n2, d, l))
    elif (n2 not in in_nodes and n2 not in out_nodes):
        reg_something_edges.append((n2, n1, d, l ))
    else:
        other_edges.append((n1, n2, d, l))

presult = Pr.create_vector()
oxresult = Ox.create_vector()

#oxtype = oxresult.copy()
#oxtype[out_nodes_ox] = 1.


for i in range(iters):
    print(f'Iter {i + 1}/{iters}')

    t1 = time.time()

    pmatrix = Pr.update_matrix(reg_reg_edges, reg_something_edges, in_edges)

    t1_5 = time.time()

    oxmatrix = Ox.update_matrix(oxresult, reg_reg_edges, reg_something_edges, other_edges)

    t2 = time.time()

    pnow = Pr.solve_equation(pmatrix, presult)

    t2_5 = time.time()

    oxnow = Ox.solve_equation(oxmatrix, oxresult)

    t3 = time.time()

    #vresult = Ve.create_vector(oxnow, oxresult)
    #vmatrix = Ve.update_matrix(vresult, reg_reg_edges, reg_something_edges)
    #vnow = Ve.solve_equation(vmatrix, vresult)

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


        #Dr.drawhist(name = f'{i//save_every:04d}.png', oxnow = oxnow, oxresult = oxresult, vnow = oxnow)
        #Dr.drawd(name = f'd{i//save_every:04d}.png', oxresult = oxresult)
        Dr.drawq(name = f'q{i//save_every:04d}.png', oxdraw = oxnow)



        """
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
        """

    reg_reg_edges, reg_something_edges, in_edges=Pr.update_graph(pnow, reg_reg_edges, reg_something_edges, in_edges)



    #reg_reg_edges, reg_something_edges, in_edges, oxresult=Ve.update_graph(vnow, oxresult, reg_reg_edges, reg_something_edges, in_edges)

    t4 = time.time()




    for e in reg_reg_edges+reg_something_edges+in_edges:
        n1, n2, d, l = e
        if (oxresult[n1] == 1 or oxresult[n2] == 1):
            if d > dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
                #oxtype[n1] = 1
                #oxtype[n2] = 1



    print(f'update mat press {t1_5-t1}, update mat ox {t2-t1_5}, solve pres eq {t2_5-t2}, solve ox eq {t3-t2_5}, update graph {t4-t3}')

