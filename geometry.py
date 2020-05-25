# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import delaunay as De


# Funkcje do znajdowania wezlow na okregu
def find_circle_nodes(G, n, R):
    pos = nx.get_node_attributes(G, 'pos')
    id_center = int(n * n / 2)
    x0, y0 = G.nodes[id_center]["pos"]

    def cylindric(p, fi0, i):
        x, y = p
        r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        if x - x0 > 0:
            if y == y0:
                fi = np.pi
            elif y - y0 < 0:
                fi = np.arctan((y - y0) / (x - x0)) + 2 * np.pi - fi0
            else:
                fi = np.arctan((y - y0) / (x - x0)) + 2 * np.pi - fi0
        elif x - x0 < 0:
            if i > np.pi * R:
                fi = -np.arcsin((y - y0) / r) + 3 * np.pi - fi0
            else:
                fi = -np.arcsin((y - y0) / r) + np.pi - fi0
        else:
            if i > np.pi * R:
                fi = np.pi * 3 / 2 - fi0 + np.pi
            else:
                fi = np.pi / 2 - fi0 + np.pi
        return (r, fi)

    def find_first_node():
        id = int(n / 2)
        p = pos[id]
        r, _ = cylindric(p, 0, 0)
        while r > R:
            id += n
            p = pos[id]
            r, _ = cylindric(p, 0, 0)
        return id

    id_first = find_first_node()
    id = id_first
    p = G.nodes[id]["pos"]
    r_act, fi0 = cylindric(p, 0, 0)
    fi_act = 0

    circle = []
    circle.append(id)
    G.nodes[id]["circle"] = True
    flag_circle = False  # flaga oznaczająca zrobienie pełnego okrążenia
    i = 0
    while flag_circle == False:
        i += 1
        # print(i)
        nei = G.nodes[id]["neigh"]
        flag_next = False  # flaga oznaczająca znalezienie kolejnego punktu sposród sąsiadów
        #        print(id, fi_act)
        for n in nei:
            p = G.nodes[n]["pos"]
            r_new, fi_new = cylindric(p, fi0, i)
            dfi = (fi_new - fi_act)
            # Czy zrobilismy już całe okrążenie?
            # print(p, x0,y0,fi_new)
            if n == id_first and i > 3:
                flag_circle = True
            # kolejny punkt ma większy kąt i miesci się w pierscieniu R-0.5 < r < R+0.5
            elif flag_next == False and dfi > 0 and r_new > (R - 0.5) and r_new < (R + 0.5):
                id = n
                circle.append(id)
                G.nodes[id]["circle"] = True
                flag_next = True
                fi_act = fi_new
            #    print(n)
    return (x0, y0, circle)
def search_nodes(G, n0):
    nodes = []

    def rec_search(n):
        if G.node[n]["visited"] == True:
            return
        elif G.node[n]["circle"] == True:
            return
        else:
            G.node[n]["visited"] = True
            nodes.append(n)
            neighbours = G.node[n]["neigh"]
            for nei in neighbours:
                if G.node[nei]["visited"] == False:
                    rec_search(nei)

    rec_search(n0)
    return nodes


def set_geometry(n, G=[], geo='rect', R=25, R_s=5, *args, **kwargs):
    def rect_default_nodes():
        in_nodes = list(range(n))
        reg_nodes = list(range(n,n * (n - 1)))
        out_nodes = list(range(n * (n - 1), n * n))
        return in_nodes, out_nodes, reg_nodes
    def cyl_default_nodes():
        id_center = n * n // 2
        if 'del' in kwargs and kwargs['del'] == True:
            id_center = De.find_center_node(G, n, xrange=n, yrange=n)
        in_nodes = [id_center]
        x0, y0 = G.nodes[id_center]["pos"]
        out_nodes = []
        reg_nodes = []
        for node in G.nodes:
            pos = G.nodes[node]["pos"]
            r = np.sqrt((pos[0] - x0) ** 2 + (pos[1] - y0) ** 2)
            if r > R and r < R + 1:
                out_nodes.append(node)
            elif r < R+1:
                reg_nodes.append(node)

        return in_nodes, out_nodes, reg_nodes
    def don_default_nodes():
        id_center = n * n // 2
        if 'del' in kwargs and kwargs['del'] == True:
            id_center = De.find_center_node(G, n, xrange=n, yrange=n)
        x0, y0 = G.nodes[id_center]["pos"]
        in_nodes = []
        out_nodes = []
        reg_nodes = []
        for node in G.nodes:
            pos = G.nodes[node]["pos"]
            r = np.sqrt((pos[0] - x0) ** 2 + (pos[1] - y0) ** 2)
            if r > R and r < R + 1:
                out_nodes.append(node)
            elif r > R_s and r < R_s + 1:
                in_nodes.append(node)
            elif r < R+1 and r > R_s:
                reg_nodes.append(node)
        return in_nodes, out_nodes, reg_nodes

    def don_default_nodes2():
        in_nodes = find_circle_nodes(G, n, R_s)[2]
        out_nodes = find_circle_nodes(G, n, R)[2]
        return in_nodes, out_nodes

    in_nodes, out_nodes, reg_nodes, in_edges = [], [], [], []
    if geo == 'rect':
        in_nodes, out_nodes, reg_nodes = rect_default_nodes()
    elif geo == 'cylindrical':
        in_nodes, out_nodes, reg_nodes = cyl_default_nodes()
    elif geo == 'donut':
        in_nodes, out_nodes, reg_nodes = don_default_nodes()
    elif geo == 'own':
        in_nodes, out_nodes = kwargs['in_nodes'], kwargs['out_nodes']
        reg_nodes = [node for node in G.nodes() if (node not in in_nodes) and (node not in out_nodes)]
    else:
        print('Wrong geometry specified')
        return 0


#    reg_nodes = [node for node in G.nodes() if (node not in in_nodes) and (node not in out_nodes)]
    for node in in_nodes:
        for neigh in G.neighbors(node):
            if neigh not in in_nodes:
                d = G[node][neigh]['d']
                l = G[node][neigh]['length']
                in_edges.append((node, neigh, d, l))
    
    other_nodes = []           
    for node in G.nodes():
        if (node not in reg_nodes) and (node not in in_nodes) and (node not in out_nodes):
            other_nodes.append(node)
                
    return in_nodes, out_nodes, reg_nodes, other_nodes, in_edges



def equidistant_geometry(G, n, R, xrange, yrange, how_many):
    id_center = De.find_center_node(G, n, xrange=xrange, yrange=yrange)

    def r_squared(node):
        # x0, y0 = G.nodes[n*n//2]["pos"]
        x0, y0 = G.nodes[id_center]['pos']
        x, y = G.nodes[node]['pos']
        r_sqr = (x - x0) ** 2 + (y - y0) ** 2
        return r_sqr

    boundary_nodes = []
    for (n1, n2) in G.edges():
        r1, r2 = r_squared(n1), r_squared(n2)
        if r1 > r2:
            n1, n2 = n2, n1
            r1, r2 = r2, r1

        n_b = n2

        if r2 >= R ** 2 and r1 <= R ** 2:
            # x, y = G.nodes[n_b]['pos'][0] - G.nodes[n**2 // 2]['pos'][0], G.nodes[n_b]['pos'][1] - G.nodes[n**2 // 2]['pos'][1]
            x, y = G.nodes[n_b]['pos'][0] - G.nodes[id_center]['pos'][0], G.nodes[n_b]['pos'][1] - \
                   G.nodes[id_center]['pos'][1]

            if x == 0: x = 0.000001
            if y == 0: y = 0.000001

            if (x >= 0 and y >= 0):
                fi = np.arctan(y / x)
            elif (x < 0 and y >= 0):
                fi = np.pi / 2 + np.arctan(-x / y)
            elif (x < 0 and y < 0):
                fi = np.pi + np.arctan(y / x)
            else:
                fi = (3 / 2) * np.pi + np.arctan(x / -y)
            boundary_nodes.append([n_b, fi])
    boundary_nodes.sort(key=lambda node: node[1])

    boundary_nodes, fis = zip(*boundary_nodes)

    num_of_out_nodes = how_many
    out_indexes = np.round(np.linspace(0, len(boundary_nodes) - 1, num_of_out_nodes + 1)).astype(int)
    out_nodes = list(np.array(boundary_nodes)[out_indexes[:-1]])
    in_nodes = [id_center]

    return in_nodes, out_nodes

def create_edgelist(G, in_nodes, out_nodes, reg_nodes):
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
       
    removetab_reg = []
    
    for index, (n1, n2, d, l) in enumerate(reg_reg_edges):
        if n1 not in reg_nodes and n2 in reg_nodes:
            removetab_reg.append(index)
        elif n1 in reg_nodes and n2 not in reg_nodes:
            removetab_reg.append(index)
    
    removetab_something = []

    for index, (n1, n2, d, l) in enumerate(reg_something_edges):
        if n1 not in reg_nodes and n2 not in in_nodes:
            removetab_something.append(index)
    
    regnew = []
    somethingnew = []
    
    for index, (n1, n2, d, l) in enumerate(reg_reg_edges):
        if index not in removetab_reg:
            regnew.append((n1, n2, d, l))
        
    for index, (n1, n2, d, l) in enumerate(reg_something_edges):
        if index not in removetab_something:
            somethingnew.append((n1, n2, d, l))
    
    reg_reg_edges = regnew
    reg_something_edges = somethingnew

    return reg_reg_edges, reg_something_edges, other_edges