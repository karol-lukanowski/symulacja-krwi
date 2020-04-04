# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 18:14:44 2020

@author: Antek
"""


def find_circle_nodes(G, n, R):
    pos = nx.get_node_attributes(G, 'pos')
    id_center = int(n*n/2)
    x0,y0 = G.node[id_center]["pos"]
    
    def cylindric(p,fi0,i):
        x,y = p
        r = np.sqrt((x-x0)**2 + (y-y0)**2) 
        if x-x0 > 0:
            if y==y0:
                fi = np.pi
            elif y - y0 < 0:
                fi = np.arctan((y-y0)/(x-x0))+2*np.pi - fi0
            else:
                fi = np.arctan((y-y0)/(x-x0)) +2*np.pi - fi0
        elif x-x0 < 0:
            if i > np.pi*R:
                fi = -np.arcsin((y-y0)/r) + 3*np.pi- fi0
            else:
                fi = -np.arcsin((y-y0)/r) + np.pi- fi0
        else:
            if i > np.pi*R:
                fi = np.pi*3/2 - fi0 + np.pi
            else:
                fi = np.pi/2 - fi0 + np.pi
        return (r,fi)
    def find_first_node():
        id = int(n/2)
        p= pos[id]
        r,_ = cylindric(p,0,0)
        while r > R:
            id += n
            p = pos[id]
            r,_ = cylindric(p,0,0)
        return id
    
    id_first = find_first_node()
    id = id_first
    p = G.node[id]["pos"]
    r_act, fi0 = cylindric(p,0,0)
    fi_act = 0
    
    circle = []
    circle.append(id)
    G.node[id]["circle"] = True
    flag_circle = False # flaga oznaczająca zrobienie pełnego okrążenia
    i = 0
    while flag_circle == False:
        i += 1
        #print(i)
        nei = G.node[id]["neigh"]
        flag_next = False # flaga oznaczająca znalezienie kolejnego punktu sposród sąsiadów
#        print(id, fi_act)
        for n in nei:
            p = G.node[n]["pos"]
            r_new, fi_new = cylindric(p,fi0,i)
            dfi = (fi_new - fi_act)
            # Czy zrobilismy już całe okrążenie?
           # print(p, x0,y0,fi_new)
            if n == id_first and i > 3:
                flag_circle = True
            # kolejny punkt ma większy kąt i miesci się w pierscieniu R-0.5 < r < R+0.5
            elif flag_next == False and  dfi > 0 and r_new > (R-0.5) and r_new < (R+0.5):
                id = n
                circle.append(id)
                G.node[id]["circle"] = True
                flag_next = True
                fi_act = fi_new
            #    print(n)
    return (x0,y0,circle)

                
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
def find_in_nodes(G):
    n0 = int(n*n/2)
    return search_nodes(G, n0)
def find_out_nodes(G):
    n0 = 1
    return search_nodes(G, n0)
