import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.spatial

def Build_delaunay_net(N, diameter_wiggle_param=1):
    points = []
    for i in range(N):
        x = np.random.rand()*np.sqrt(N)
        y = np.random.rand()*np.sqrt(N)        
        points.append((x,y)) 
    
    # make a Delaunay triangulation of the point data
    delTri = scipy.spatial.Delaunay(points)
    
    # create a set for edges that are indexes of the points
    edges = set()
    # for each Delaunay triangle
    for n in range(delTri.nsimplex):
        # for each edge of the triangle
        # sort the vertices
        # (sorting avoids duplicated edges being added to the set)
        # and add to the edges set
        edge = sorted([delTri.vertices[n,0], delTri.vertices[n,1]])
        edges.add((edge[0], edge[1]))
        edge = sorted([delTri.vertices[n,0], delTri.vertices[n,2]])
        edges.add((edge[0], edge[1]))
        edge = sorted([delTri.vertices[n,1], delTri.vertices[n,2]])
        edges.add((edge[0], edge[1]))
    
    
    # make a graph based on the Delaunay triangulation edges
    G = nx.Graph(list(edges))
    
    for node in G.nodes:
        G.node[node]["pos"]= points[node]

    def find_edges_lengths_and_diameters():
        length_avr = 0
        for edge in G.edges():
            node, neigh = edge
            pos1 = G.nodes[node]['pos']
            pos2 = G.nodes[neigh]['pos']
        
            l = np.linalg.norm(np.array(pos1)-np.array(pos2))
            G[node][neigh]['length'] = l
            length_avr += l
            G[node][neigh]['d'] = np.random.rand() * diameter_wiggle_param + 1
        length_avr /= len(G.edges())
        #Usunięcie zbyt długich krawędzi (szczególnie tych po brzegu)
        Gcopy = G.copy()
        for edge in Gcopy.edges():
            node, neigh = edge
            l = G[node][neigh]['length']
            if l > 3*length_avr:
                G.remove_edge(node, neigh)
    find_edges_lengths_and_diameters()
    '''
    plt.figure(figsize=(10, 10))
    nx.draw(G,pos = points,node_size=30)
    plt.savefig("xd1", dpi=150)
    
    plt.show()
    '''
    return G
'''
N = 1000 
Build_delaunay_net(N)
'''