import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.spatial


def Build_delaunay_net(n, diameter_wiggle_param=1):
    N = n ** 2
    points = np.random.uniform(0, n, (N,2))

    # sortowanie wezlow = najpierw po "kolumnie" od 0 do n, potem po wierszu - tak jak nasza prostokatna geometria
    points = sorted(points, key = lambda elem: (elem[0]//1, elem[1]))

    # make a Delaunay triangulation of the point data
    delTri = scipy.spatial.Delaunay(points)

    # create a set for edges that are indexes of the points
    edges = set()
    # for each Delaunay triangle
    for node in range(delTri.nsimplex):
        # for each edge of the triangle
        # sort the vertices
        # (sorting avoids duplicated edges being added to the set)
        # and add to the edges set
        edge = sorted([delTri.vertices[node, 0], delTri.vertices[node, 1]])
        edges.add((edge[0], edge[1]))
        edge = sorted([delTri.vertices[node, 0], delTri.vertices[node, 2]])
        edges.add((edge[0], edge[1]))
        edge = sorted([delTri.vertices[node, 1], delTri.vertices[node, 2]])
        edges.add((edge[0], edge[1]))

    # make a graph based on the Delaunay triangulation edges
    G = nx.Graph(list(edges))

    for node in G.nodes:
        G.nodes[node]["pos"] = points[node]

    def find_edges_lengths_and_diameters():
        length_avr = 0
        for node, neigh in G.edges():
            pos1 = G.nodes[node]['pos']
            pos2 = G.nodes[neigh]['pos']

            l = np.linalg.norm(np.array(pos1) - np.array(pos2))
            G[node][neigh]['length'] = l
            length_avr += l
            G[node][neigh]['d'] = np.random.rand() * diameter_wiggle_param + 1
            G[node][neigh]['q'] = 0

        length_avr /= len(G.edges())

        # Usunięcie zbyt długich krawędzi (szczególnie tych po brzegu)
        Gcopy = G.copy()
        for node, neigh in Gcopy.edges():
            l = G[node][neigh]['length']
            if l > 3 * length_avr:
                G.remove_edge(node, neigh)

    def find_center_node():
        x0 = y0 = n / 2
        pos0 = (x0, y0)
        Min = 100 * n
        for node in G.nodes:
            pos = G.nodes[node]["pos"]
            r = np.linalg.norm(np.array(pos) - np.array(pos0))
            if r < Min:
                Min = r
                id_center = node
        return id_center

    find_edges_lengths_and_diameters()

    """
    plt.figure(figsize=(7, 7))
    nx.draw(G,pos = points,node_size=10,with_labels=False)
    id_center = find_center_node()
    nx.draw_networkx(G, pos = points, nodelist = [id_center], node_size=50, node_color='r', with_labels = False)
    plt.show()
    """

    return G
