import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.spatial

def Build_delaunay_net(n, diameter_wiggle_param=1):
    N = n**2


    points = np.random.uniform(0, n, (N, 2))
    points = np.array(sorted(points, key = lambda elem: (elem[0]//1, elem[1])))

    points_above = points.copy() + np.array([0, n])
    points_below = points.copy() + np.array([0, -n])


    all_points = np.concatenate([points, points_above, points_below])

    delTri = scipy.spatial.Delaunay(all_points)

    # create a set for edges that are indexes of the points
    edges = set()
    # for each Delaunay triangle
    for node in range(delTri.nsimplex):
        # for each edge of the triangle
        # sort the vertices
        # (sorting avoids duplicated edges being added to the set)
        # and add to the edges set
        edge = sorted([delTri.vertices[node, 0], delTri.vertices[node, 1]])
        edges.add((int(edge[0]), int(edge[1])))
        edge = sorted([delTri.vertices[node, 0], delTri.vertices[node, 2]])
        edges.add((int(edge[0]), int(edge[1])))
        edge = sorted([delTri.vertices[node, 1], delTri.vertices[node, 2]])
        edges.add((int(edge[0]), int(edge[1])))

    edges = list(edges)
    edges_lengths = []

    for edge in edges:
        n1, n2 = edge
        pos1, pos2 = all_points[n1], all_points[n2]
        l = np.linalg.norm(np.array(pos1) - np.array(pos2))
        edges_lengths.append(l)

    # now choose edges between "points" and take care of the boundary conditions (edges between points and points_above)
    # points are indexes 0:(N-1), points_above are N:(2N-1)

    
    final_edges = []
    dontdraw_edges = []

    final_edges_lengths = []
    for edge, l in zip(edges, edges_lengths):
        n1, n2 = edge
        if n2 < n1:
            n1, n2 = n2, n1

        if (n1 < N) and (n2 < N):
            final_edges.append((n1, n2))
            final_edges_lengths.append(l)

        elif (n1 < N) and (n2 >= N) and (n2 < 2*N):
            final_edges.append((n1, n2-N))

            dontdraw_edges.append((n1, n2-N))

            final_edges_lengths.append(l)


    G = nx.Graph()
    G.add_nodes_from(list(range(N)))

    G.add_edges_from(final_edges)

    for node in G.nodes:
        G.nodes[node]["pos"] = points[node]

    length_avr = 0
    for edge, l in zip(final_edges, final_edges_lengths):
        node, neigh = edge
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

    return G, dontdraw_edges, "de"






def find_center_node(G, n, xrange, yrange):
    x0, y0 = xrange / 2, yrange / 2
    pos0 = (x0, y0)
    Min = 100 * n
    for node in G.nodes:
        pos = G.nodes[node]["pos"]
        r = np.linalg.norm(np.array(pos) - np.array(pos0))
        if r < Min:
            Min = r
            id_center = node
    return id_center