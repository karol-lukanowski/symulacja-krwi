import networkx as nx
import numpy as np
import scipy.spatial

import utils as Ut


def Build_delaunay_net(n, periodic = 'top', diameterStats = ["uniform", 1, 1]):

    
    nkw = n**2

    points = np.random.uniform(0, n, (nkw, 2))
    points = np.array(sorted(points, key = lambda elem: (elem[0]//1, elem[1])))
    
    points_above = points.copy() + np.array([0, n])
    points_below = points.copy() + np.array([0, -n])
    points_right =  points.copy() + np.array([n, 0])
    points_left = points.copy() + np.array([-n, 0])

    if periodic == 'none':
       all_points = points
    elif periodic == 'top': 
        all_points = np.concatenate([points, points_above, points_below])
    elif periodic == 'side':
        all_points = np.concatenate([points, points_right, points_left])
    elif periodic == 'all':
        all_points = np.concatenate([points, points_above, points_below, points_right, points_left])


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
    boundary_edges = []

    final_edges_lengths = []
    for edge, l in zip(edges, edges_lengths):
        n1, n2 = edge
        if n2 < n1:
            n1, n2 = n2, n1

        if (n1 < nkw) and (n2 < nkw):
            final_edges.append((n1, n2))
            final_edges_lengths.append(l)
        elif (n1 < nkw) and (n2 >= nkw) and (n2 < 2*nkw):
            final_edges.append((n1, n2-nkw))
            boundary_edges.append((n1, n2-nkw))
            final_edges_lengths.append(l)
        elif (n1 < nkw) and (n2 >= 3 * nkw) and (n2 < 4*nkw):
            final_edges.append((n1, n2-3 * nkw))
            boundary_edges.append((n1, n2-3 * nkw))
            final_edges_lengths.append(l)


    nxGraph = nx.Graph()
    nxGraph.add_nodes_from(list(range(nkw)))

    nxGraph.add_edges_from(final_edges)

    for node in nxGraph.nodes:
        nxGraph.nodes[node]["pos"] = points[node]

    length_avr = 0
    for edge, l in zip(final_edges, final_edges_lengths):
        node, neigh = edge
        nxGraph[node][neigh]['length'] = l
        length_avr += l
        if diameterStats[0] == "uniform":
            nxGraph[node][neigh]['d'] = np.random.rand() * diameterStats[2] + diameterStats[1]
        elif diameterStats[0] == "gaussian":
            nxGraph[node][neigh]['d'] = Ut.PosGauss(diameterStats[1], diameterStats[2])
        elif diameterStats[0] == "lognormal":
            nxGraph[node][neigh]['d'] = np.random.lognormal(diameterStats[1], diameterStats[2])
        nxGraph[node][neigh]['q'] = 0

    length_avr /= len(nxGraph.edges())

    # Usunięcie zbyt długich krawędzi (szczególnie tych po brzegu)
    nxGraphCopy = nxGraph.copy()
    for node, neigh in nxGraphCopy.edges():
        l = nxGraph[node][neigh]['length']
        if l > 3 * length_avr:
            nxGraph.remove_edge(node, neigh)

    return nxGraph, boundary_edges






def find_center_node(nxGraph, n, xrange, yrange):
    x0, y0 = xrange / 2, yrange / 2
    pos0 = (x0, y0)
    Min = 100 * n
    for node in nxGraph.nodes:
        pos = nxGraph.nodes[node]["pos"]
        r = np.linalg.norm(np.array(pos) - np.array(pos0))
        if r < Min:
            Min = r
            id_center = node
    return id_center


def find_node(nxGraph, pos):
    def r_squared(node):
        x, y = nxGraph.nodes[node]['pos']
        r_sqr = (x - pos[0]) ** 2 + (y - pos[1]) ** 2
        return r_sqr
    r_min = len(nxGraph.nodes())
    n_min = 0
    for node in nxGraph.nodes():
        r = r_squared(node)
        if r < r_min:
            r_min = r
            n_min = node
    return n_min
    