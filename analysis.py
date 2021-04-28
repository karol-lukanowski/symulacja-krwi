import numpy as np
import matplotlib.pyplot as plt
from build import dth
import networkx as nx
from morphopy.neurontree import NeuronTree as Nt
from morphopy.computation import file_manager as Fm


def getDiGraph(G, pnow, oxresult, in_nodes):
    pos = nx.get_node_attributes(G, 'pos')
    DiG = nx.DiGraph()
    for edge in G.edges:
        n1, n2 = edge
        d = G[n1][n2]['d']
        if oxresult[n1] and oxresult[n2] and d >= dth:
            deltaP = pnow[n1] - pnow[n2]
            if deltaP >= 0:
                DiG.add_edge(n1, n2)
            else:
                DiG.add_edge(n2, n1)

            DiG.nodes[n1]["pos"] = pos[n1]
            DiG.nodes[n2]["pos"] = pos[n2]

    if list(nx.simple_cycles(DiG)):
        print("petleeeee")
    return DiG


def createTestTree():
    G = nx.DiGraph()
    G.add_edges_from([(1, 2), (2, 3), (3,4), (4,5), (5,6), (7,2), (8,9), (9,3), (10,9), (11,4), (12,6), (13,14), (14,15), (15,5), (16,14), (17,18), (18, 15), (19,18), (18,14), (20,1)])
    pos = {20:(0,4), 1:(0,5), 2:(1,6), 3:(2,6), 4:(3,5), 5:(4,4), 6:(5,3), 7:(0,6), 8:(1,8), 9:(2,7), 10:(2,8), 11:(2,5), 12:(4,3), 13:(3,8), 14:(4,7), 15:(4,5), 16:(4,8), 17:(6,6), 18:(5,5), 19:(6,5)}
    
    for node in G.nodes:
        G.nodes[node]["pos"] = pos[node]
        G.nodes[node]["strahler"] = 0
    G = G.reverse()
    return G

def getStartingNodes(G):
    startingNodes = []
    for node in G.nodes:
        predecessors = G.predecessors(node)
        P = []
        for pre in predecessors:
            P.append(pre)
        if len(P) == 0:
            startingNodes.append(node)
    return startingNodes


def performStrahlerAlghorithm(node, G):
    # print(node)
    degree = G.out_degree(node)
    if degree == 0:
        G.nodes[node]["strahler"] = 1
        return G
    else:
        for subNode in G.successors(node):
            performStrahlerAlghorithm(subNode, G)
        subStrahlers = []
        for subNode in G.successors(node):
            subStrahlers.append(G.nodes[subNode]["strahler"])
            # print(subNode, G.nodes[subNode]["strahler"])

        maxStrahler = max(subStrahlers)
        if len(subStrahlers) == 2:            
            avrStrahlers = sum(subStrahlers)/2
            if maxStrahler == avrStrahlers:
                G.nodes[node]["strahler"] = maxStrahler + 1
            else:
                G.nodes[node]["strahler"] = maxStrahler
        elif len(subStrahlers) > 2:
            howmany = 0
            for sub in subStrahlers:
                if sub == maxStrahler:
                    howmany += 1
            if howmany > 1:
                G.nodes[node]["strahler"] = maxStrahler + 1
            else:
                G.nodes[node]["strahler"] = maxStrahler
        else:
            G.nodes[node]["strahler"] = maxStrahler
    return G

def strahlerOrder(G):
    startingNodes = getStartingNodes(G)
    # print(startingNodes)
    for node in startingNodes:
        G = performStrahlerAlghorithm(node, G)
    return G

def testStrahler():

    G = createTestTree()
    G = strahlerOrder(G)

    pos = nx.get_node_attributes(G, 'pos')
    strahler = nx.get_node_attributes(G, 'strahler')
    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_labels(G, pos, labels= strahler)
    # nx.draw_networkx_labels(G, pos, labels= {1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,11:11,12:12,13:13,14:14,15:15,16:16,17:17,18:18,19:19})
    plt.show()
def plotStrahlerHistogram(DiG, dirname):
    strahler = nx.get_node_attributes(DiG, 'strahler')

    S = []
    for value in strahler.values():
        S.append(value)
    d = 0.5
    left_of_first_bin = min(S) - float(d)/2
    right_of_last_bin = max(S) + float(d)/2
    fig, ax = plt.subplots()
    _ = ax.hist(S, np.arange(left_of_first_bin, right_of_last_bin + d, d))
    ax.set_xticks(np.arange(1,max(S)+1,1))
    
    plt.xlabel('Strahler number')
    plt.ylabel('Count')
    plt.savefig(dirname + "/strahler_histogram.png")

def plotStrahlerGraph(DiG, dirname):
    pos = nx.get_node_attributes(DiG, 'pos')
    strahler = nx.get_node_attributes(DiG, 'strahler')
    plt.figure(figsize=(100, 100))
    plt.axis('off')
    nx.draw_networkx_edges(DiG, pos)
    nx.draw_networkx_labels(DiG, pos, labels= strahler)
    plt.savefig(dirname + "/strahler_graph.png")
    # plt.show()


def getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname):
    DiG = getDiGraph(G, pnow, oxresult, in_nodes)
    DiG = strahlerOrder(DiG)
    
    plotStrahlerHistogram(DiG, dirname)
    
    plotStrahlerGraph(DiG, dirname)

# NETWORKTREE 

def setTypes(G, DiG, oxresult):
    soma = 0
    for n in DiG.nodes:
        neigh = G.neighbors(n)
        for n2 in neigh:
            d = G[n][n2]['d']
            neighbors = []
            if oxresult[n2] and d >= dth: 
                neighbors.append(n2)
            if len(neighbors) == 1:
                DiG.nodes[n]["type"] = 1
                soma +=1 
            else:
                DiG.nodes[n]["type"] = 2
    print("soma", soma)
    return DiG

def getDiGraphNetworkTree(G, pnow, oxresult, in_nodes):
    pos = nx.get_node_attributes(G, 'pos')
    DiG = nx.DiGraph()
    soma = 0
    for edge in G.edges:
        n1, n2 = edge
        d = G[n1][n2]['d']
        l = G[n1][n2]['length']
        if oxresult[n1] and oxresult[n2] and d >= dth:
            deltaP = pnow[n1] - pnow[n2]
            if deltaP >= 0:
                DiG.add_edge(n1, n2)
                DiG[n1][n2]["euclidean_dist"] = l
                DiG[n1][n2]["path_length"] = l
            else:
                DiG.add_edge(n2, n1)
                DiG[n2][n1]["euclidean_dist"] = l
                DiG[n2][n1]["path_length"] = l

            DiG.nodes[n1]["pos"] = [pos[n1][0],pos[n1][1]]
            DiG.nodes[n2]["pos"] = [pos[n2][0],pos[n2][1]]
            # if n1 in in_nodes:
            #     DiG.nodes[n1]["type"] = 1
            #     soma += 1
            # else:
            #     DiG.nodes[n1]["type"] = 2
            # if n2 in in_nodes:
            #     DiG.nodes[n2]["type"] = 1
            #     soma += 1
            # else:
            #     DiG.nodes[n2]["type"] = 2
            DiG.nodes[n1]["radius"] = 1
            DiG.nodes[n2]["radius"] = 1

    DiG = setTypes(G, DiG,oxresult)
    
                        
            
    posNew = nx.get_node_attributes(DiG, 'pos')
    print("ile nodesow", len(posNew))
    # print("ile soma", soma)

    # print(posNew)
    if list(nx.simple_cycles(DiG)):
        print("petleeeee")
    # plt.figure(figsize=(100, 100))
    # plt.axis('off')
    # nx.draw_networkx_edges(DiG, posNew, arrows = True)
    # plt.savefig("arrows.png")
    return DiG

def createNeuronTree(DiG):
    # tree = Fm.load_swc_file("sample_three_soma_points.swc")
    tree = Nt.NeuronTree(graph=DiG)
    print(len(tree.get_node_attributes("type")))
    from morphopy.neurontree.plotting import show_threeview
    fig = plt.figure(figsize=(10,10))
    show_threeview(tree, fig)
    # plt.show()
    hist, edges = tree.get_histogram(key='strahler_order')
    print(hist, edges)
    strahler = tree.get_strahler_order()
    print(strahler)
    print(len(strahler))
