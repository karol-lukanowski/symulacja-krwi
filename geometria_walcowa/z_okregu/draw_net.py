import matplotlib.pyplot as plt
import networkx as nx

def drawq(G, n, name, qdrawconst = 5, normalize = True):
    """
    rysowanie przepływów
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    #nx.draw(G, pos, node_size=0, with_labels=False)
    
    qmax = 0
    for edge in G.edges(data='q'):
        if (edge[2] > qmax):
            qmax=edge[2]

    if (normalize == False):
        qmax = 1
        qdrawconst = 1
    
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            nx.draw_networkx_edges(G, pos, edgelist=[edge], width=qdrawconst * edge[2] / qmax)
    
    plt.axis('equal')
    plt.savefig(name)
    plt.close()

def drawd(G, n, name, ddrawconst = 3, normalize = True):
    """
    rysowanie srednic
    """
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    #nx.draw(G, pos, node_size=0)
    
    dmax = 0
    for edge in G.edges(data='d'):
        if (edge[2] > dmax):
            dmax=edge[2]

    if (normalize == False):
        dmax = 1
        ddrawconst = 1

    for edge in G.edges(data='d'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            nx.draw_networkx_edges(G, pos, edgelist=[edge], width=ddrawconst * edge[2] / dmax)
    
    plt.axis('equal')
    plt.savefig(name)
    plt.close()
