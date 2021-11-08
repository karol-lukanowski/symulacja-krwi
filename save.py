import networkx as nx
import dill
                                                                            
def save_all(name, sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges):
    pos = nx.get_node_attributes(G,'pos')
    All = [sid, edges, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, oxresult, boundary_edges, pos]

    with open(sid.dirname+name, 'wb') as file:
        dill.dump(All, file)


def load(name):
    with open(name, 'rb') as file:
        All= dill.load(file)
    
    sid = All[0]
    edges = All[1]
    in_nodes = All[2]
    out_nodes = All[3]
    in_nodes_ox = All[4]
    out_nodes_ox = All[5]
    oxresult = All[6]
    boundary_edges = All[7]
    pos = All[8]

    def reproduct():
        G1 = nx.Graph()
        new_pos = {}
        for key, value in pos.items():
            new_pos[int(key)] =  value
        
        for node in new_pos:
            G1.add_node(node, pos = new_pos[node])
        for n1, n2, d, l, t in edges:
            G1.add_edge(n1, n2, d = d, q = 0, length = l)

        return G1
    
    G1 = reproduct()
    return  sid, G1, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges