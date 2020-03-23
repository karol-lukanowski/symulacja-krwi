import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin


qin = 1  # ilosć wpływającej krwi
presout = 0  # cisnienie na wyjsciu
n = 10  # rozmiar siatki
mu = 0.0035  # współczynnik lepkosci
l = 1  # długosć krawędzi
c1 = np.pi / (128 * mu)  # stała przepływu
c2 = 64 * mu / (np.pi)  # stała siły
N = 1000  # liczba iteracji
length_wiggle_param = 0.3
diameter_wiggle_param = 0
qdrawconst = 10
ddrawconst = 10

# rysowanie przepływów
def drawq(G, name):
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, node_size=0, with_labels=False)
    
    qmax = 0
    for edge in G.edges(data='q'):
        if (edge[2] > qmax):
            qmax=edge[2]
    
    for edge in G.edges(data='q'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            nx.draw_networkx_edges(G, pos, edgelist=[edge], width=qdrawconst * edge[2] / qmax)
    
    plt.axis('equal')
    plt.savefig(name)
    plt.close()

# rysowanie srednic
def drawd(G, name):
    plt.figure(figsize=(10, 10))
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, node_size=0)
    
    dmax = 0
    for edge in G.edges(data='d'):
        if (edge[2] > dmax):
            dmax=edge[2]
    
    for edge in G.edges(data='d'):
        if not ((edge[0] % n == 0 and edge[1] % n == n - 1) or (edge[1] % n == 0 and edge[0] % n == n - 1)):
            nx.draw_networkx_edges(G, pos, edgelist=[edge], width=ddrawconst * edge[2] / dmax)
    #       print (edge[2])
    
    plt.axis('equal')
    plt.savefig(name)
    plt.close()

# czysto techniczne funkcje do warunków brzegowych
def boundaryr(i):
    if (i < n - 1):
        return i + 1
    else:
        return 0
def boundaryl(i):
    if (i > 0):
        return i - 1
    else:
        return n - 1


# sąsiedztwo: pierwsza i ostatnia kolumna siatki ma tylko po jednym połączeniu,
# żeby łatwiej trzymać stały wpływ i cisnienie; poza tym standardowe sąsiedztwo
# dla trójkatnej siatki
def neighbours(index):
    j = index % n
    i = int((index - j) / n)
    neightab = []
    if (i > 1) and (i < n - 2):
        neightab.append(n * i + boundaryr(j))
        neightab.append(n * i + boundaryl(j))
        if (i % 2 == 0):
            if (i > 0):
                neightab.append(n * (i - 1) + j)
                neightab.append(n * (i - 1) + boundaryl(j))
            if (i < n - 1):
                neightab.append(n * (i + 1) + j)
                neightab.append(n * (i + 1) + boundaryl(j))
        else:
            if (i > 0):
                neightab.append(n * (i - 1) + j)
                neightab.append(n * (i - 1) + boundaryr(j))
            if (i < n - 1):
                neightab.append(n * (i + 1) + j)
                neightab.append(n * (i + 1) + boundaryr(j))
    elif (i == 0):
        neightab.append(n * (i + 1) + j)
    elif (i == 1):
        neightab.append(n * (i - 1) + j)
        neightab.append(n * i + boundaryr(j))
        neightab.append(n * i + boundaryl(j))
        neightab.append(n * (i + 1) + j)
        neightab.append(n * (i + 1) + boundaryr(j))
    elif (i == n - 2):
        neightab.append(n * (i + 1) + j)
        neightab.append(n * i + boundaryr(j))
        neightab.append(n * i + boundaryl(j))
        if (i % 2 == 0):
            neightab.append(n * (i - 1) + j)
            neightab.append(n * (i - 1) + boundaryl(j))
        else:
            neightab.append(n * (i - 1) + j)
            neightab.append(n * (i - 1) + boundaryr(j))
    elif (i == n - 1):
        neightab.append(n * (i - 1) + j)
    return neightab




def create_pressure_flow_vector(qin,presout):
    """
    Wektor przepływów/cisnień, będący wynikiem działania macierzy na wektor cisnień
    pierwsze n elementów jest takie same, co odpowiada stałemu wpływowi do sieci
    dalej mamy elementy równe 0, co odpowiada sumie wpływów do każdego węzła równej 0
    ostatnie n elementów jest 0, bo ostatnie równania odpowiadają trzymaniu cisnienia
    wyjsciowego stale równego 0 (inaczej w układzie szybko pojawiają się ujemne cisnienia)
    """
    presult=np.zeros(n*n)
    for i in range(n*n):
        if (i<n):
            presult[i]=-qin
        elif (i>=n*(n-1)):
            presult[i]=presout
    return presult

def create_matrix(G):
    """
    macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
    z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
    srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
    """
    matrix=np.zeros((n*n,n*n))
    for node in G.nodes:
        if (node>=n and node<n*(n-1)):
            matrix[node][node]=0
            for ind in G.nodes[node]["neigh"]:
                d = G[node][ind]["d"]
                l = G[node][ind]['length']
                matrix[node][ind] = c1*d**4 /l
                matrix[node][node] -= c1*d**4 / l
        elif (node<n):
            matrix[node][node]=0
            for ind in range(n):
                d = G[ind][ind+n]["d"]
                l = G[ind][ind+n]['length']
                matrix[node][ind+n]=c1*d**4 /l
                matrix[node][node]-=c1*d**4 /l
        else:
            matrix[node][node]=1
    return matrix

def solve_equation_for_pressure(matrix, presult):
    """
    Zamieniamy macierz w równaniu na formę macierzy rzadkiej
    w celu usprawnienia obliczeń
    """
    S = spr.csc_matrix(matrix)
    pnow = sprlin.spsolve(S, presult)
    return pnow

def update_pressure_in_nodes(G, pnow):
    for node in G.nodes:
        G.nodes[node]["p"] = pnow[node]
    return G
def update_flow_in_edge(G, node1, node2, q):
    G[node][ind]['q'] = q
    return G

def update_diameter_in_edge(node1, node2, F):
    """
    Sposób, w jaki siła działa na srednice jest największym problemem, wynik bardzo mocno od tego zależy
    """
    if (F > 0.00002):
        if (F < 0.0002):
            G[node][ind]["d"] += 1000 * F
        else:
            G[node][ind]["d"] += 0.2
    return G
def find_edge_parameters(G, node1, node2):
    p1 = G.nodes[node]["p"]
    p2 = G.nodes[ind]["p"]
    d = G[node][ind]["d"]
    l = G[node][ind]['length']
    return (p1, p2, d, l)

def Flow(p1,p2,l,d):
    """
    Przepływ pomiędzy dwoma węzłami
    """
    dp = p1 - p2
    q = c1 * d**4 * dp/l
    return q

def Force(q,d):
    """
    Siła z jaką ciecz działa na scianki kanału
    """
    f = c2 * np.abs(q) / d**3
    return f

def add_nodes(G, n):   
    h = 3**0.5 / 2           # wysokosc trojkata równobocznego
    for i in range(n):
        for j in range(n):
            index=n*i+j
            if (i%2==0):
                pos=(i*h,j)
            else:
                pos=(i*h,j+0.5)
            if (i==0):
                G.add_node(index,p=1,pos=pos,neigh=neighbours(index))
            else:
                G.add_node(index,p=0,pos=pos,neigh=neighbours(index))
    return G
def wiggle_nodes(G, n, length_wiggle_param):
    """
    Przesuwanie wezlow zgodnie z length_wiggle_param, oprócz węzłów brzegowych
    """
    for node in G.nodes:
        if (node >= n and node < n*(n-1) and (node%n != 0) and ((node+1)%n != 0)):
            r = np.random.ranf() * 0.5
            fi = np.random.ranf() * 2 * np.pi
            dx = length_wiggle_param * r * np.cos(fi)
            dy = length_wiggle_param * r * np.sin(fi)
            pos = G.nodes[node]['pos']
            new_x, new_y = pos[0] + dx, pos[1] + dy
            G.nodes[node]['pos'] = (new_x, new_y)
    return G
def add_edges(G, n, diameter_wiggle_param, l):
    """
    Deklaracja połączeń: węzeł początkowy, węzeł końcowy, grubosć (z wylosowanym szumem), przepływ
    """
    for node in G.nodes:
        for ind in G.nodes[node]["neigh"]:
            if (ind>node):
                d0 = np.random.rand() * diameter_wiggle_param + 1
                G.add_edge(node,ind, d=d0,q=0, length=l)
    return G
def find_edges_lengths(G):
    for edge in G.edges():
        node, neigh = edge
        pos1 = G.nodes[node]['pos']
        pos2 = G.nodes[neigh]['pos']
    
        len = np.linalg.norm(np.array(pos1)-np.array(pos2))**0.5
        if (len <= 2.0): # czyli nie zmieniaj dlugosci miedzy brzegowymi wezlami, ona jest stale 1
            G[node][neigh]['length'] = len
    return G
    
def Build_triangular_net(n):
    G = nx.Graph()
    G = add_nodes(G, n)
    G = wiggle_nodes(G, n, length_wiggle_param)
    G = add_edges(G, n, diameter_wiggle_param, l)
    G = find_edges_lengths(G)
    return G
'''
W każdym kroku czasowym: 
-> obliczamy cisnienia na podstawie wzoru na przepływy między węzłami
oraz warunków brzegowych tj. stałe cisnienie na brzegach siatki + zerowanie się przepływów w węźle

-> Na podstawie nowych cisnień obliczamy przepływy w każdej krawędzi

-> Z przepływów obliczamy siłę działającą na scianki kanałów

-> Następnie rozszerzamy scianki na podstawie wyliczonej siły
'''
G = Build_triangular_net(n)
presult = create_pressure_flow_vector(qin,presout)

for i in range(N):
    print(f'{i+1}/{N}')
    matrix = create_matrix(G)    
    pnow = solve_equation_for_pressure(matrix,presult)
    G = update_pressure_in_nodes(G, pnow)
            
    for node in G.nodes:
        for ind in G.nodes[node]["neigh"]:
            if (ind > node):
                p1, p2, d, l = find_edge_parameters(G, node, ind)
                q = Flow(p1,p2,l,d)
                G = update_flow_in_edge(G,node,ind,q)             
                F = Force(q,d)
                G = update_diameter_in_edge(node,ind,F)
                
drawq(G, "q.png")
drawd(G, "d.png")

