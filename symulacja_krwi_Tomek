import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.sparse as spr
import scipy.sparse.linalg as sprlin

G=nx.Graph()

qin=1 #ilosć wpływającej krwi
presout=0 #cisnienie na wyjsciu
n=100 #rozmiar siatki
mu=0.0035 #współczynnik lepkosci
l=1 #długosć krawędzi - póki co trzymamy stałą długosć, a szum wprowadzamy zmianą początkowych grubosci
c1=np.pi/(128*mu*l) #stała przepływu
c2=64*mu/(np.pi) #stała siły
N=10 #liczba iteracji

#czysto techniczne funkcje do warunków brzegowych
def boundaryr(i):
    if (i<n-1):
        return i+1
    else:
        return 0

def boundaryl(i):
    if (i>0):
        return i-1
    else:
        return n-1
        

#sąsiedztwo: pierwsza i ostatnia kolumna siatki ma tylko po jednym połączeniu,
#żeby łatwiej trzymać stały wpływ i cisnienie; poza tym standardowe sąsiedztwo
#dla trójkatnej siatki        
def neighbours(index):
    j=index%n
    i=int((index-j)/n)
    neightab=[]
    if (i>1) and (i<n-2):
        neightab.append(n*i+boundaryr(j))
        neightab.append(n*i+boundaryl(j))
        if (i%2==0):
            if (i>0):
                neightab.append(n*(i-1)+j)
                neightab.append(n*(i-1)+boundaryl(j))
            if (i<n-1):
                neightab.append(n*(i+1)+j)
                neightab.append(n*(i+1)+boundaryl(j))
        else:
            if (i>0):
                neightab.append(n*(i-1)+j)
                neightab.append(n*(i-1)+boundaryr(j))
            if (i<n-1):
                neightab.append(n*(i+1)+j)
                neightab.append(n*(i+1)+boundaryr(j))
    elif (i==0):
        neightab.append(n*(i+1)+j)
    elif (i==1):
        neightab.append(n*(i-1)+j)
        neightab.append(n*i+boundaryr(j))
        neightab.append(n*i+boundaryl(j))
        neightab.append(n*(i+1)+j)
        neightab.append(n*(i+1)+boundaryr(j))
    elif (i==n-2):
        neightab.append(n*(i+1)+j)
        neightab.append(n*i+boundaryr(j))
        neightab.append(n*i+boundaryl(j))    
        if (i%2==0):
            neightab.append(n*(i-1)+j)
            neightab.append(n*(i-1)+boundaryl(j))
        else:
            neightab.append(n*(i-1)+j)
            neightab.append(n*(i-1)+boundaryr(j))
    elif (i==n-1):
        neightab.append(n*(i-1)+j)
    return neightab
        
#deklaracja węzłów: kolejno indeks węzła, cisnienie, pozycja, sąsiedztwo
for i in range(n):
    for j in range(n):
        index=n*i+j
        if (i%2==0):
            pos=(i,j)
        else:
            pos=(i,j+0.5)
        if (i==0):
            G.add_node(index,p=1,pos=pos,neigh=neighbours(index))
        else:
            G.add_node(index,p=0,pos=pos,neigh=neighbours(index))

#deklaracja połączeń: węzeł początkowy, węzeł końcowy, grubosć (z wylosowanym szumem), przepływ
for node in G.node:
    for ind in G.node[node]["neigh"]:
        if (ind>node):
            d0=np.random.rand()*0.1+1
            G.add_edge(node,ind, d=d0,q=0)

#wektor przepływów/cisnień, będący wynikiem działania macierzy na wektor cisnień
#pierwsze n elementów jest takie same, co odpowiada stałemu wpływowi do sieci
#dalej mamy elementy równe 0, co odpowiada sumie wpływów do każdego węzła równej 0
#ostatnie n elementów jest 0, bo ostatnie równania odpowiadają trzymaniu cisnienia
#wyjsciowego stale równego 0 (inaczej w układzie szybko pojawiają się ujemne cisnienia)
presult=np.zeros(n*n)
for i in range(n*n):
    if (i<n):
        presult[i]=-qin
    elif (i>=n*(n-1)):
        presult[i]=presout
       
#macierz przepływów/cisnień; pierwsze n wierszy odpowiada za utrzymanie stałego wpływu:
#z kolei ostatnie n wierszy utrzymuje wyjsciowe cisnienie równe 0
#srodkowe wiersze utrzymują sumę wpływów/wypływów w każdym węźle równa 0
matrix=np.zeros((n*n,n*n))
for i in range(N):
    print (i)
    for node in G.node:
        if (node>=n and node<n*(n-1)):
            matrix[node][node]=0
            for ind in G.node[node]["neigh"]:
                matrix[node][ind]=c1*(G[node][ind]["d"])**4
                matrix[node][node]-=c1*(G[node][ind]["d"])**4
        elif (node<n):
            matrix[node][node]=0
            for ind in range(n):
                matrix[node][ind+n]=c1*(G[ind][ind+n]["d"])**4
                matrix[node][node]-=c1*(G[ind][ind+n]["d"])**4
        else:
            matrix[node][node]=1

    #do odkomentowania liczenia na macierzach rzadkich, póki co dobrze działa normalne rozwiązywania równania macierzowego
    S = spr.csc_matrix(matrix)
    pnow=sprlin.spsolve(S,presult)
#    pnow=np.linalg.solve(matrix,presult)

    #aktualizujemy cisnienia w siatce
    for node in G.node:
        G.node[node]["p"]=pnow[node]
        #na podstawie nowych cisnień obliczamy przepływy w każdej krawędzi, a z nich siłę
    for node in G.node:
        for ind in G.node[node]["neigh"]:
            if (ind>node):
                G[node][ind]["q"]=c1*(G[node][ind]["d"])**4*(G.node[node]["p"]-G.node[ind]["p"])
                F=c2*np.abs(G[node][ind]["q"])/(G[node][ind]["d"])**3
                #sposób, w jaki siła działa na srednice jest największym problemem, wynik bardzo mocno od tego zależy    
                if (F>0.00002):
                    if (F<0.0002):
                        G[node][ind]["d"]+=1000*F
                    else:
                        G[node][ind]["d"]+=0.2
                   
                        
                
                


#for edge in G.edges(data='q'):
#    if (edge[1]>=90):
#        print (edge[2])

#for edge in G.edges(data='d'):
#    print (edge[2])

#rysowanie przepływów - ich wartosci są małe, więc grubosci rysowanych krawędzi trzeba przemnożyć
plt.figure(figsize=(10,10)) 
pos = nx.get_node_attributes(G,'pos')
nx.draw(G,pos,node_size=0)

for edge in G.edges(data='q'):
    if not((edge[0]%n==0 and edge[1]%n==n-1) or (edge[1]%n==0 and edge[0]%n==n-1)):
        nx.draw_networkx_edges(G, pos, edgelist=[edge], width=100*edge[2])
#        print (edge[2])
    
plt.savefig("graphq.png")

#rysowanie srednic
plt.figure(figsize=(10,10)) 
pos = nx.get_node_attributes(G,'pos')
nx.draw(G,pos,node_size=0)

for edge in G.edges(data='d'):
    if not((edge[0]%n==0 and edge[1]%n==n-1) or (edge[1]%n==0 and edge[0]%n==n-1)):
        nx.draw_networkx_edges(G, pos, edgelist=[edge], width=edge[2])
#       print (edge[2])
    
plt.savefig("graphd.png")
