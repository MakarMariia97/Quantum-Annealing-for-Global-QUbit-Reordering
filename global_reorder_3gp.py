import networkx as nx
import math
import sys
from collections import Counter
import time
from utils import solve_QUBO, checkBalance, makePermutation


def graphPart(G, nparts, vdegree, total_ordering, trials,sim=False):
    print_graph(G)
    k=3
    if (nparts>1):
        nNodes = len(G.nodes)
        is_ok = [-1]*nNodes
        trials -= 1
        unbal=-1
        while (sum(is_ok)!=0 and unbal!=0):
            trials += 1
            sample = solve_QUBO(G, vdegree, k, sim)
            for i in range(nNodes):
                if (sample[i]+sample[nNodes+i]+sample[2*nNodes+i] == 1):
                    is_ok[i] = 0
                else:
                    is_ok[i] = -1
            part3 = makePermutation(G, nNodes, k, sample)
            unbal = checkBalance(G, nNodes, part3, k)
            
        # partition graph into 3 parts
        lG = nx.Graph()
        mG = nx.Graph()
        rG = nx.Graph()
        for u in range(len(G.nodes)):
            if part3[u] == 0:
                lG.add_node(list(G.nodes)[u])
            if part3[u] == 1:
                mG.add_node(list(G.nodes)[u])
            if part3[u] == 2:
                rG.add_node(list(G.nodes)[u])
        for u, v, d in G.edges(data=True):

            if part3[list(G.nodes).index(u)] == 0 and part3[list(G.nodes).index(v)] == 0:
                lG.add_edge(u, v, weight = d["weight"])
            if part3[list(G.nodes).index(u)] == 1 and part3[list(G.nodes).index(v)] == 1:
                mG.add_edge(u, v,  weight = d["weight"])
            if part3[list(G.nodes).index(u)] == 2 and part3[list(G.nodes).index(v)] == 2:
                rG.add_edge(u, v,  weight = d["weight"])

        if (len(list(lG.nodes)) < k and len(list(mG.nodes)) < k and len(list(rG.nodes)) < k):
            return (list(lG.nodes) + list(mG.nodes) + list(rG.nodes),trials)

        if (len(list(lG.nodes())) >= k):
            total_ord_from_part_lG, trial_lG = graphPart(
                lG, math.floor(nparts/k), list(dict(lG.degree(weight="weight")).values()), total_ordering, trials)
            total_ordering = total_ordering + total_ord_from_part_lG
            trials = trial_lG
        if (len(list(mG.nodes())) >= k):
            total_ord_from_part_mG, trial_mG = graphPart(
                mG, math.floor(nparts/k), list(dict(mG.degree(weight="weight")).values()), total_ordering, trials)
            total_ordering = total_ordering + total_ord_from_part_mG
            trials = trial_mG
        if (len(rG.nodes) >= k):
            total_ord_from_part_rG, trial_rG = graphPart(
                rG, math.floor(nparts/k), list(dict(rG.degree(weight="weight")).values()), [], trials)
            total_ordering = total_ordering + total_ord_from_part_rG
            trials = trial_rG

    return (total_ordering, trials)

def print_graph(G):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    print("nodes", G.nodes)
    print("edges", G.edges)

def calculateSWAP(total_ordering, gates):
    sum=0
    for u, v, w in gates:
        ind_v = list(total_ordering).index(v)
        ind_w = list(total_ordering).index(w)
        if (u != -1):
            ind_u = list(total_ordering).index(u)
            a = ind_u
            b = ind_v
            ind_u = min(a, b)
            ind_v = max(a, b)
            if (ind_w < ind_u & ind_u < ind_v):
                if (ind_u-ind_w>1):
                    sum+= abs(ind_u - ind_w) - 1 + abs(ind_v - ind_u) - 1
                elif (ind_v - ind_u>1):
                        sum+=ind_v - ind_u - 1
            elif (ind_u < ind_w & ind_w < ind_v):
                if (ind_w-ind_u>1):
                    sum+= abs(ind_w - ind_u) - 1
                if (ind_v - ind_w>1):
                        sum+=ind_v - ind_w - 1
            elif (ind_u < ind_v & ind_v < ind_w):
                if (ind_w - ind_v >1):
                    sum+= abs(ind_w - ind_v) - 1 + abs(ind_v - ind_u) - 1
                elif (ind_v - ind_u>1):
                    sum+= abs(ind_v - ind_u) - 1
        else:
            sum+=abs(ind_v-ind_w) - 1
    return sum

def main():
    start_time = time.time()
    
  #  orig_stdout = sys.stdout
  #  f = open('SA_MULT_19_13_out_graph_4part_gamma96.txt', 'w')
  #  sys.stdout = f

    nq=6 # dTof
   # nq=8 # ham7
  #  nq=4
  #  nq=12
    sim=True

    G = nx.Graph()
    G.add_nodes_from(range(nq))
    #G.add_weighted_edges_from([(0, 1, 1), (1, 2, 1), (1, 3, 2)])
    #G.add_weighted_edges_from([(0, 5, 1), (1, 6, 1), (1, 7, 2), (2, 5, 3), (3, 7, 4), (4, 7, 1), (7, 8, 1)])
 #   G.add_weighted_edges_from([(3,4,1), (1,4,1),(0,5,1),(3,5,1),(0,6,1),(3,6,2),(2,5,1),(5,6,2),(2,6,1), (1,7,1)]) # GRAPH FROM PAPER WITH SWAP COUNT 15
 #   gates = [(-1, 3, 4), (-1, 1, 4), (0, 3, 5), (0, 3, 6), (-1, 3, 6), (2, 6, 5), (2, 5, 6), (-1, 1, 7)]

    # rd_84 G.add_weighted_edges_from([(0,1,1),(0,8,1),(1,2,1),(1,8,2),(2,3,1),(2,8,2),(2,9,1),(3,4,1,),(3,8,2),(3,9,1),(3,10,1),(4,5,1),(4,8,2,),(4,9,1),(4,10,1),(4,11,1),(5,6,1),(5,8,2),(5,9,1),(5,10,1),(5,12,1),(6,7,1),(6,8,2),(6,9,1),(6,10,1),(6,13,1),(7,8,1),(7,10,1),(7,14,1),(8,9,5),(9,10,5),(10,11,1),(11,12,1),(12,13,1),(13,14,1)])
    #gates=[(0,1,8),(-1,0,1),(2,8,9),(1,2,8),(-1,1,2),(3,9,10),(3,8,9),(2,3,8),(-1,2,3),(4,10,11),(4,9,10),(4,8,9),(3,4,8),(-1,3,4),(5,11,12),(5,9,10),(5,8,9),(4,5,8),(-1,4,5),(6,12,13),(6,9,10),(6,8,9),(5,6,8),(-1,5,6),(7,13,14),(7,9,10),(6,7,8),(-1,6,7)]
 
    # mult_19_13
  #  G.add_weighted_edges_from([(3,8,1),(5,8,1),(2,8,1),(6,8,1),(1,8,1),(7,8,1),(3,9,1),(6,9,1),(2,9,1),(7,9,1),(3,10,1),(7,10,1),(10,11,1),(9,10,1),(8,9,1),(3,11,1),(4,11,1),(2,11,1),(5,11,1),(1,11,1),(6,11,1),(0,11,1),(7,11,1),(2,10,1),(4,10,1),(1,10,1),(5,10,1),(0,10,1),(6,10,1),(1,9,1),(4,9,1),(0,9,1),(5,9,1),(0,8,1),(4,8,1)])
  #  gates = [(3,5,8),(2,6,8),(1,7,8),(3,6,9),(2,7,9),(3,7,10),(-1,10,11),(-1,9,10),(-1,8,9),(3,4,11),(2,5,11),(1,6,11),(0,7,11),(2,4,10),(1,5,10),(0,6,10),(1,4,9),(0,5,9),(0,4,8)]
    
    # ham7
  #  G.add_weighted_edges_from([(4,5,2),(6,1,1),(6,2,2),(4,6,3),(5,6,2),(3,2,3),(3,0,4),(3,1,3),(5,3,3),(4,2,1),(4,3,1),(6,3,2),(1,4,1),(2,1,1),(2,5,1),(0,2,1)])
 #   gates = [(-1,4,5),(-1,6,1),(-1,6,2),(4,5,6),(-1,6,4),(-1,6,5),(-1,3,2),(-1,3,0),(-1,3,1),(2,5,3),(4,6,2),(2,5,3),(-1,3,1),(5,6,4),(-1,3,0),(-1,3,5),(4,6,3),(-1,1,4),(-1,2,1),(-1,3,6),(-1,2,5),(-1,1,3),(-1,3,0),(-1,0,3),(-1,0,2)]
    
        #dTof
    G.add_weighted_edges_from([(2,1,2),(1,0,2),(0,3,2),(1,3,1)])
    gates = [(-1,2,1),(-1,1,0),(-1,0,3),(-1,1,0),(-1,0,3),(-1,1,3),(-1,2,1)]

        #decod24
 #   G.add_weighted_edges_from([(2,1,3),(3,1,2),(3,0,1),(0,2,2)])
 #   gates = [(-1,2,1),(-1,3,1),(-1,3,0),(-1,0,2),(-1,2,1),(-1,1,2),(-1,2,0),(-1,1,3)] 

    vdegree = list(dict(G.degree(weight="weight")).values())
    nparts = nq
 
    N = 1
    res=[]
    swapCounts = []
    trials = [0]*N
    for i in range(N):
        total_ordering = [-1]
        while (total_ordering.count(-1) > 0 or len(total_ordering) != len(G.nodes)):
            total_ordering, tr = graphPart(G, nparts, vdegree, [], trials[i],sim)
            trials[i] = tr
        swapCount = calculateSWAP(total_ordering, gates)
        res.append({"total_ordering": total_ordering, "swapCount": swapCount})
        swapCounts.append(swapCount)

    uniqueSwapCounts = list(Counter(swapCounts).items())

    for i in range(len(res)):
        res[i]["swapFrequency"] = [item[1] for item in uniqueSwapCounts if item[0] == res[i]["swapCount"]][0]/N
    
    res_sorted = sorted(res, key= lambda x: x["swapFrequency"], reverse=True)

    for i in range(len(res_sorted)):
        print(*list(res_sorted[i].values()), sep=";")
        print("\n")
    
  #  sys.stdout = orig_stdout
 #   f.close()

    return 0


try:
    main()
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)
