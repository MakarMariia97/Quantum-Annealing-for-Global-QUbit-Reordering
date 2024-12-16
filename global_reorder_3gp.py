import networkx as nx
import math
import sys
from collections import Counter
import time
from utils import solve_QUBO, checkBalance, makePermutation, print_graph, calculateSWAP

# function to perform balanced partitioning of graph G into 3 parts recursively
def graphPart(G, nparts, vdegree, total_ordering, trials,sim=False):
    print_graph(G)
    k=3
    if (nparts>1):
        nNodes = len(G.nodes)
        is_ok = [-1]*nNodes
        trials -= 1
        unbal=-1
        while (sum(is_ok)!=0 and unbal!=0): # sample new solution until balanced one is appeared
            trials += 1
            sample = solve_QUBO(G, vdegree, k, sim)
            for i in range(nNodes):
                if (sample[i]+sample[nNodes+i]+sample[2*nNodes+i] == 1): # Check if each node is at only one part
                    is_ok[i] = 0
                else:
                    is_ok[i] = -1
            part3 = makePermutation(G, nNodes, k, sample)
            unbal = checkBalance(G, nNodes, part3, k) # Check if partition balanced
            
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
                
        # Check if new graphs cannot be partitioned again (not enough number of vertices)
        if (len(list(lG.nodes)) < k and len(list(mG.nodes)) < k and len(list(rG.nodes)) < k):
            return (list(lG.nodes) + list(mG.nodes) + list(rG.nodes),trials)

        # Partition new graphs
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

def main():
    start_time = time.time()
    
    nq=6 # number of qubits needed to optimize double Toffoli or 2-4 decoder (4 + 2 additional qubits)
    # nq=9 # number of qubits needed to optimize Hamming gate (7 + 2 additional qubits) or modified circuit from [31] (8 + 1 additional qubit)
    # nq=18 # number of qubits needed to optimize multiplier (12 + 6 additional qubits)

    sim=True # simulation mode, if True

    G = nx.Graph()
    G.add_nodes_from(range(nq))
    # modified circuit from [31]
    # G.add_weighted_edges_from([(3,4,1), (1,4,1),(0,5,1),(3,5,1),(0,6,1),(3,6,2),(2,5,1),(5,6,2),(2,6,1), (1,7,1)])
    # set circuit's gates in a form of (control1, control2, target), if there is only one control, then control1=-1 and control2=control    
    # gates = [(-1, 3, 4), (-1, 1, 4), (0, 3, 5), (0, 3, 6), (-1, 3, 6), (2, 6, 5), (2, 5, 6), (-1, 1, 7)]
    # real_nq=8 # number of circuit within the circuit originally
    
    # multiplier gate
    # G.add_weighted_edges_from([(3,8,1),(5,8,1),(2,8,1),(6,8,1),(1,8,1),(7,8,1),(3,9,1),(6,9,1),(2,9,1),(7,9,1),(3,10,1),(7,10,1),(10,11,1),(9,10,1),(8,9,1),(3,11,1),(4,11,1),(2,11,1),(5,11,1),(1,11,1),(6,11,1),(0,11,1),(7,11,1),(2,10,1),(4,10,1),(1,10,1),(5,10,1),(0,10,1),(6,10,1),(1,9,1),(4,9,1),(0,9,1),(5,9,1),(0,8,1),(4,8,1)])
    # gates = [(3,5,8),(2,6,8),(1,7,8),(3,6,9),(2,7,9),(3,7,10),(-1,10,11),(-1,9,10),(-1,8,9),(3,4,11),(2,5,11),(1,6,11),(0,7,11),(2,4,10),(1,5,10),(0,6,10),(1,4,9),(0,5,9),(0,4,8)]
    # real_nq = 12
    
    # Hamming gate 
    # G.add_weighted_edges_from([(4,5,2),(6,1,1),(6,2,2),(4,6,3),(5,6,2),(3,2,3),(3,0,4),(3,1,3),(5,3,3),(4,2,1),(4,3,1),(6,3,2),(1,4,1),(2,1,1),(2,5,1),(0,2,1)])
    # gates = [(-1,4,5),(-1,6,1),(-1,6,2),(4,5,6),(-1,6,4),(-1,6,5),(-1,3,2),(-1,3,0),(-1,3,1),(2,5,3),(4,6,2),(2,5,3),(-1,3,1),(5,6,4),(-1,3,0),(-1,3,5),(4,6,3),(-1,1,4),(-1,2,1),(-1,3,6),(-1,2,5),(-1,1,3),(-1,3,0),(-1,0,3),(-1,0,2)]
    # real_nq = 7
    
    # double Toffoli date
    G.add_weighted_edges_from([(2,1,2),(1,0,2),(0,3,2),(1,3,1)])
    gates = [(-1,2,1),(-1,1,0),(-1,0,3),(-1,1,0),(-1,0,3),(-1,1,3),(-1,2,1)]
    real_nq = 4
    
    # 2-4 decoder
    # G.add_weighted_edges_from([(2,1,3),(3,1,2),(3,0,1),(0,2,2)])
    # gates = [(-1,2,1),(-1,3,1),(-1,3,0),(-1,0,2),(-1,2,1),(-1,1,2),(-1,2,0),(-1,1,3)] 
    # real_nq = 4
    
    vdegree = list(dict(G.degree(weight="weight")).values())
    nparts = nq
 
    N = 1 # number of repetitions
    res=[]
    swapCounts = []
    trials = [0]*N # count unsuccessful trials during each repetition
    for i in range(N):
        total_ordering = [-1]
        while (total_ordering.count(-1) > 0 or len(total_ordering) != len(G.nodes)):
            total_ordering, tr = graphPart(G, nparts, vdegree, [], trials[i],sim)
            trials[i] = tr
        for i in range(real_nq,nq): # exclude additional qubits
            if (total_ordering.count(i) > 0):
                total_ordering.remove(i)
        swapCount = calculateSWAP(total_ordering, gates)
        res.append({"total_ordering": total_ordering, "swapCount": swapCount})
        swapCounts.append(swapCount)

    uniqueSwapCounts = list(Counter(swapCounts).items())

    for i in range(len(res)):
        res[i]["swapFrequency"] = [item[1] for item in uniqueSwapCounts if item[0] == res[i]["swapCount"]][0]/N
    
    res_sorted = sorted(res, key= lambda x: x["swapFrequency"], reverse=True)
    
    # ------- Print set of qubit permutations with the corresponding SWAP count and its frequency -------
    for i in range(len(res_sorted)):
        print(*list(res_sorted[i].values()), sep=";")
        print("\n")
        
    return 0

try:
    main()
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)
