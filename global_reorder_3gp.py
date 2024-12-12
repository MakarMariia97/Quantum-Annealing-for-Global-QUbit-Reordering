import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math
from neal import SimulatedAnnealingSampler
import sys
import dwave_networkx as dnx
import dimod
from collections import Counter
import time
from scipy import stats
import statistics 

def solve_QUBO(G, nparts, vdegree):
    #print("vdegree ", vdegree)
    # ------- Set up our QUBO dictionary -------
    print("nparts ", nparts)
    # Initialize our Q matrix
    Q = defaultdict(int)

    nNodes = len(G.nodes)
 
    alpha_lagr = 4

    for i in range(4):
    #    print("i ",i)
        for u in range(len(G.nodes)):
    #        print("u ", u)
            Q[(i*nNodes+u,i*nNodes+u)] += vdegree[u] + (1-2*nNodes/4)*alpha_lagr - alpha_lagr
        for v, w in combinations(G.nodes, 2):
            Q[(i*nNodes+list(G.nodes).index(v), i*nNodes+list(G.nodes).index(w))] += 2*alpha_lagr
        for v, w, d in G.edges(data=True):
            Q[(i*nNodes+list(G.nodes).index(v), i*nNodes+list(G.nodes).index(w))] += -2*d["weight"]
    for i, j in combinations(range(4), 2):
        for u in range(len(G.nodes)):
            Q[(i*nNodes+u, j*nNodes+u)] += 2*alpha_lagr
    # Run the QUBO on the solver from your config file

    #print("Q ", Q)

  #  sampler = EmbeddingComposite(DWaveSampler())

    sampler = SimulatedAnnealingSampler()
    response = sampler.sample_qubo(Q,
                                   #chain_strength=100,
                                   num_reads=num_reads,
                                   label='Example - Graph Partitioning')#,return_embedding=True)

    sample = response.record.sample[0]
    print("length of sample ", len(sample))
  #  embed_info =  response.info['embedding_context']['embedding']
  #  print(f"Number of logical variables: {len(embed_info.keys())}")
  #  print(f"Number of physical qubits used in embedding: {sum(len(chain) for chain in embed_info.values())}")
    #print(sample)
    return sample


def graphPart(G, nparts, vdegree, total_ordering, trials):
    print_graph(G)
    print("nparts ",nparts)

    if (nparts>1):
        nNodes = len(G.nodes)
        is_ok = [-1]*nNodes
        indices_ones = [-1]
        duples_ones=[-1]
        indices_col = [-1]
        duples_col=[-1]
        trials -= 1
        unbal=-1
        while (sum(is_ok)!=0 and unbal!=0):
            trials += 1
            sample = solve_QUBO(G, nparts, vdegree)
            #print("sample in while ", sample)
            #print("sample ", sample)
            for i in range(nNodes):
                if (sample[i]+sample[nNodes+i]+sample[2*nNodes+i] == 1):
                    is_ok[i] = 0
                else:
                    is_ok[i] = -1
            part3 = [-1]*nNodes
            # print("part3 ", part3)
            for i in range(3):
                for u in range(len(G.nodes)):
                    if sample[i*nNodes+u] == 1:
        #            print("u ", u)
                        part3[u] = i
            print("part3 ", part3)

            dd=[0]*3
            for i in range(3):
          #      print(i, "partition")
                dd[i] = len([o for o, e in enumerate(part3) if e == i])
          #      print(dd[i], "elements")
                if (dd[i] in [math.floor(len(G.nodes)/3), math.ceil(len(G.nodes)/3)]):
                  #  print(i, "balanced")
                    unbal=0
                else:
               #     print(i, "unbalanced")
                    unbal=-1
        # partition graph into 3 parts
        lG = nx.Graph()
        mG = nx.Graph()
        rG = nx.Graph()
        print("balanced:", part3)


        if (sum(part3) % 2 != 0):
           # print("HERE!!!")
            if (part3.count(-1)>0):
                part3[part3.index(-1)] = nNodes - sum(part3)
            duples = [x for n, x in enumerate(part3) if x in part3[:n]]
            if (duples != []):
                part3[part3.index(duples[0])] = abs(nNodes - sum(part3))
       # print("part3 after removing duplicates ", part3)
        for u in range(len(G.nodes)):
            if part3[u] == 0:
                lG.add_node(list(G.nodes)[u])
            if part3[u] == 1:
                mG.add_node(list(G.nodes)[u])
            if part3[u] == 2:
                rG.add_node(list(G.nodes)[u])

        for u, v, d in G.edges(data=True):
            #print("u, v", u, v)
            if part3[list(G.nodes).index(u)] == 0 and part3[list(G.nodes).index(v)] == 0:
                lG.add_edge(u, v, weight = d["weight"])
            if part3[list(G.nodes).index(u)] == 1 and part3[list(G.nodes).index(v)] == 1:
                mG.add_edge(u, v,  weight = d["weight"])
            if part3[list(G.nodes).index(u)] == 2 and part3[list(G.nodes).index(v)] == 2:
                rG.add_edge(u, v,  weight = d["weight"])



        if (len(list(lG.nodes)) <= 3 and len(list(mG.nodes)) <= 3 and len(list(rG.nodes)) <= 3):
      #      print("partial ord ", list(llG.nodes) + list(lG.nodes) + list(rG.nodes) + list(rrG.nodes))
            return (list(lG.nodes) + list(mG.nodes) + list(rG.nodes),trials)

        if (len(list(lG.nodes())) >= 3):
            total_ord_from_part_lG, trial_lG = graphPart(
                lG, math.floor(nparts/3), list(dict(lG.degree(weight="weight")).values()), total_ordering, trials)
        #    print("partial ordering AFTER PARTITIONING llG ",
        #          total_ord_from_part_llG)
        #    print("total ordering before ", total_ordering)
            total_ordering = total_ordering + total_ord_from_part_lG
            trials = trial_lG
         #   print("total ordering after ", total_ordering)
        if (len(list(mG.nodes())) >= 3):
            total_ord_from_part_mG, trial_mG = graphPart(
                mG, math.floor(nparts/3), list(dict(mG.degree(weight="weight")).values()), total_ordering, trials)
         #   print("partial ordering AFTER PARTITIONING lG ",
        #          total_ord_from_part_lG)
        #    print("total ordering before ", total_ordering)
            total_ordering = total_ordering + total_ord_from_part_mG
            trials = trial_mG
        #    print("total ordering after ", total_ordering)
        if (len(rG.nodes) >= 3):
            total_ord_from_part_rG, trial_rG = graphPart(
                rG, math.floor(nparts/3), list(dict(rG.degree(weight="weight")).values()), [], trials)
       #     print("partial ordering AFTER PARTITIONING rG ",
       #           total_ord_from_part_rG)
        #    print("total ordering before ", total_ordering)
            total_ordering = total_ordering + total_ord_from_part_rG
            trials = trial_rG
        #    print("total ordering after ", total_ordering)

        #    print("total ordering after ", total_ordering)

    #print(total_ordering)
    return (total_ordering, trials)

def print_graph(G):
 #   print("Graph on {} nodes created with {} out of {} possible edges.".format(
 #       len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    print("nodes", G.nodes)
 #   print("edges", G.edges)

# ------- Set tunable parameters -------
num_reads = 1000
#gamma = 80

# ------- Set up our graph -------


def calculateSWAP(total_ordering, gates):
    sum=0
    for u, v, w in gates:
    #    print ("u v w ", u, v, w)
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
                    # modified calculation
                    sum+= abs(ind_u - ind_w) - 1 + abs(ind_v - ind_u) - 1
                elif (ind_v - ind_u>1):
                        sum+=ind_v - ind_u - 1
       #         print("sum ", sum)
            elif (ind_u < ind_w & ind_w < ind_v):
                if (ind_w-ind_u>1):
                    sum+= abs(ind_w - ind_u) - 1
                if (ind_v - ind_w>1):
                        sum+=ind_v - ind_w - 1
         #       print("sum ", sum)
            elif (ind_u < ind_v & ind_v < ind_w):
                if (ind_w - ind_v >1):
                    sum+= abs(ind_w - ind_v) - 1 + abs(ind_v - ind_u) - 1
                elif (ind_v - ind_u>1):
                    sum+= abs(ind_v - ind_u) - 1
          #      print("sum ", sum)
        else:
            sum+=abs(ind_v-ind_w) - 1
          #  print("sum ", sum)
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
    
    # RELABEL the graph to count the new distance
    #total_ordering = graphPart(G, nparts, vdegree, [])
    #print(total_ordering)
  #  mapping={}
  #  for i in range(len(total_ordering)):
  #      mapping[total_ordering[i]]=i

  #  print("total order")

  #  for i in range(nparts):
   #     print("part ", i)
   #     for j in range(nq):
   #         print(total_ordering[i*nq+j], end=" ")
   #     print("\n")
    
  #  Gnew=nx.relabel_nodes(G,mapping)
    
  #  print_graph(Gnew)
  #  distance(Gnew)
    
    #print("total ordering", graphPart(G, nparts, []))
    N = 1
    res=[]
    swapCounts = []
    trials = [0]*N
    for i in range(N):
        #print("iteration i ", i)
        total_ordering = [-1]
        while (total_ordering.count(-1) > 0 or len(total_ordering) != len(G.nodes)):
            total_ordering, tr = graphPart(G, nparts, vdegree, [], trials[i])
            trials[i] = tr
            #print(total_ordering)
        if (total_ordering.count(12) > 0):
            total_ordering.remove(12)
        if (total_ordering.count(13) > 0):
            total_ordering.remove(13)
        if (total_ordering.count(14) > 0):
            total_ordering.remove(14)
        if (total_ordering.count(15) > 0):
            total_ordering.remove(15)
        #if (total_ordering.count(7) > 0):
        #    total_ordering.remove(7)
            #print(total_ordering)
        swapCount = calculateSWAP(total_ordering, gates)
        res.append({"total_ordering": total_ordering, "swapCount": swapCount})
        swapCounts.append(swapCount)
    #res_sorted = sorted(res, key= lambda x: x["swapCount"])

    uniqueSwapCounts = list(Counter(swapCounts).items())
    #swapSum=Counter(swapCounts).total()

    for i in range(len(res)):
        res[i]["swapFrequency"] = [item[1] for item in uniqueSwapCounts if item[0] == res[i]["swapCount"]][0]/N
    
    res_sorted = sorted(res, key= lambda x: x["swapFrequency"], reverse=True)

    for i in range(len(res_sorted)):
        print(*list(res_sorted[i].values()), sep=";")
        print("\n")

    skew = stats.skew(swapCounts)
    mean = statistics.mean(swapCounts)
    mode = statistics.mode(swapCounts)
    median = statistics.median(swapCounts)

    print("skewness ", skew)
    print("mean ", mean)
    print("mode ", mode)
    print("median ", median)
    print("trials ", sum(trials))

    end_time = time.time()
    print("time ", end_time - start_time)
    
  #  sys.stdout = orig_stdout
 #   f.close()

    return 0


try:
    main()
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)