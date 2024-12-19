from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math
from neal import SimulatedAnnealingSampler

# ------- Set tunable parameters -------
num_reads = 1000

def solve_QUBO(G, vdegree, k, sim):
    # ------- Set up our QUBO dictionary -------
    # Initialize our Q matrix
    Q = defaultdict(int)

    nNodes = len(G.nodes)
 
    alpha_lagr = 4
    for i in range(k):
        for u in range(len(G.nodes)):
            Q[(i*nNodes+u,i*nNodes+u)] += vdegree[u] + (1-2*nNodes/k)*alpha_lagr - alpha_lagr
        for v, w in combinations(G.nodes, 2):
            Q[(i*nNodes+list(G.nodes).index(v), i*nNodes+list(G.nodes).index(w))] += 2*alpha_lagr
        for v, w, d in G.edges(data=True):
            Q[(i*nNodes+list(G.nodes).index(v), i*nNodes+list(G.nodes).index(w))] += -2*d["weight"]
    for i, j in combinations(range(k), 2):
        for u in range(len(G.nodes)):
            Q[(i*nNodes+u, j*nNodes+u)] += 2*alpha_lagr

    # Run the QUBO on the solver from your config file

    if sim:
        sampler = SimulatedAnnealingSampler()
        response = sampler.sample_qubo(Q,
                                   num_reads=num_reads,
                                   label='Example - Graph Partitioning')
    else:
        sampler = EmbeddingComposite(DWaveSampler())
        response = sampler.sample_qubo(Q,
                                    num_reads=num_reads,
                                   label='Example - Graph Partitioning',return_embedding=True)
        embed_info =  response.info['embedding_context']['embedding']
        print(f"Number of logical variables: {len(embed_info.keys())}")
        print(f"Number of physical qubits used in embedding: {sum(len(chain) for chain in embed_info.values())}")
    sample = response.record.sample[0]

    return sample

def checkBalance(G, part, k):
    dd=[0]*k
    unbal=-1
    for i in range(k):
        dd[i] = len([o for o, e in enumerate(part) if e == i])
        if (dd[i] in [math.floor(len(G.nodes)/k), math.ceil(len(G.nodes)/k)]):
            unbal=0
        else:
            unbal=-1
    return unbal

def makePermutation(G, nNodes, k, sample):
    part = [-1]*nNodes
    for i in range(k):
        for u in range(len(G.nodes)):
            if sample[i*nNodes+u] == 1:
                part[u] = i
    return part

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

def print_graph(G):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))

def select_circuit(pr,pr1):
    if pr == "1":
        cnotbased = True
        dToffoli = False
        multipler = False
        hamming = False
        decod = False
    elif pr == "2":
        cnotbased = False
        dToffoli = True
        multipler = False
        hamming = False
        decod = False
    elif pr == "3":
        cnotbased = False
        dToffoli = False
        multipler = True
        hamming = False
        decod = False
    elif pr == "4":
        cnotbased = False
        dToffoli = False
        multipler = False
        hamming = True
        decod = False
    elif pr == "5":
        cnotbased = False
        dToffoli = False
        multipler = False
        hamming = False
        decod = True
    else:
        print("[ERROR] string is not valid, exiting...")
        exit(2)
    
    if pr1 == "1":
        sim=False
    elif pr1 == "0":
        sim=True
    else:
        print("[ERROR] string is not valid, exiting...")
        exit(2)
    return cnotbased, dToffoli, multipler, hamming, decod, sim
