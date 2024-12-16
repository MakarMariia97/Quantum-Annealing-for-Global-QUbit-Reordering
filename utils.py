from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math
from neal import SimulatedAnnealingSampler

# ------- Set tunable parameters -------
num_reads = 1000

# function to partition graph G into k parts using Quantum Annealing or Simulated Annealing
def solve_QUBO(G, vdegree, k, sim):
    # ------- Set up our QUBO dictionary -------
    # Initialize our Q matrix
    Q = defaultdict(int)

    nNodes = len(G.nodes)
 
    alpha_lagr = 4 # penalty constant
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

# function to check if partition is balanced 
def checkBalance(G, nNodes, part, k):
    dd=[0]*k
    unbal=-1
    for i in range(k):
        dd[i] = len([o for o, e in enumerate(part) if e == i])
        if (dd[i] in [math.floor(len(G.nodes)/k), math.ceil(len(G.nodes)/k)]):
            unbal=0
        else:
            unbal=-1
    return unbal

# function to convert binary sample to qubit permutation
def makePermutation(G, nNodes, k, sample):
    part = [-1]*nNodes
    for i in range(k):
        for u in range(len(G.nodes)):
            if sample[i*nNodes+u] == 1:
                part[u] = i
    return part

# function to print info about the graph G
def print_graph(G):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    print("nodes", G.nodes)
    return 0
    print("edges", G.edges)
