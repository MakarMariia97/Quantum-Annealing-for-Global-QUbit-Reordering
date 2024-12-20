# Copyright 2019 D-Wave Systems, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ------ Import necessary packages ----
import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math
from neal import SimulatedAnnealingSampler
import sys
from collections import Counter
from utils import calculateSWAP, select_circuit

def solve_QUBO(G, trials, sim):
    # ------- Set up our QUBO dictionary -------
    # Initialize our Q matrix
    Q = defaultdict(int)

    # Fill in Q matrix
    for u, v, d in G.edges(data=True):
        Q[(u, u)] += d["weight"]
        Q[(v, v)] += d["weight"]
        Q[(u, v)] += -2*d["weight"]

    g_par = 4 # penalty constant
    for i in G.nodes():
        Q[(i, i)] += g_par*(1-len(G.nodes))
    for i, j in combinations(G.nodes, 2):
        Q[(i, j)] += 2*g_par

    # ------- Run our QUBO on the QPU -------

    # Run the QUBO on the solver

    sample = []
    sample.append(-1)
    trials -= 1
    while (sample[0] == -1):
        trials+=1
        if sim: # use Simulated Annealing in simulation mode
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

        # See if the best solution found is feasible, and if so print the number of cut edges.
        sample = response.record.sample[0]
        
        if sum(sample) in [math.floor(len(G.nodes)/2), math.ceil(len(G.nodes)/2)]:
            print("Valid partition found.")
            sample = sample
        else:
            print("Invalid partition.")
            sample = [-1]
    print(sample)
    return (sample, trials)


def graphPart(G, nparts, total_ordering, trials, sim=False):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    
    if (nparts > 1):
        sample, trial = solve_QUBO(G, trials, sim)

        # In the case when n is odd, the set may have one more or one fewer nodes
        if not(sum(sample) in [math.floor(len(G.nodes)/2), math.ceil(len(G.nodes)/2)]):
            return [-1]

        # partition graph into 2 parts
        lG = nx.Graph()
        rG = nx.Graph()
        for u in G.nodes:
            if sample[list(G.nodes).index(u)] == 0:
                lG.add_node(u)
            if sample[list(G.nodes).index(u)] == 1:
                rG.add_node(u)
        for u, v, d in G.edges(data=True):
            if sample[list(G.nodes).index(u)] == 0 and sample[list(G.nodes).index(v)] == 0:
                lG.add_edge(u, v, weight = d["weight"])
            if sample[list(G.nodes).index(u)] == 1 and sample[list(G.nodes).index(v)] == 1:
                rG.add_edge(u, v,  weight = d["weight"])
                
        # Check if new graphs cannot be partitioned balanced again
        if (len(list(lG.nodes)) == 1 & len(list(rG.nodes)) == 1):
            return (list(lG.nodes) + list(rG.nodes), trial)
            
        # Partition new graphs
        if (len(list(lG.nodes())) >= 2):
            total_ord_from_part_lG, trial_lG = graphPart(
                lG, math.floor(nparts/2), total_ordering, trial,sim)
            total_ordering = total_ordering + total_ord_from_part_lG
            trials = trial_lG
        if (len(list(rG.nodes())) >= 2):
            total_ord_from_part_rG, trial_rG = graphPart(rG, math.floor(nparts/2), [], trials,sim)
            total_ordering = total_ordering + total_ord_from_part_rG

    return (total_ordering, trials)

# ------- Set tunable parameters -------
num_reads = 1000

def main():
    # Create graph with nq vertices and edges from the set of circuit's gates
    G = nx.Graph()

    pr = input("Which quantum circuit would you like to optimize (CNOT-based, double Toffoli, mutilplier gate, Hamming or 2-4 decoder)? Press the corresponding index (starting with 1).")
    pr1 = input("Which mode would you like to run: Quantum Annealing or Simulated Annealing? Press 1 for Quantum Annealing and 0 otherwise.")
    cnotbased, dToffoli, multipler, hamming, decod, sim = select_circuit(pr,pr1)

    if dToffoli:
        # double Toffoli gate
        real_nq = 4 # original number of qubits within the circuit
        nq = 4 # total number of qubits (including additional ones)
        # set circuit's gates in a form of (control1, control2, target), if there is only one control, then control1=-1 and control2=control
        gates = [(-1,2,1),(-1,1,0),(-1,0,3),(-1,1,0),(-1,0,3),(-1,1,3),(-1,2,1)]
        for i in range(nq):
            G.add_node(i)
        G.add_weighted_edges_from([(2,1,2),(1,0,2),(0,3,2),(1,3,1)]) # weighted edge between control and target
    elif cnotbased:
        # modified circuit [31]
        real_nq = 8
        nq = 8
        gates = [(-1, 3, 4), (-1, 1, 4), (0, 3, 5), (0, 3, 6), (-1, 3, 6), (2, 6, 5), (2, 5, 6), (-1, 1, 7)]
        for i in range(nq):
            G.add_node(i)
        G.add_weighted_edges_from([(3,4,1), (1,4,1),(0,5,1),(3,5,1),(0,6,1),(3,6,2),(2,5,1),(5,6,2),(2,6,1), (1, 7, 1)])
    elif multipler:
        # multiplier
        real_nq = 12
        nq = 16 # +4 additional qubits
        gates = [(3,5,8),(2,6,8),(1,7,8),(3,6,9),(2,7,9),(3,7,10),(-1,10,11),(-1,9,10),(-1,8,9),(3,4,11),(2,5,11),(1,6,11),(0,7,11),(2,4,10),(1,5,10),(0,6,10),(1,4,9),(0,5,9),(0,4,8)]
        for i in range(nq):
            G.add_node(i)
        G.add_weighted_edges_from([(3,8,1),(5,8,1),(2,8,1),(6,8,1),(1,8,1),(7,8,1),(3,9,1),(6,9,1),(2,9,1),(7,9,1),(3,10,1),(7,10,1),(10,11,1),(9,10,1),(8,9,1),(3,11,1),(4,11,1),(2,11,1),(5,11,1),(1,11,1),(6,11,1),(0,11,1),(7,11,1),(2,10,1),(4,10,1),(1,10,1),(5,10,1),(0,10,1),(6,10,1),(1,9,1),(4,9,1),(0,9,1),(5,9,1),(0,8,1),(4,8,1)])        
    elif hamming:
        # Hamming gate
        real_nq = 7 
        nq = 8 # +1 additional qubit
        gates = [(-1,4,5),(-1,6,1),(-1,6,2),(4,5,6),(-1,6,4),(-1,6,5),(-1,3,2),(-1,3,0),(-1,3,1),(2,5,3),(4,6,2),(2,5,3),(-1,3,1),(5,6,4),(-1,3,0),(-1,3,5),(4,6,3),(-1,1,4),(-1,2,1),(-1,3,6),(-1,2,5),(-1,1,3),(-1,3,0),(-1,0,3),(-1,0,2)]
        for i in range(nq):
            G.add_node(i)
        G.add_weighted_edges_from([(4,5,2),(6,1,1),(6,2,2),(4,6,3),(5,6,2),(3,2,3),(3,0,4),(3,1,3),(5,3,3),(4,2,1),(4,3,1),(6,3,2),(1,4,1),(2,1,1),(2,5,1),(0,2,1)])       
    else:
        # 2-4 decoder
        real_nq = 4
        nq = 4
        gates = [(-1,2,1),(-1,3,1),(-1,3,0),(-1,0,2),(-1,2,1),(-1,1,2),(-1,2,0),(-1,1,3)] 
        for i in range(nq):
            G.add_node(i)
        G.add_weighted_edges_from([(2,1,3),(3,1,2),(3,0,1),(0,2,2)])
    
    nparts = nq # number of parts to be made

    N = 1 # number of repetitions of the algorithm to study the sensitivity to penalty costant
    res=[]
    swapCounts = []
    trials=[0]*N # number of unsuccessful trials during i-th run
    for i in range(N):
        total_ordering, tr = graphPart(G, nparts, [], trials[i],sim)
        trials[i] = tr
        while total_ordering.count(-1) > 0: # exclude invalid partition
            total_ordering, tr = graphPart(G, nparts, [], trials[i],sim)
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
