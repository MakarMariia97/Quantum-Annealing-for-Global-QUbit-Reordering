This repository implements optimization of quantum circuits by minimizing swap count. It supplements PhD thesis titled as "Optimizing Quantum Circuit
Layout using Quantum Annealing" written by Makarova Mariia (PhD student, University of Trento, DISI). Quantum circuit optimization is performed by global reordering of qubits by constructing balanced partitioning of the corresponding graph into $k$ parts using Quantum Annealng (on D-Wave) ($k=2;3;4$). Also, simulation mode is available. To test, at each file implementing $k$-partitioning, there are several quantum circuits (Hamming gate, multiplier, 2-4 decoder, double Toffoli gate) initialized as a qubit adjacency graph $G$ and set of gates $gates$ within the circuit.
In order to perform balanced partitioning, some additional qubits are added. Thus, variable $nq$ denotes the total number of qubits (including additional ones) and $real\_nq$ denotes original number of qubits within the circuit.

The approach is based on [Chakrabarti, Amlan & Sur Kolay, Susmita & Chaudhury, Ayan. (2011). Linear Nearest Neighbor Synthesis of Reversible Circuits by Graph
Partitioning].
