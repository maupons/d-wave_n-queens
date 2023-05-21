# 
# Grafo, 2 entradas 2 salidas - Ising Model
# 
#

from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
import neal
from time import time
import networkx as nx


def n_graph(n):
    """Returns a potential solution to the n-queens problem in a dictionary, where
    the keys are qubit id and value on or off

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. 
        Can use: 1) LeapHybridSampler
                 2) DWaveSampler - QPU
                 3) SimulatedAnnealingSampler
    """

    bqm = BinaryQuadraticModel({}, {}, 0, 'SPIN')
    bqm.offset = 0
    itr = 200
    w = 1
    # cs = 0.5

    # Build a graph with 5 nodes and 2 outputs each
    bqm.add_variable(1, 0)
    bqm.add_variable(2, 0)
    bqm.add_variable(3, 0)
    bqm.add_variable(4, 0)
    bqm.add_variable(5, 0)

    bqm.add_interaction(1, 2, w)
    bqm.add_interaction(1, 3, -1*w)
    bqm.add_interaction(2, 3, w)
    bqm.add_interaction(2, 4, w)
    bqm.add_interaction(3, 4, -1*w)
    bqm.add_interaction(3, 5, w)
    bqm.add_interaction(4, 5, w)
    bqm.add_interaction(4, 1, w)
    bqm.add_interaction(5, 1, -1*w)
    bqm.add_interaction(5, 2, w)

    print(bqm)

    print('solver started...')
    start_time = time()

    
    # QPU
    sampler = EmbeddingComposite(DWaveSampler())
    # sampler = EmbeddingComposite(DWaveSampler(solver=dict(topology__type='pegasus')))
    # sampler = DWaveSampler(solver=dict(topology__type='zephyr'))
    # embedded = EmbeddingComposite(sampler)
    # sampler = FixedEmbeddingComposite(DWaveSampler(),embedding)
    # sampler = FixedEmbeddingComposite(DWaveSampler(solver={'qpu': True}), embedding)
    # print(DWaveSampler().properties)
    # sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n} QD - QPU Fix Embbeded')    
    sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n} node Graph')

    # Hybrid Solver
    # sampler = LeapHybridSampler()                             # Hybrid Solver
    # sampleset = sampler.sample(bqm, label=f'{n}-q; {d}-d config')
    # print('quota_conversion_rate', sampler.properties['quota_conversion_rate'])

    # CPU
    # sampler = neal.SimulatedAnnealingSampler()                  # CPU
    # sampleset = sampler.sample(bqm, num_reads=itr)
    # print('category', sampler.properties['category'])
      
    solver_time = f'finished in {time()-start_time} seconds'
    print(solver_time)
    print('sampleset.info', sampleset.info)

    f = open("sampleset.txt", "w")
    for sample in sampleset:
        f.write(str(sample)+'\n')
    f.close()

    f = open("sampleset2.txt", "w")
    f.write(str(sampleset))
    f.close()

    sample = sampleset.first.sample

    H = 0
    H += -1*sample[1]*sample[2]
    H += -1*sample[1]*sample[3]
    H += -1*sample[2]*sample[3]
    H += -1*sample[2]*sample[4]
    H += -1*sample[3]*sample[4]
    H += -1*sample[3]*sample[5]
    H += -1*sample[4]*sample[5]
    H += -1*sample[4]*sample[1]
    H += -1*sample[5]*sample[1]
    H += -1*sample[5]*sample[2]

    print("Sample:\n", sample)
    f = open("outGraph.txt", "a")
    f.write(f'{n} node graph\n')
    f.write(str(sample)+'\n')
    f.write(f'Energy is {H}'+'\n')
    f.write(solver_time + '\n')
    f.write(str(sampleset.info)+'\n\n')
    f.close()

    
    # Inspect
    dwave.inspector.show(sampleset)
    return sample


def is_valid_solution(n, solution):
    """Check energy

    Args:
        n: Number of queens in the problem.

        solution: A dictionary of qubits, key = qubit id, value = 0 or 1
    """

    H = 0
    H += 1*solution[1]*solution[2]
    H += 1*solution[1]*solution[3]
    H += 1*solution[2]*solution[3]
    H += 1*solution[2]*solution[4]
    H += 1*solution[3]*solution[4]
    H += 1*solution[3]*solution[5]
    H += 1*solution[4]*solution[5]
    H += 1*solution[4]*solution[1]
    H += 1*solution[5]*solution[1]
    H += 1*solution[5]*solution[2]

    print('Energy found ', H)

    # Calculate Energy with H = - sum sum J_ik s_i s_k
    H = 0
    H += -1*solution[1]*solution[2]
    H += -1*solution[1]*solution[3]
    H += -1*solution[2]*solution[3]
    H += -1*solution[2]*solution[4]
    H += -1*solution[3]*solution[4]
    H += -1*solution[3]*solution[5]
    H += -1*solution[4]*solution[5]
    H += -1*solution[4]*solution[1]
    H += -1*solution[5]*solution[1]
    H += -1*solution[5]*solution[2]

    print('Energy Hans found ', H)

    return True, H


def plot_chessboard(n, queens):
    """Create a chessboard with queens using matplotlib. Image is saved
    in the root directory. Returns the image file name.
    """
        # Build graph for visualization
    G = nx.Graph()

    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    G.add_node(4)
    G.add_node(5)

    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(2, 3)
    G.add_edge(2, 4)
    G.add_edge(3, 4)
    G.add_edge(3, 5)
    G.add_edge(4, 5)
    G.add_edge(4, 1)
    G.add_edge(5, 1)
    G.add_edge(5, 2)

    pos = {
        1:(3,1),
        2:(4,3),
        3:(2,4),
        4:(0,3),
        5:(1,1),
    }
    
    nx.draw(G,pos=pos,with_labels=True,node_color="red",node_size=3000,
            font_color="white",font_size=20,font_family="Times New Roman",
            font_weight="bold",width=5)
    plt.margins(0.2)
    # plt.show()

    file_name = 'grafo5.png'
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    n = 5

    print("Building graph with {n} nodes.".format(n=n))
    solution = n_graph(n)

    if is_valid_solution(n, solution):
        write = "Solution is valid."
    else:
        write = "Solution is invalid."

    file_name = plot_chessboard(n, solution)
    print("Graph created. See: {}".format(file_name))
