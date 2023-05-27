# Copyright 2020 D-Wave Systems Inc.
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


# Descripcion:
# 
# Stats for embedding time
# 
# Ejemplo de find_embedding() en:
#  https://docs.dwavesys.com/docs/latest/c_qpu_timing.html
#  sección: Estimating the QPU Access Time for Problems


import numpy as np
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from minorminer import busclique
from minorminer import find_embedding
import dwave.inspector
import neal
import sys
from time import time


def n_queens(n,n_emb):

    dp = []
    dm = []
    d = len(dp) + len(dm)

    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    w = 2

    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
        if (ri+ci) in dp:
            bqm.add_variable(i, w)
        elif (ri-ci) in dm:
            bqm.add_variable(i, w)
        else:
            bqm.add_variable(i, -w)
            
        for j in range(i):
            rj = j // n
            cj = j-n*rj
            if rj == ri:
                bqm.add_interaction(i, j, w)
            if cj == ci:
                bqm.add_interaction(i, j, w)
            if abs(ri-rj) == abs(ci-cj):
                bqm.add_interaction(i, j, w)

    # print(bqm)
    # print(list(bqm.quadratic.keys()))
    print(f'Embedding {n_emb} times...')

    f = open(f"data/embedding/find/chimera/time_n{n}_r{n_emb}.txt", "a")
    f.write(f'n\tqpu_t\tfind_embedding_t\tFixedEmbedding_t\tnum_qubits\ttopology\n')

    for e in range(n_emb):
        tick = time()
        qpu = DWaveSampler(solver={'topology__type': 'chimera'})
        # qpu = DWaveSampler(solver={'topology__type': 'pegasus'})
        # qpu = DWaveSampler(solver={'topology__type': 'zephyr'})
        qpu_t = time() - tick
        
        tick = time()
        embedding = find_embedding(list(bqm.quadratic.keys()), qpu.to_networkx_graph())
        find_embedding_t = time() - tick
        
        tick = time()
        sampler = FixedEmbeddingComposite(qpu, embedding)
        FixedEmbedding_t = time() - tick

        num_qubits = sum(len(chain) for chain in embedding.values())
        topology = qpu.properties['topology']['type']
        print(embedding)
        print('num_qubits', num_qubits)

        f.write(f'{n}\t{qpu_t}\t{find_embedding_t}\t{FixedEmbedding_t}\t{num_qubits}\t{topology}\n')
        f.flush()
    f.close()

    # print(qpu)
    # print(qpu.properties['topology'])
    # print(qpu.parameters)
    # print(qpu.to_networkx_graph())    
    # print('embedding.values', embedding.values())

    # print('embedding', embedding)


if __name__ == "__main__":

    tick = time()
    n = int(sys.argv[1])
    # n = 5
    n_emb = 1

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    n_queens(n,n_emb)

    tock = time() - tick
    print('Exe(s)', tock)

