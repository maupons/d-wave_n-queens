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


#
# 
# Stats for embedding time
# 



import numpy as np
from dimod import BinaryQuadraticModel
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
import neal
from time import time


def n_queens(n,n_emb):

    dp = []
    dm = []
    d = len(dp) + len(dm)

    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    # bqm.offset = 2*n
    itr = 500
    w = 2
    # cs = 89

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
    
    print(f'Embedding {n_emb} times...')
    avg_emb_time = 0
    for e in range(n_emb):
        start_time = time()
        sampler = EmbeddingComposite(DWaveSampler())
        embedding_time = time()-start_time
        avg_emb_time += embedding_time
    print('sum =', avg_emb_time)
    avg_emb_time = avg_emb_time / n_emb
    print('avg =', avg_emb_time)

    f4 = open("data/stats/embedding.txt", "a")
    line = f'{n}   {avg_emb_time}\n'
    f4.write(line)
    f4.close()


if __name__ == "__main__":

    tick = time()
    n = 4
    n_emb = 500

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    n_queens(n,n_emb)

    tock = time() - tick
    print('Exe(s)', tock)

