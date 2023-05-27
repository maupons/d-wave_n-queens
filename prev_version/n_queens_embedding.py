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
    emb_time_list = []
    for e in range(n_emb):
        start_time = time()
        sampler = EmbeddingComposite(DWaveSampler())
        embedding_time = time()-start_time
        emb_time_list.append(embedding_time)
    


    f4 = open("data/stats/embedding.txt", "a")
    line = f'{n}   {d}   '
    for t in emb_time_list:
        f4.write(line + f'{t}\n')
    f4.close()


if __name__ == "__main__":

    n = 7
    n_emb = 10

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    n_queens(n,n_emb)

