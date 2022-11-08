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
# N Queen Solution using cell's row & columns
# 
# - block diagonals
# - solution verifier
# - draw board

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


def n_queens(n,dp,dm, sampler=None):
    """Returns a potential solution to the n-queens problem in a dictionary, where
    the keys are qubit id and value on or off

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. 
        Can use: 1) LeapHybridSampler
                 2) DWaveSampler - QPU
                 3) SimulatedAnnealingSampler
    """

    d = len(dp) + len(dm)

    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    bqm.offset = 2*n
    itr = 1000
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
    # print(bqm)
    

    print('solver started...')
    start_time = time()

    
    # QPU
    # sampler = EmbeddingComposite(DWaveSampler())
    # sampler = EmbeddingComposite(DWaveSampler(solver=dict(topology__type='pegasus')))
    # sampler = DWaveSampler(solver=dict(topology__type='zephyr'))
    # embedded = EmbeddingComposite(sampler)
    # sampler = FixedEmbeddingComposite(DWaveSampler(),embedding)
    # sampler = FixedEmbeddingComposite(DWaveSampler(solver={'qpu': True}), embedding)
    # print(DWaveSampler().properties)
    # sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n} QD - QPU Fix Embbeded')

    # tp = sampler.properties["topology"]['type']
    # sampleset = embedded.sample(bqm, num_reads=itr, chain_strength=cs, label=f'{n} QD, w:{w}, cs:{cs}, {tp}')
    # sampleset = embedded.sample(bqm, num_reads=itr, label=f'{n}Q, w:{w}, {tp}')

    # Hybrid Solver
    sampler = LeapHybridSampler()                             # Hybrid Solver
    sampleset = sampler.sample(bqm, label=f'{n}-q; {d}-d config')
    # print('quota_conversion_rate', sampler.properties['quota_conversion_rate'])

    # CPU
    # sampler = neal.SimulatedAnnealingSampler()                  # CPU
    # sampleset = sampler.sample(bqm, num_reads=itr, label=f'{n}-q; {d}-d config Classical')
    # print('category', sampler.properties['category'])
      
    solver_time = f'finished in {time()-start_time} seconds'
    print(solver_time)
    print('sampleset.info', sampleset.info)
    # print('sampler.properties', sampler.properties)
    # sampleset_iterator = sampleset.samples(num_iter)
    # f = open("sampleset.txt", "w")
    # for sample in sampleset:
    #     # print(sample)
    #     f.write(str(sample)+'\n')
    # f.close()

    sample = sampleset.first.sample
    print("Sample:\n", sample)
    f = open("outHybrid.txt", "a")
    f.write(f'{n}-q; {d}-d config\n')
    f.write(str(sample)+'\n')
    f.write(solver_time + '\n')
    f.write(str(sampleset.info)+'\n')
    f.close()

    
    # Inspect
    # dwave.inspector.show(sampleset)
    return sample,dp,dm


def is_valid_solution(n, solution):
    """Check that solution is valid by making sure all constraints were met.

    Args:
        n: Number of queens in the problem.

        solution: A dictionary of qubits, key = qubit id, value = 0 or 1
    """
    board = np.zeros((n,n))
    ldp = []                            # Keep track of queens on D+
    ldm = []                            # Keep track of queens on D-

    # Build board & Check diag/anti-diag constraints
    for qb,v in solution.items():
        r = qb // n
        c = qb-n*r
        if v:
            board[r,c] = v
            dp = r+c
            dm = r-c
            # ldp.append(dp)
            # ldm.append(dm)
            if dp not in ldp:
                ldp.append(dp)
            else:
                print(f"D+ diagonal {dp} has more than 1 queen")
                return False
            if dm not in ldm:
                ldm.append(dm)
            else:
                print(f"D- diagonal {dm} has more than 1 queen")
                return False

    # print(board)

    # Check row/col constraints
    for i in range(n):
        sum_row = sum(board[i,:])
        if sum_row != 1:
            print(f"Row {i} has {sum_row} queens.")
            return False
        sum_col = sum(board[:,i])
        if sum_col != 1:
            print(f"Column {i} has {sum_col} queens.")
            return False
    
    # I think this doesnt work becuase I cant tell if all n queens are in same row or col
    # print(sum(board))
    # print(sum(sum(board)))
    # return False if sum(sum(board)) != n else True

    # Check diag/anti-diag constraints
    # print(ldp)
    # print(ldm)
    # if len(ldp) != len(set(ldp)):
    #     print(f"A diagonal in D+ set has more than 1 queen.")
    #     return False
    # elif len(ldm) != len(set(ldm)):
    #     print(f"A diagonal in D- set has more than 1 queen.")
    #     return False

    return True


def plot_chessboard(n, queens, dp, dm):
    """Create a chessboard with queens using matplotlib. Image is saved
    in the root directory. Returns the image file name.
    """
    chessboard = np.zeros((n,n))
    chessboard[1::2,0::2] = 1
    chessboard[0::2,1::2] = 1

    # Adjust fontsize for readability
    if n <= 10:
        fontsize = 30
    elif n <= 20:
        fontsize = 10
    else:
        fontsize = 5

    plt.xticks(np.arange(n))
    plt.yticks(np.arange(n))

    plt.imshow(chessboard, cmap='binary')

    # Place queens
    for qb,v in solution.items():
        y = qb // n
        x = qb-n*y
        if v:
            plt.text(x, y, u"\u2655", fontsize=fontsize, ha='center',
                        va='center', color='black' if (x - y) % 2 == 0 else 'white')

    # Draw blocked diagonals
    dp_dict = {}
    dm_dict = {}
    for i in range(n**2):
        y = i // n
        x = i-n*y
        if (y+x) in dp:
            if y+x not in dp_dict:
                dp_dict[y+x] = [[],[]]
                dp_dict[y+x][0] = [x]
                dp_dict[y+x][1] = [y]
            else:
                dp_dict[y+x][0].append(x)
                dp_dict[y+x][1].append(y)
        if (y-x) in dm:
            if y-x not in dm_dict:
                dm_dict[y-x] = [[],[]]
                dm_dict[y-x][0] = [x]
                dm_dict[y-x][1] = [y]
            else:
                dm_dict[y-x][0].append(x)
                dm_dict[y-x][1].append(y)
    # print(dp_dict)
    # print(dm_dict)

    for e in dp_dict:
        if len(dp_dict[e][0]) == 1:
            plt.plot(dp_dict[e][0],dp_dict[e][1],'bs')    
        else:
            plt.plot(dp_dict[e][0],dp_dict[e][1],'b')
    
    for e in dm_dict:
        if len(dm_dict[e][0]) == 1:
            plt.plot(dm_dict[e][0],dm_dict[e][1],'rs')
        else:
            plt.plot(dm_dict[e][0],dm_dict[e][1],'r')

    # Save file in root directory
    file_name = "{}-queens-solution.png".format(n)
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    n = 23
    dp = [39,42,24,12,1,40,20,11,7,27,10,44,38,0]
    dm = [22,16,19,1,-20,-4,14,13,-22,-15,-11,-1,5,-17,-18,18]

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    solution,dp,dm = n_queens(n,dp,dm)

    if is_valid_solution(n, solution):
        write = "Solution is valid."
    else:
        write = "Solution is invalid."
    
    print(write)
    f = open("outHybrid.txt", "a")
    f.write(write+'\n\n')
    f.close()

    file_name = plot_chessboard(n, solution,dp,dm)
    print("Chessboard created. See: {}".format(file_name))
