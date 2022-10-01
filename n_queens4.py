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
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector
import neal


def n_queens(n, sampler=None):
    """Returns a potential solution to the n-queens problem in a dictionary, where
    the keys are qubit id and value on or off

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. 
        Can use: 1) LeapHybridSampler
                 2) DWaveSampler - QPU
                 3) SimulatedAnnealingSampler
    """
    # dplus  = [20,19,5,7,3,2,17]
    dp  = [18,17,3,5,1,0,15]
    dm = [9,-1,3,5,0,8,-9,-8]

    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    bqm.offset = 2*n

    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
        if (ri+ci) in dp:
            # print('i',i,'r+c',ri+ci)
            # print('i',i,'r',ri,'c',ci,'r+c',ri+ci)
            bqm.add_variable(i, 2)
        elif (ri-ci) in dm:
            # print('i',i,'r-c',ri-ci)
            # print('i',i,'r',ri,'c',ci,'r-c',ri-ci)
            bqm.add_variable(i, 2)
        else:
            bqm.add_variable(i, -2)
            
        for j in range(i):
            rj = j // n
            cj = j-n*rj
            if rj == ri:
                bqm.add_interaction(i, j, 2)
            if cj == ci:
                bqm.add_interaction(i, j, 2)
            if abs(ri-rj) == abs(ci-cj):
                bqm.add_interaction(i, j, 2)
    # print(bqm)

    # sampler = EmbeddingComposite(DWaveSampler())              # QPU
    # sampler = LeapHybridSampler()                             # Hybrid Solver
    sampler = neal.SimulatedAnnealingSampler()                  # CPU

    # sampleset = sampler.sample(bqm, label=f'{n} QueensD')
    sampleset = sampler.sample(bqm, num_reads=1, label=f'{n} QueensD')
    sample = sampleset.first.sample

    print("Sample:\n", sample)
    # print("Sampleset:\n", sampleset)

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
            # print(x,y)
        # Try to show blocked diagonals
        # if x+y == 3:
        #     plt.plot([0,1,2,3])

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
    print(dp_dict)
    print(dm_dict)
    # for e in dp_dict:
    #     plt.plot(dp_dict[e][0],dp_dict[e][1],'bs')
    
    # for e in dm_dict:
    #     plt.plot(dm_dict[e][0],dm_dict[e][1],'rs')

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


    # x2 = list(range(1,n))
    # print(x2)
    # plt.plot(x2,dm_dict[-1],'r')
    
    # plt.plot(dp_dict[0],'b')
    # plt.plot(dp_dict[1],'b')



    # Save file in root directory
    file_name = "{}-queens-solution.png".format(n)
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    n = 10

    if n > 20:
        print("Solution image is large and may be difficult to view.")
        print("Plot settings in plot_chessboard() may need adjusting.")

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    solution,dp,dm = n_queens(n)

    if is_valid_solution(n, solution):
        print("Solution is valid.")
    else:
        print("Solution is invalid.")

    file_name = plot_chessboard(n, solution,dp,dm)
    print("Chessboard created. See: {}".format(file_name))
