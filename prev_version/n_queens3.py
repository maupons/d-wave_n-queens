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
    """Returns a potential solution to the n-queens problem in a list of sets,
    each containing constraint IDs representing a queen's location.

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. Defaults to dwave-system's LeapHybridSampler.
    """
    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    bqm.offset = 2*n

    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
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

    bqm.add_variable(9, 4)
    bqm.add_variable(3, 4)
    # print(bqm)

    # sampler = EmbeddingComposite(DWaveSampler())
    # sampler = LeapHybridSampler()
    sampler = neal.SimulatedAnnealingSampler()

    # sampleset = sampler.sample(bqm, label='Example - N Queens')
    # sampleset = sampler.sample(bqm, label=f'{n} QueensD')
    sampleset = sampler.sample(bqm, num_reads=1, label=f'{n} QueensD - Simulated Annealing')
    sample = sampleset.first.sample
    print("Sample:")
    print(sample)
    return sample

    # print("Sampleset:")
    # print(sampleset)
    # print(type(sampleset))

    # Inspect
    # dwave.inspector.show(sampleset)

    # return [subsets[i] for i in sample if sample[i]]
    # return list()

def is_valid_solution(n, solution):
    """Check that solution is valid by making sure all constraints were met.

    Args:
        n: Number of queens in the problem.

        solution: A list of sets, each containing constraint IDs that represent
                  a queen's location.
    """
    board = np.zeros((n,n))
    ldp = []                            # Keep track of queens on D+
    ldm = []                            # Keep track of queens on D-
    # solution[1] = 1
    # print(solution)

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

    # board = np.zeros((n,n))
    # board[0,0] = 1
    # board[1,0] = 1
    # board[2,0] = 1
    # board[3,0] = 1
    print(board)

    # sum_cols = np.sum(board,axis=0)
    # print(sum_cols)
    # sum_rows = np.sum(board,axis=1)
    # print(sum_rows)
    # for i in range(len(sum_cols)):
    #     ncol_queens = sum_cols[i]
    #     nrow_queens = sum_rows[i]
    #     if ncol_queens != 1:
    #         print(f"Column {i} has {ncol_queens} queens.")
    #         # return False
    #     elif nrow_queens != 1:
    #         print(f"Row {i} has {nrow_queens} queens.")
    #         # return False
    # # return True
    
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


def plot_chessboard(n, queens):
    """Create a chessboard with queens using matplotlib. Image is saved
    in the root directory. Returns the image file name.
    """
    chessboard = np.zeros((n,n))
    chessboard[1::2,0::2] = 1
    chessboard[0::2,1::2] = 1
    print(chessboard)

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
    # plt.imshow(chessboard, cmap='gray')
    # plt.imshow(chessboard)

    # Place queens
    for qb,v in solution.items():
        y = qb // n
        x = qb-n*y
        if v:
            plt.text(x, y, u"\u2655", fontsize=fontsize, ha='center',
                        va='center', color='black' if (x - y) % 2 == 0 else 'white')
            print(x,y)
        # Try to show blocked diagonals
        if x+y == 3:
            plt.plot([0,1,2,3])

    # for subset in solution:
    #     x = y = -1
    #     for constraint in subset:
    #         if constraint < n:
    #             x = constraint
    #         elif constraint >= n and constraint < 2*n:
    #             y = np.abs(constraint - (2*n - 1)) # Convert constraint ID to row index

    #     if x != -1 and y != -1:
    #         plt.text(x, y, u"\u2655", fontsize=fontsize, ha='center',
    #                  va='center', color='black' if (x - y) % 2 == 0 else 'white')

    # Save file in root directory
    file_name = "{}-queens-solution.png".format(n)
    plt.savefig(file_name)

    return file_name

def get_sanitized_input():
    while True:
        print("Enter the number of queens to place (n > 0):")
        n = input()

        try:
            n = int(n)
            if n <= 0:
                print("Input must be greater than 0.")
                continue
            if n >= 200:
                # Run but give a warning
                print("Problems with large n will run very slowly.")

        except ValueError:
            print("Input type must be int.")
            continue

        return n

if __name__ == "__main__":
    # n = get_sanitized_input()
    n = 4

    if n > 20:
        print("Solution image is large and may be difficult to view.")
        print("Plot settings in plot_chessboard() may need adjusting.")

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    solution = n_queens(n)

    if is_valid_solution(n, solution):
        print("Solution is valid.")
    else:
        print("Solution is invalid.")

    file_name = plot_chessboard(n, solution)
    print("Chessboard created. See: {}".format(file_name))
