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
# 
# - removed qubits for blocked diagonals
# 

from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt
from dimod import BinaryQuadraticModel, ExactSolver
from dwave.samplers import SimulatedAnnealingSampler, SteepestDescentSolver \
, TabuSampler, TreeDecompositionSolver, RandomSampler
from time import time
import sys
import datetime as dt


def n_queens(n,dp,dm,itr, sampler=None):
    """Returns a potential solution to the n-queens problem in a dictionary, where
    the keys are qubit id and value on or off

    Args:
        n: Number of queens to place.

        sampler: A binary quadratic model sampler. 
        Can use: 1) LeapHybridSampler
                 2) DWaveSampler - QPU
                 3) SimulatedAnnealingSampler
    """

    # d = len(dp) + len(dm)
    bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')
    w = 2

    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
        if (ri+ci) in dp:
            continue
        elif (ri-ci) in dm:
            continue
        else:
            bqm.add_variable(i, -w)
            
        for j in range(i):
            if j not in bqm.linear:
                continue
            rj = j // n
            cj = j-n*rj
            if rj == ri:
                bqm.add_interaction(i, j, w)
            if cj == ci:
                bqm.add_interaction(i, j, w)
            if abs(ri-rj) == abs(ci-cj):
                bqm.add_interaction(i, j, w)
    # print(bqm)
    # print('LINEAR')
    # print(bqm.linear)
    # print(len(bqm.linear))
    # print('QUADRATIC')
    # print(bqm.quadratic)
    # print(len(bqm.quadratic))


    print(f'Classical solver started with {itr} reads...')
    start_time = time()

    # Exact Solver
    # sampler = ExactSolver()
    # sampleset = sampler.sample(bqm)

    # Heuristic  Solvers
    sampler = SimulatedAnnealingSampler()
    # sampler = SteepestDescentSolver()
    # sampler = TabuSampler()
    # sampler = TreeDecompositionSolver()
    # sampler = RandomSampler()
    print(sampler)
    sampleset = sampler.sample(bqm, num_reads=itr)

    py_time = time()-start_time
    return sampleset, sampler, py_time


def is_valid_solution(n, solution):
    """Check that solution is valid by making sure all constraints were met.

    Args:
        n: Number of queens in the problem.

        solution: A dictionary of qubits, key = qubit id, value = 0 or 1
    """
    board = np.zeros((n,n))
    ldp = []                            # Keep track of queens on D+
    ldm = []                            # Keep track of queens on D-

    for qb,v in solution.items():       # Build & Check diag/anti-diag constraints
        r = qb // n
        c = qb-n*r
        if v:
            board[r,c] = v
            dp = r+c
            dm = r-c
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

    for i in range(n):                              # Check row/col constraints
        sum_row = sum(board[i,:])
        if sum_row != 1:
            print(f"Row {i} has {sum_row} queens.")
            return False
        sum_col = sum(board[:,i])
        if sum_col != 1:
            print(f"Column {i} has {sum_col} queens.")
            return False

    return True


def plot_chessboard(n,queens, dp, dm, ruta):
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
    for qb,v in queens.items():
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
    file_name = f"{ruta}figs/{n}-queens-solution.png"
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    ruta = 'data/c_data/'
    # n = int(sys.argv[1])
    n = 50
    itr = 30000

    # n 10
    # dp = [0, 1, 3, 6, 8, 12, 16, 17, 18]
    # dm = [-9, -8, -7, -4, 0, 5, 6, 8, 9]

    # n 11
    # dp = [0, 1, 3, 6, 9, 12, 14, 16, 19, 20]
    # dm = [-10, -9, -8, -7, -4, 5, 6, 8, 9, 10]

    # n 12
    # dp = [0, 1, 2, 7, 9, 10, 12, 18, 19, 21, 22]
    # dm = [-11, -10, -8, -7, -6, -3, 7, 8, 9, 10, 11]

    # n 20
    # dp = [0, 1, 2, 3, 4, 5, 6, 15, 18, 23, 25, 26, 28, 29, 31, 34, 36, 37, 38]
    # dm = [-19, -18, -16, -15, -14, -13, -12, -9, -6, -4, 3, 11, 13, 14, 15, 16, 17, 18, 19]

    # n 30
    # dp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 14, 19, 22, 23, 27, 33, 42, 43, 44, 45, 47, 50, 51, 53, 54, 55, 56, 57, 58]
    # dm = [-29, -28, -27, -26, -25, -24, -22, -20, -19, -17, -15, -13, -11, -7, 3, 4, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 27, 28, 29]

    # n 40
    # dp = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 14, 15, 16, 18, 20, 22, 26, 28, 30, 32, 45, 53, 54, 55, 59, 60, 63, 65, 66, 67, 69, 70, 71, 73, 74, 75, 76, 77, 78]
    # dm = [-39, -38, -37, -35, -34, -33, -30, -29, -28, -27, -26, -24, -23, -22, -20, -16, -15, -14, -11, -8, -4, 2, 13, 21, 22, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

    # n 50
    dp = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 20, 22, 23, 24, 26, 35, 40, 45, 53, 55, 58, 59, 67, 72, 73, 74, 77, 81, 82, 84, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98]
    dm = [-49, -48, -47, -46, -45, -44, -43, -42, -41, -40, -39, -38, -37, -35, -33, -27, -26, -25, -22, -20, -9, -8, -5, -3, -2, 0, 1, 11, 18, 22, 25, 27, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 42, 43, 44, 45, 47, 48, 49]

    # n 60
    # dp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15, 17, 18, 19, 20, 21, 22, 25, 27, 30, 31, 32, 33, 35, 36, 39, 41, 42, 43, 60, 72, 80, 81, 84, 85, 88, 93, 94, 95, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109, 110, 111, 112, 113, 114, 115, 117, 118]
    # dm = [-59, -58, -57, -56, -55, -54, -53, -52, -51, -50, -49, -48, -47, -45, -44, -42, -41, -39, -38, -37, -35, -33, -28, -25, -23, -21, -7, -2, 1, 4, 6, 7, 10, 11, 23, 25, 26, 31, 35, 36, 38, 40, 41, 42, 43, 44, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]

    # n 70
    # dp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 21, 24, 25, 28, 30, 33, 37, 38, 40, 41, 45, 47, 53, 57, 60, 63, 65, 78, 86, 87, 90, 92, 94, 95, 97, 99, 103, 104, 106, 108, 110, 114, 116, 117, 119, 120, 122, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138]
    # dm = [-69, -68, -67, -66, -65, -64, -63, -62, -61, -60, -59, -58, -57, -56, -55, -54, -53, -52, -51, -48, -44, -43, -42, -41, -40, -38, -36, -35, -28, -27, -20, -17, -13, -10, -2, 11, 22, 25, 26, 34, 35, 36, 37, 38, 39, 40, 41, 43, 45, 46, 47, 48, 49, 50, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 68, 69]


    d = len(dp) + len(dm)

    for ix in range(1):
        print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
        sampleset, sampler, py_time = n_queens(n,dp,dm,itr)

        # Extract run info
        sample = sampleset.first.sample
        energy = sampleset.first.energy
        num_itr = len(sampleset)
        nvars = len(sample)
        df = sampleset.to_pandas_dataframe()
        nsols = df[df["energy"] == -2*n]['num_occurrences'].sum()
        p_sol = psol = nsols / num_itr
        sp_name = str(sampler).split(".")[4].split()[0]
        ID = time()
        
        if is_valid_solution(n,sample):
            solved = "YES"
        else:
            solved = "NO"

        # print('Solution\n', sample)
        print('sampleset.info: ',sampleset.info)
        print('sampler.properties: ',sampler.properties)
        print('sampler',sp_name)
        print('py_time', py_time)
        print('nsols', nsols)
        print('num_itr', num_itr)
        print('p_sol', p_sol)
        print('energy', energy)
        print('nvars', nvars)        
        print('- Solved - ', solved)

        # Write sampleset and solutions to file
        # f1 = open(f"{ruta}sp/{n}_sols_{ID}.txt", "w")
        # for spl in sampleset:
        #     f1.write(str(spl)+'\n')
        # f1.close()

        # f22 = open(f"{ruta}sp/{n}_samplesetPD_{ID}.txt", "w")
        # f22.write(df.to_string())
        # f22.close()

        # f3 = open(f"{ruta}logs_vsH.txt", "a")
        # f3.write(f'{n}-q; {d}-d config\n')
        # f3.write(str(sample)+'\n')
        # f3.write('py_time ' + str(py_time) + '\n')
        # f3.write(str(sampleset.info)+'\n')
        # f3.write('Solution - ' + solved + '\n\n')
        # f3.close()

        line = f'{n}   {d}   {nvars}   {num_itr}   {p_sol}   {sp_name}   '\
               f'{py_time*10**3}   {energy}   {solved}\n'
        f4 = open(f"{ruta}time.txt", "a")
        f4.write(line)
        f4.close()

        # f2 = open(f"{ruta}sp/{n}_sampleset_{ID}.txt", "w")
        # f2.write(str(sampleset))
        # f2.close()

        # file_name = plot_chessboard(n,sample,dp,dm,ruta)
        # print("Chessboard created. See: {}".format(file_name))
