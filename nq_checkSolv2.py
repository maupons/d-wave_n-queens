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

def get_energy_v1(n,dp,dm,solution):
    e = 0
    w = 2

    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
        if (ri+ci) in dp:
            e += w*solution[i]
            # bqm.add_variable(i, w)
        elif (ri-ci) in dm:
            e += w*solution[i]
            # bqm.add_variable(i, w)
        else:
            # bqm.add_variable(i, -w)
            e -= w*solution[i]
            
        for j in range(i):
            rj = j // n
            cj = j-n*rj
            if rj == ri:
                # bqm.add_interaction(i, j, w)
                e += w*solution[i]*solution[j]
            if cj == ci:
                # bqm.add_interaction(i, j, w)
                e += w*solution[i]*solution[j]
            if abs(ri-rj) == abs(ci-cj):
                # bqm.add_interaction(i, j, w)
                e += w*solution[i]*solution[j]
    return e

def get_energy_v2(n,dp,dm,solution):
    w = 2

    x = list(solution.values())
    x = np.array(x)
    x.shape = (n**2,1) 
    
    bqm = np.zeros((n**2, n**2))
    for i in range(n**2):
        ri = i // n
        ci = i-n*ri
        if (ri+ci) in dp:
            bqm[i,i] = w
        elif (ri-ci) in dm:
            bqm[i,i] = w
        else:
            bqm[i,i] = -w
            
        for j in range(i):
            rj = j // n
            cj = j-n*rj
            if rj == ri:
                bqm[j,i] = w
            if cj == ci:
                bqm[j,i] = w
            if abs(ri-rj) == abs(ci-cj):
                bqm[j,i] = w
    # np.set_printoptions(threshold=np.inf)
    # print(bqm)

    # row = ''
    # for i in bqm:
    #     for j in i:
    #         row += f'{j} '
    #     print(row)
    #     row = ''

    e = np.matmul(x.transpose(),np.matmul(bqm,x))[0,0]
    return e

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

    return True


def plot_chessboard(n,s, dp, dm):
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
    for qb,v in s.items():
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

    n = 6
    dp = []
    dm = []
    s = {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 1, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 0, 19: 0, 20: 0, 21: 1, 22: 0, 23: 0, 24: 0, 25: 1, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    print(s)

    if is_valid_solution(n, s):
        write = "Solution is valid."
    else:
        write = "Solution is invalid."
    
    e1 = get_energy_v1(n,dp,dm,s)
    print('energy_v1:', e1)

    e2 = get_energy_v2(n,dp,dm,s)
    print('energy_v2:', e2)


    print(write)
    file_name = plot_chessboard(n,s,dp,dm)
    print("Chessboard created. See: {}".format(file_name))
