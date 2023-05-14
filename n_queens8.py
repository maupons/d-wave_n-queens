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
import sys
import datetime as dt


def n_queens(n,dp,dm,itr,ruta, sampler=None):
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
    # bqm.offset = 2*n
    # itr = 500
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

    print(f'solver started with {itr} reads...')
    start_time = time()
    
    # QPU
    sampler = EmbeddingComposite(DWaveSampler())
    embedding_time = time()-start_time
    sampleset = sampler.sample(bqm, num_reads=itr, label=f'n {n}, time {start_time}')
    
    py_time = time()-start_time
    print('py_time', py_time)
    print('sampleset.done', sampleset.done())
    cnt = 0
    while not sampleset.done():
        cnt = cnt + 1
        continue
    print('wait itr: ', cnt)
    print('sampleset.done', sampleset.done())
    if sampleset.done():
        done_time = time()-start_time
    print('done_time', done_time)

    # Extract run info
    # print('sampleset.info:\n', sampleset.info)
    embedding = sampleset.info['embedding_context']['embedding']
    chain_break_method = sampleset.info['embedding_context']['chain_break_method']
    chain_strength = sampleset.info['embedding_context']['chain_strength']
    n_vars = len(embedding.keys())
    n_qb = sum(len(chain) for chain in embedding.values())
    energy = sampleset.first.energy
    print(f"Number of logical variables: {n_vars}")
    print(f"Number of physical qubits used in embedding: {n_qb}")
    print('energy', energy)

    # ruta = 'data/n_data/'
    f1 = open(f"{ruta}/sp/{n}_sols_{start_time}.txt", "w")
    for sample in sampleset:
        f1.write(str(sample)+'\n')
    f1.close()

    f2 = open(f"{ruta}/sp/{n}_sampleset_{start_time}.txt", "w")
    f2.write(str(sampleset))
    f2.close()

    df = sampleset.to_pandas_dataframe()
    f22 = open(f"{ruta}/sp/{n}_samplesetPD_{start_time}.txt", "w")
    f22.write(df.to_string())
    f22.close()    
    nsols = df[df["energy"] == -2*n]['num_occurrences'].sum()
    psol = nsols / itr

    sample = sampleset.first.sample
    f3 = open(f"{ruta}logsQPU.txt", "a")
    f3.write(f'{n}-q; {d}-d config\n')
    f3.write('embedding_time ' + str(embedding_time) + '\n')
    f3.write(str(sample)+'\n')
    f3.write('py_time ' + str(py_time) + '\n')
    f3.write(str(sampleset.info)+'\n')
    f3.close()

    f4 = open(f"{ruta}time.txt", "a")
    line = f'{n}   {d}   {itr}   {psol}   {energy}   {n_vars}   {n_qb}   ' \
           f'{chain_strength}   {chain_break_method}   ' \
           f'{embedding_time*10**3}   {py_time*10**3}   {done_time*10**3}'
    for v in sampleset.info["timing"]:
        line = line + '   ' + str(sampleset.info["timing"][v]*10**-3)
    f4.write(line + '\n')
    f4.close()

    # print('Sampleset:')
    # print(sampleset)
    
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


def plot_chessboard(n,s,dp,dm,ruta):
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
    # file_name = "{}-queens-solution.png".format(n)
    file_name = f"{ruta}/figs/{n}-queens-solution.png"
    plt.savefig(file_name)

    return file_name

if __name__ == "__main__":

    ruta = 'data/n_data/'
    itr = 10
    n = int(sys.argv[1])
    dp = []
    dm = []

    print("Trying to place {n} queens on a {n}*{n} chessboard.".format(n=n))
    s,dp,dm = n_queens(n,dp,dm,itr,ruta)

    if is_valid_solution(n,s):
        write = "Solution - VALID"
    else:
        write = "Solution - NOT VALID"
    print(write)

    f = open(f"{ruta}logsQPU.txt", "a")
    f.write(write +'\n\n')
    f.close()

    file_name = plot_chessboard(n,s,dp,dm,ruta)
    print("Chessboard created. See: {}".format(file_name))
