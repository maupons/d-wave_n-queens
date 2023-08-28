# 
# Converte diagonals from column to list 
# 

import numpy as np

dp = list(np.genfromtxt('d.txt', skip_header=1, dtype=int, usecols=(0)))
dm = list(np.genfromtxt('d.txt', skip_header=1, dtype=int, usecols=(1)))

print('d =', len(dp)+len(dm))
print(f'dp = {dp}')
print('dm =', dm)