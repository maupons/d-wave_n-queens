import numpy as np

n = 4
n2 = n**2

bqm = np.zeros((n2,n2))
for i in range(n2):
    ri = i // n
    ci = i-n*ri
    bqm[i,i] = -2
    for j in range(i):
        rj = j // n
        cj = j-n*rj
        if rj == ri:
            bqm[j,i] = 2
        if cj == ci:
            bqm[j,i] = 2
        if abs(ri-rj) == abs(ci-cj):
            bqm[j,i] = 2

print(bqm)

        
        
