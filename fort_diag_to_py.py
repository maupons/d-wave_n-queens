# D + needs to be reduced by 2
# D - does not need conversion


a = '8          19          16          22          -1         -17           3         -23          14           7         -15         -16          -8         -24           1           7          -9          18          23          21          23           1          12         -10          18           0           0           0           0           0           0           0'
n = 12

a = a.split()[0:n]
res = ','.join(a)
print(res)

a = [35,32,42,18,39,50,45,33,28,9,46,6,3,8,17,4,43,7,44,2]
b = [str(i-2) for i in a]
b = ','.join(b)
print(b)