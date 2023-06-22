# D + needs to be reduced by 2
# D - does not need conversion

ndp = 8
dp = '3          17           4          22           9          13          12          18          13          16           4           9          19          13          16           0           0'

ndm = 9
dm = '4          -4          -7         -10           9          -9           8          10           7          -1           2          -5           6          -5           5           0           0'

dp = dp.split()[0:ndp]
dp = [int(i) for i in dp]
dp = [i-2 for i in dp]


dm = dm.split()[0:ndm]
dm = [int(i) for i in dm]

print('dp', dp)
print('dm', dm)


# b = [str(i-2) for i in dp]
# b = ','.join(b)
# print(b)

# dp = dp.split()[0:ndp]
# # res = ','.join(dp)
# # print(res)
# dp = [int(x) for x in dp]
# print(dp)

# # dp = [16,10,9,2,8]
# # b = [str(i-2) for i in dp]
# dp = [i-2 for i in dp]
# print(dm)
# # b = [str(i-2) for i in dp]
# # b = ','.join(b)
# # print(b)