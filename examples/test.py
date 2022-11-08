from dimod import BinaryQuadraticModel
import numpy as np
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt

# bqm = BinaryQuadraticModel({}, {}, 0, 'BINARY')

# bqm.add_variable(0, -1)
# bqm.add_variable(0, -1)
# bqm.add_variable(0, -1)

# bqm.add_variable(1, -1)


# bqm.add_interaction(0, 1, 2)
# bqm.add_interaction(1, 0, 1)
# bqm.add_interaction(6, 7, 1)
# print(bqm)



# n = 4
# board = np.zeros((n,n))
# print(board)
# for c in range(n):
#     for r in range(c):
#         board[r,c] = 1
# print(board)

# a = np.arange(n**2).reshape((n, n))
# print(a)


# n = 5
# n2 = n**2
# for i in range(n2):
#     row = i // n
#     col = i-n*row
#     print(i,":",row,col)



# mylist = [5, 3, 5, 2, 1, 6, 6, 4] # the original list of integers with duplicates
# newlist = [] # empty list to hold unique elements from the list
# duplist = [] # empty list to hold the duplicate elements from the list
# for i in mylist:
#     newlist.append(i)
#     if i not in newlist:
#         print(f'{i} already in list')
#     else:
#         duplist.append(i) # this method catches the first duplicate entries, and appends them to the list
# # The next step is to print the duplicate entries, and the unique entries
# print("List of duplicates", duplist)
# print("Unique Item List", newlist)



n = 10
a = np.arange(n**2).reshape((n, n))
print(a)

# Solution with 15 diagonals blocked
# dplus  = [20,19,5,7,3,2,17]
# dminus = [9,-1,3,5,0,8,-9,-8]
b = np.zeros((n,n))
xlist = [2,7,5,8,0,9,4,6,1,3]
b[0,2] = 1
b[1,7] = 1
b[2,5] = 1
b[3,8] = 1
b[4,0] = 1
b[5,9] = 1
b[6,4] = 1
b[7,6] = 1
b[8,1] = 1
b[9,3] = 1
print(b)

def plot_chessboard(n, xlist):
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
    y = -1
    for x in xlist:
        y += 1
        print(x,y)
        plt.text(x, y, u"\u2655", fontsize=fontsize, ha='center',
                    va='center', color='black' if (x - y) % 2 == 0 else 'white')
            # print(x,y)
        # Try to show blocked diagonals
        # if x+y == 3:
        #     plt.plot([0,1,2,3])

    # Save file in root directory
    file_name = f"{n}-queens-solution-15d.png"
    plt.savefig(file_name)

plot_chessboard(n, xlist)


{0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0, 36: 1, 37: 0, 38: 0, 39: 0, 40: 0, 41: 0, 42: 0, 43: 0, 44: 0, 45: 0, 46: 0, 47: 1, 48: 0, 49: 0, 50: 0, 51: 0, 52: 1, 53: 0, 54: 0, 55: 0, 56: 1, 57: 0, 58: 0, 59: 0, 60: 0, 61: 0, 62: 0, 63: 0}