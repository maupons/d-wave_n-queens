import numpy as np

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

    print(board)

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





a = {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0, 36: 1, 37: 0, 38: 0, 39: 0, 40: 0, 41: 0, 42: 0, 43: 0, 44: 0, 45: 0, 46: 0, 47: 1, 48: 0, 49: 0, 50: 0, 51: 0, 52: 1, 53: 0, 54: 0, 55: 0, 56: 1, 57: 0, 58: 0, 59: 0, 60: 0, 61: 0, 62: 0, 63: 0}
print('valid sol?: ', is_valid_solution(8, a))