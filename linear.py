# Script to solve system of linear equations

import matrix
from matrix import Matrix

def solve_system(A : Matrix, B : Matrix) -> Matrix :
    n, m = A.dimensions
    d, e = B.dimensions
    if not n == d: raise ValueError(f"Dimensions of A and B are incompatible, got ({n}x{m}) and ({d}x{e})")

    P, L, U = matrix.LU_decompose(A)
    # AX = B,   PA = LU,    PAX = PB = LUX, UX = Y, LY = PB
    # Solve LY = B and UX = Y using forward and backward substitution respectively
    Y = matrix.zeroes(d,e)
    X = matrix.zeroes(d,e)

    # LY = PB = B'
    B = P@B
    for k in range(e):
        for i in range(n):
            Y[i,k] = B[i,k]
            for j in range(i):
                Y[i,k] -= L[i,j]*Y[j,k]

    # UX = Y
    for k in range(e):
        for i in range(1, n+1):
            X[-i,k] = Y[-i,k]
            for j in range(1, i):
                X[-i,k] -= U[-i,-j]*X[-j,k]
            X[-i,k] /= U[-i,-i]
    
    return X