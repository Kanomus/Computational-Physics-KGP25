# linear equations solver using matrices

'''
[b1]     [a11, a12, a13, a14, ...]   [x1]
[b2]     [a21, a22, a23, a24, ...]   [x2]
[b3]  =  [a31, a32, a33, a34, ...]   [x3]
[b4]     [a41, a42, a43, a44, ...]   [x4]
[..]     [..., ..., ..., ..., ...]   [..]

[B] = [A][X]
[x] = [A]^-1[B]
'''

class Matrix:
    def __init__(self, *rows) :
        if len(rows) == 1 and isinstance(rows[0], list) and all(isinstance(row, list) for row in rows[0]):
            rows = tuple(rows[0])
        n = len(rows)
        m = len(rows[0])
        for row in rows:
            if len(row) != m: raise ValueError("Matrix dimensions not consistent")
        self.dimensions = (n, m)
        self.vals = [*rows]
        self.isSquare = (self.dimensions[0] == self.dimensions[1])
    
    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            if len(idx) != 2: raise ValueError(f"Expected 2 positional arguments, got {len(idx)}")
            return self.vals[idx[0]][idx[1]]
        return self.vals[idx]
    
    def __setitem__(self, idx, val):
        if isinstance(idx, tuple):
            if len(idx) != 2: raise ValueError(f"Expected 2 positional arguments, got {len(idx)}")
            self.vals[idx[0]][idx[1]] = val
        elif isinstance(val, list):
            if len(val) != self.dimensions[1]: raise ValueError(f"Matrix row length mismatch, expected {self.dimensions[1]}, got {len(val)}")
            self.vals[idx] = val

    
    def __str__(self) :
        dimension_message = f"{self.dimensions[0]} x {self.dimensions[1]} matrix"
        matrix_representation = ""
        for i in range(self.dimensions[0]):
            matrix_representation += "\n["
            for j in range(self.dimensions[1]):
                matrix_representation += str(self.vals[i][j])
                if j != self.dimensions[1] - 1 : matrix_representation += ", "
            matrix_representation += "]"
        return matrix_representation
    
    def __iter__(self):
        for row in self.vals:
            yield row
    
    def __add__(self, other) :
        if not isinstance(other, Matrix): raise TypeError(f"Cannot add matrix with {type(other)}")
        if self.dimensions != other.dimensions: raise TypeError(f"Cannot add matrices with dimensions {self.dimensions} & {other.dimensions}")
        n,m = self.dimensions
        return Matrix([[self.vals[i][j] + other.vals[i][j] for j in range(m)] for i in range(n)])
    
    def __sub__(self, other) :
        if not isinstance(other, Matrix): raise TypeError(f"Cannot subtract matrix from {type(other)}")
        if self.dimensions != other.dimensions: raise TypeError(f"Cannot subtract matrices with dimensions {self.dimensions} & {other.dimensions})")
        n,m = self.dimensions
        return Matrix([[self.vals[i][j] - other.vals[i][j] for j in range(m)] for i in range(n)])
    
    def __mul__(self, other) :
        if isinstance(other, (int, float)):
            n,m = self.dimensions
            return Matrix([[self.vals[i][j]*other for j in range(m)] for i in range(n)])
        if isinstance(other, Matrix):
            n,m = self.dimensions
            o,p = other.dimensions
            if (n,m) != (o,p): raise TypeError("Cannot element-wise multiply matrices of different sizes (Use @ for matrix multiplication)")
            return Matrix([[self.vals[i][j] * other.vals[i][j] for j in range(m)] for i in range(n)])
    
    def __rmul__(self, other) :
        if isinstance(other, (int, float)):
            n,m = self.dimensions
            return Matrix([[self.vals[i][j] * other for j in range(m)] for i in range(n)])
    
    def __matmul__(self, other) :
        if not isinstance(other, Matrix) : raise TypeError(f"Cannot do matrix multiplication with {type(other)}")
        return multiply(self, other)

    def fill(self, val):
        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                self.vals[i][j] = val
    
    def trace(self):
        if not self.isSquare: raise TypeError("Trace is only defined for square matrices")
        tr = 0
        for i in range(self.dimensions[0]):
            tr += self.vals[i][i]
        return tr
    
    def minor(self, row, col):
        return Matrix(*[self.vals[i][:col-1]+self.vals[i][col:] for i in range(self.dimensions[0]) if i!=(row-1)])

    def transpose(self):
        return Matrix(*[[self.vals[i][j] for i in range(self.dimensions[0])] for j in range(self.dimensions[1])])
    
    def swap_rows(self, row1, row2):
        temp = self.vals[row1-1]
        self.vals[row1-1] = self.vals[row2-1]
        self.vals[row2-1] = temp

    def swap_columns(self, column1, column2):
        temp = [self.vals[i][column1 - 1] for i in range(self.dimensions[0])]
        for i in range(self.dimensions[0]):
            self.vals[i][column1 - 1] = self.vals[i][column2 - 1]
            self.vals[i][column2 - 1] = temp[i]
    
    def inverse(self) :
        return NotImplementedError

def zeroes(n : int, m : int) -> Matrix :
    return Matrix(*[[0 for _ in range(m)] for _ in range(n)])

def identity(n : int) -> Matrix :
    return Matrix(*[[1 if i==j else 0 for j in range(n)] for i in range(n)])

def multiply(A : Matrix, B : Matrix) -> Matrix :
    n = A.dimensions[0]
    m = B.dimensions[1]
    if A.dimensions[1] != B.dimensions[0]:
        raise TypeError(f"Cannot multiply matrices with dimensions {n}x{A.dimensions[1]} & {B.dimensions[0]}x{m}  (Use * for element-wise multiplication")
    o = A.dimensions[1]

    C = zeroes(n, m)
    for i in range(n):
        for j in range(m):
            for k in range(o):
                C.vals[i][j] += A.vals[i][k] * B.vals[k][j]
    return C

def LU_decompose(matrix : Matrix) -> tuple[Matrix, Matrix, Matrix] :
    if(matrix.dimensions[0] != matrix.dimensions[1]): raise TypeError("Cannot decompose non-square matrices into LU form")
    size = matrix.dimensions[0]
    
    # Find L & U such that A = L@U, and L has all diagonal elements = 1

    # Initialize L as identity
    # For each column:
        # Find multiplier to be put in L to make element of U in that row = 0
        # Find the corresponding row of U by multiplying with element of L and subtracting from A
        # Uk,j = Ak,j - sum(0, k-1)Lk,s Us,j, j>=k
        # Li,k = (Ai,k - sum(0, k-1)Li,s Us,k)/Uk,k, i>K
    
    # P is the permutation matrix to keep track of row swaps when pivoting
    P = identity(size)
    L = identity(size)
    U = Matrix(*[row[:] for row in matrix.vals])

    for i in range(size):
        pivot = max(range(i,size), key = lambda r: U[r][i])
        if U[pivot][i] == 0: raise ValueError("The matrix is singular and thus cannot be decomposed")

        # swap rows to put pivot on the top row
        U[i], U[pivot] = U[pivot], U[i]
        P[i], P[pivot] = P[pivot], P[i]
        if i>0:
            L[i][:i], L[pivot][:i] = L[pivot][:i], L[i][:i]
        
        # Eliminate terms from U
        for j in range(i+1, size):
            factor = U[j][i]/U[i][i]
            L[j][i] = factor
            
            for k in range(i, size):
                U[j][k] -= factor*U[i][k]
        
    return P, L, U
