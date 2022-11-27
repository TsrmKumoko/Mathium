from Complex import Complex

class Matrix:

    def __init__(self, *args):
        self.r = len(args)
        self.c = len(args[0])
        for i in range(1, self.r):
            if len(args[i]) > self.c:
                self.c = len(args[i])
        self.Elem = [[0 for _ in range(self.c)] for _ in range(self.r)]
        for i in range(self.r):
            self.Elem[i] = args[i] + [0] * (self.c - len(args[i]))

    def __add__(self, other):
        if self.r != other.r or self.c != other.c:
            print('Cannot add two matrices together with different dimensions')
            return Matrix([0])
        else:
            Elem = [[0 for _ in range(self.c)] for _ in range(self.r)]
            for i in range(self.r):
                for j in range(self.c):
                    Elem[i][j] = self.Elem[i][j] + other.Elem[i][j]
            return Matrix(*Elem)

    def __sub__(self, other):
        if self.r != other.r or self.c != other.c:
            print('Cannot add two matrices together with different dimensions')
            return Matrix([0])
        else:
            Elem = [[0 for _ in range(self.c)] for _ in range(self.r)]
            for i in range(self.r):
                for j in range(self.c):
                    Elem[i][j] = self.Elem[i][j] - other.Elem[i][j]
            return Matrix(*Elem)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            if self.c != other.r:
                print('Cannot multiply these two matrices')
                return Matrix([0])
            else:
                Elem = [[0 for _ in range(other.c)] for _ in range(self.r)]
                for i in range(self.r):
                    for j in range(other.c):
                        for k in range(self.c):
                            Elem[i][j] += self.Elem[i][k] * other.Elem[k][j]
                return Matrix(*Elem)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, Complex):
            Elem = [[0 for _ in range(self.c)] for _ in range(self.r)]
            for i in range(self.r):
                for j in range(self.c):
                    Elem[i][j] = other * self.Elem[i][j]
            return Matrix(*Elem)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, Complex):
            return self * (1 / other)
        else:
            return NotImplemented

    # def __pow__(self, other):
    #     if isinstance(other, int):
    #         return pow(self, other)
    #     else:
    #         return NotImplemented

    def sub(self, m, n):
        Elem = [[0 for _ in range(self.c - 1)] for _ in range(self.r - 1)]
        for i in range(self.r - 1):
            if i < m:
                p = i
            else:
                p = i + 1
            for j in range(self.c - 1):
                if j < n:
                    q = j
                else:
                    q = j + 1
                Elem[i][j] = self.Elem[p][q]
        return Matrix(*Elem)

def I(n):
    Elem = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        Elem[i][i] = 1
    return Matrix(*Elem)

def det(M):
    if M.r != M.c:
        print('Connot compute determinant')
        return 0
    else:
        Ans = 0
        if M.r == 1:
            Ans = M.Elem[0][0]
        elif M.r == 2:
            Ans = M.Elem[0][0] * M.Elem[1][1] - M.Elem[1][0] * M.Elem[0][1]
        else:
            for i in range(M.r):
                Ans += (-1) ** i * M.Elem[i][0] * det(M.sub(i, 0))
        return Ans

def adjoint(M):
    Elem = [[0 for _ in range(M.c)] for _ in range(M.r)]
    for i in range(M.r):
        for j in range(M.c):
            Elem[i][j] = (-1) ** (i + j) * det(M.sub(j, i))
    return Matrix(*Elem)

def inv(M):
    if det(M) == 0:
        print('Do not have inverse matrix')
        return Matrix([0])
    else:
        return adjoint(M) / det(M)

def pow(M, n):
    if M.r != M.c:
        print('Cannot compute power')
        return Matrix([0])
    if n == -1:
        return inv(M)
    elif n == 0:
        return I(M.r)
    elif n > 0:
        Ans = M
        for i in range(n - 1):
            Ans = Ans * M
        return Ans
    elif n < -1:
        return pow(inv(M), -n)

