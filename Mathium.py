from math import sqrt
from copy import copy
from typing import Union
import IPython.display as ipy

ascending = False
auto_reduce = True

class MathiumObject(object):
    
    def __init__(self) -> None:
        pass

    def sign(self) -> bool:
        return True

    def LaTeX(self) -> str:
        return str(self)

class Real(MathiumObject):

    def __init__(self, value: Union[int, float] = 0) -> None:
        if isinstance(value, (int, float)):
            self.__value = value
        elif isinstance(value, Real):
            self = value
        else:
            raise TypeError(f'Cannot initialize a Real number with {value.__class__}.')

    def __eq__(self, other) -> bool:
        if isinstance(other, Real):
            return self.__value == other.__value
        elif isinstance(other, (int, float)):
            return self.__value == other
        else:
            return NotImplemented

    def __pos__(self):
        return self

    def __neg__(self):
        return Real(-self.__value)

    def __abs__(self):
        return Real(abs(self.__value))

    def __add__(self, other):
        if isinstance(other, Real):
            return Real(self.__value + other.__value)
        elif isinstance(other, (int, float)):
            return Real(self.__value + other)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, Real):
            return Real(self.__value * other.__value)
        elif isinstance(other, (int, float)):
            return Real(self.__value * other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __imul__(self, other):
        return self * other

    # It may be better if a Fraction can be returned
    def __truediv__(self, other):
        if isinstance(other, Real):
            return Real(self.__value / other.__value)
        elif isinstance(other, (int, float)):
            return Real(self.__value / other)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        return Real(1) / self * other

    def __str__(self):
        return str(self.value)

    def factors(self) -> list[int]:
        out = []
        value = abs(self.__value)
        if isinstance(self.__value, int):
            factor = 2
            while factor <= value:
                if value % factor == 0:
                    out.append(factor)
                    value /= factor
                else:
                    factor += 1
        return out

    def sign(self) -> bool:
        if self.__value < 0:
            return False
        else:
            return True
        
    def LaTeX(self) -> str:
        return str(self.__value)

    def value(self) -> Union[int, float]:
        return self.__value

class Complex(MathiumObject):

    def __init__(self, Re, Im):
        Re = toMathium(Re)
        Im = toMathium(Im)
        if isinstance(Re, (Real, Fraction, Polynomial)):
            self.__Re = Re
        else:
            raise TypeError(f'Cannot initialize the real part with {Re.__class__}.')
        if isinstance(Im, (Real, Fraction, Polynomial)):
            self.__Im = Im
        else:
            raise TypeError(f'Cannot initialize the imaginary part with {Im.__class__}.')

    def Re(self):
        return self.__Re

    def Im(self):
        return self.__Im

    def conj(self):
        return Complex(self.__Re, -self.__Im)

    def norm(self):
        return sqrt(self.__Re * self.__Re + self.__Im * self.__Im)

    def __pos__(self):
        return self

    def __neg__(self):
        return Complex(-self.__Re, -self.__Im)

    def __add__(self, other):
        if isinstance(other, Complex):
            return Complex(self.__Re + other.__Re, self.__Im + other.__Im)
        elif isinstance(other, (int, float, Real)):
            return Complex(self.__Re + other, self.__Im)
        elif isinstance(other, (Matrix, Tensor)):
            return TypeError(f'\'+\' not supported between instances of \'Complex\' and {other.__class__}.')
        else:
            return NotImplemented
            
    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
        return self - other

    def __mul__(self, other):
        if isinstance(other, Complex):
            Re = self.__Re * other.__Re - self.__Im * other.__Im
            Im = self.__Re * other.__Im + self.__Im * other.__Re
            return Complex(Re, Im)
        elif isinstance(other, (int, float, Real)):
            Re = self.__Re * other
            Im = self.__Im * other
            return Complex(Re, Im)
        elif isinstance(other, (Matrix, Tensor)):
            return TypeError(f'\'*\' not supported between instances of \'Complex\' and {other.__class__}.')
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Complex):
            c = self * other.conj()
            r2 = other.__Re * other.__Re + other.__Im * other.__Im
            return Complex(c.__Re / r2, c.__Im / r2)
        elif isinstance(other, (int, float, Real)):
            return self * (1 / other)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, Complex):
            return other / self
        elif isinstance(other, (int, float, Real)):
            return other * (Complex(1, 0) / self)
        else:
            return NotImplemented

    def sign(self) -> bool:
        return 

    def LaTeX(self) -> str:
        outRe = ''
        outIm = ''
        if self.__Re != 0 or self.__Im == 0:
            outRe += self.__Re.LaTeX()
        if self.__Im.sign():
            outIm += '+'
        if self.__Im != 0:
            if self.__Im == -1:
                outIm += '-'
            elif self.__Im == 1:
                pass
            else:
                outIm += self.__Im.LaTeX()
            outIm += '\mathrm{i}'
        return outRe + outIm

class Fraction(MathiumObject):

    def __init__(self, numerator, denominator = 1) -> None:
        if denominator == 0 and numerator == 0:
            raise ZeroDivisionError()
        if isinstance(numerator, (int, float)):
            self.__num = Real(numerator)
        elif isinstance(numerator, (Real, Polynomial)):
            self.__num = numerator
        else:
            raise TypeError(f'Cannot initialize the numerator with {numerator.__class__}.')
        if isinstance(denominator, (int, float)):
            self.__den = Real(denominator)
        elif isinstance(denominator, (Real, Polynomial)):
            self.__den = denominator
        else:
            raise TypeError(f'Cannot initialize the denominator with {denominator.__class__}.')
        if isinstance(self.__num, Real) and isinstance(self.__den, Real):
            if isinstance(self.__num.value(), int) and isinstance(self.__den.value(), int):
                if auto_reduce:
                    self.reduce()
            elif isinstance(self.__num.value(), float) and self.__den == 1:
                while self.__num.value() % 1 not in [0, 0.0]:
                    self.__num *= 10
                    self.__den *= 10
                self.__num = Real(int(self.__num.value()))
                if auto_reduce:
                    self.reduce()

    def reduce(self):
        if isinstance(self.__num, Real) and isinstance(self.__den, Real):
            factor_num = self.__num.factors()
            factor_den = self.__den.factors()
            factor_list = self.__num.factors()
            for factor in factor_list:
                if factor in factor_den:
                    factor_num.remove(factor)
                    factor_den.remove(factor)
            num = 1 if self.__num.sign() else -1
            den = 1 if self.__den.sign() else -1
            for factor in factor_num:
                num *= factor
            for factor in factor_den:
                den *= factor
            self.__num = Real(num)
            self.__den = Real(den)
        else:
            TypeError(f'Cannot reduce {self.__num.__class__} / {self.__den.__class__} fraction.')
        return self
    
    def __pos__(self):
        return self
    
    def __neg__(self):
        if self.__den.sign():
            return Fraction(-self.__num, self.__den)
        else:
            return Fraction(self.__num, -self.__den)

    def __add__(self, other: MathiumObject):
        if isinstance(other, (int, float)):
            other = Real(other)
        if isinstance(other, (Real, Complex, Polynomial)):
            try:
                num = self.__num + self.__den * other
                den = self.__den
                return Fraction(num, den)
            except:
                return NotImplemented
        elif isinstance(other, Fraction):
            num = self.__num * other.__den + self.__den * other.__num
            den = self.__den * other.__den
            return Fraction(num, den)
        else:
            return NotImplemented

    def __radd__(self, other: MathiumObject):
        return self + other

    def __iadd__(self, other: MathiumObject):
        return self + other
    
    def __sub__(self, other: MathiumObject):
        return self + (-other)

    def __rsub__(self, other: MathiumObject):
        return -self + other

    def __isub__(self, other: MathiumObject):
        return self - other
    
    def __mul__(self, other: object):
        if isinstance(other, (int, float)):
            other = Real(other)
        if isinstance(other, (Real, Complex, Polynomial)):
            try:
                num = self.__num * other
                den = self.__den
                return Fraction(num, den)
            except:
                return NotImplemented
        elif isinstance(other, Fraction):
            num = self.__num * other.__num
            den = self.__den * other.__den
            return Fraction(num, den)
        else:
            return NotImplemented
    
    def __rmul__(self, other: MathiumObject):
        return self * other

    def __imul__(self, other: MathiumObject):
        return self * other
    
    def __truediv__(self, other: MathiumObject):
        if isinstance(other, (int, float)):
            other = Real(other)
        if isinstance(other, (Real, Complex, Polynomial)):
            try:
                other_inv = Fraction(1, other)
                return self * other_inv
            except:
                return NotImplemented
        elif isinstance(other, Fraction):
            num = self.__num * other.__den
            den = self.__den * other.__num
            return Fraction(num, den)
        else:
            return NotImplemented

    def __itruediv__(self, other: MathiumObject):
        return self / other

    def sign(self) -> bool:
        if self.__num.sign() ^ self.__den.sign():
            return False
        else:
            return True

    def LaTeX(self) -> str:
        out = ''
        if not self.sign():
            out += '-'
        out += '\\frac{' + abs(self.__num).LaTeX() + '}{' + abs(self.__den).LaTeX() + '}'
        return out

    def numerator(self) -> MathiumObject:
        return self.__num

    def denominator(self) -> MathiumObject:
        return self.__den

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

    def LaTeX(self) -> str:
        out = ''
        mtrx = self.Elem
        out += '\\begin{bmatrix}\n'
        for i in range(len(mtrx)):
            for j in range(len(mtrx[0])):
                out += toMathium(mtrx[i][j]).LaTeX()
                if j < len(mtrx[0]) - 1:
                    out += ' & '
            if i < len(mtrx) - 1:
                out += ' \\\\\n'
            else:
                out += '\n\\end{bmatrix}'
        return out

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

class Tensor(object):
    
    def __init__(self, elem: list = [0]) -> None:
        self.__deg = []
        self.__deg.append(len(elem))
        temp = elem
        while True:
            try:
                temp = temp[0]
                self.__deg.append(len(temp))
            except:
                break
        
        self.__volume = 1
        for item in self.__deg:
            self.__volume *= item

        stack = [elem]
        self.__elem: list[MathiumObject] = []
        while len(stack) != 0:
            top = stack.pop(0)
            if isinstance(top, list):
                for item in top:
                    stack.append(item)
            else:
                self.__elem.append(toMathium(top))

    def __eq__(self, other) -> bool:
        if other == 0:
            for item in self.__elem:
                if item != 0:
                    return False
            return True
        elif isinstance(other, Tensor):
            if self.__elem == other.__elem:
                return True
            else:
                return False

    def __getitem__(self, key: Union[int, tuple, list]):
        if isinstance(key, int):
            length = int(self.__volume / self.__deg[0])
            elem = self.__elem[key * length : (key + 1) * length]
            out = Tensor()
            out.__elem = elem
            out.__deg = self.__deg[1:]
            out.__volume = len(elem)
            if out.__volume == 1:
                return out.__elem[0]
            else:
                return out
        elif isinstance(key, (tuple, list)):
            if len(key) == 0:
                return self
            elif len(key) == 1:
                return self[key[0]]
            else:
                return self[key[0]][key[1:]]
        
    def __pos__(self):
        return self
    
    def __neg__(self):
        out = Tensor()
        out.__elem = []
        out.__deg = copy(self.__deg)
        out.__volume = self.__volume
        for i in range(self.__volume):
            out.__elem.append(-self.__elem[i])
        return out

    def __len__(self) -> int:
        self.__volume = 1
        for item in self.__deg:
            self.__volume *= item
        return self.__volume

    def volume(self) -> int:
        self.__volume = 1
        for item in self.__deg:
            self.__volume *= item
        return self.__volume

    def degree(self) -> list[int]:
        return copy(self.__deg)
        
    def LineToCoor(self, idx) -> list[int]:
        den = self.__volume
        out = []
        for i in range(len(self.__deg)-1):
            den = int(den / self.__deg[i])
            out.append(idx // den)
            idx %= den
        out.append(idx)
        return out

    @staticmethod
    def zeros(degree):
        out = Tensor()
        out.__deg = degree
        depth = len(degree) - 1
        out.__elem = [Real(0) for _ in range(len(out))]
        return out

    @staticmethod
    def PolyConv(ten0, ten1, map0, map1):
        degrees = max(max(map0), max(map1)) + 1
        deg = [0 for _ in range(degrees)]
        for i in range(len(map0)):
            deg[map0[i]] = ten0.__deg[i]
        for i in range(len(map1)):
            if deg[map1[i]] == 0:
                deg[map1[i]] = ten1.__deg[i]
            else:
                deg[map1[i]] += ten1.__deg[i] - 1

        volume = 1
        for item in deg:
            volume *= item
        elem = [Real(0) for _ in range(volume)]
        for i in range(len(ten0)):
            for j in range(len(ten1)):
                coor = [0 for _ in range(degrees)]
                coor0 = ten0.LineToCoor(i)
                coor1 = ten1.LineToCoor(j)
                for k in range(len(coor0)):
                    coor[map0[k]] += coor0[k]
                for k in range(len(coor1)):
                    coor[map1[k]] += coor1[k]
                idx = coor[0]
                for k in range(1, degrees):
                    idx = idx * deg[k] + coor[k]
                elem[idx] += ten0.__elem[i] * ten1.__elem[j]
        
        out = Tensor()
        out.__elem = elem
        out.__deg = deg
        out.__volume = volume
        return out

    @staticmethod
    def PolyPlus(ten0, ten1, map0, map1):
        degrees = max(max(map0), max(map1)) + 1
        deg = [0 for _ in range(degrees)]
        for i in range(len(map0)):
            deg[map0[i]] = ten0.__deg[i]
        for i in range(len(map1)):
            if deg[map1[i]] == 0:
                deg[map1[i]] = ten1.__deg[i]
            else:
                deg[map1[i]] += ten1.__deg[i] - 1
        
        volume = 1
        for item in deg:
            volume *= item
        elem = [Real(0) for _ in range(volume)]
        
        for i in range(ten0.__volume):
            coor = [0 for _ in range(degrees)]
            coor0 = ten0.LineToCoor(i)
            for k in range(len(coor0)):
                coor[map0[k]] = coor0[k]
            idx = coor[0]
            for k in range(1, degrees):
                idx = idx * deg[k] + coor[k]
            elem[idx] += ten0.__elem[i]

        for j in range(ten1.__volume):
            coor = [0 for _ in range(degrees)]
            coor1 = ten1.LineToCoor(j)
            for k in range(len(coor1)):
                coor[map1[k]] = coor1[k]
            idx = coor[0]
            for k in range(1, degrees):
                idx = idx * deg[k] + coor[k]
            elem[idx] += ten1.__elem[j]

        out = Tensor()
        out.__elem = elem
        out.__deg = deg
        out.__volume = volume
        return out

class Variable(object):

    def __init__(self, name: str) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if self.name == other.name:
            return True
        else:
            return False
    
    def __ne__(self, other: object) -> bool:
        if self.name != other.name:
            return True
        else:
            return False
        
    def __add__(self, other: object):
        if isinstance(other, int) or isinstance(other, float):
            return Polynomial([other, 1], self)
        else:
            raise TypeError()

class Polynomial(MathiumObject):

    def __init__(self, Coefficients: list, Variables: list[str] = ['x']) -> None:
        if isinstance(Variables[0], str):
            self.__Variable: list[Variable] = []
            for item in Variables:
                self.__Variable.append(Variable(item))
        elif isinstance(Variables[0], Variable):
            self.__Variable = Variables
        else:
            raise TypeError()
        if isinstance(Coefficients, list):
            self.__Coefficients = Tensor(Coefficients)
        elif isinstance(Coefficients, Tensor):
            self.__Coefficients = Coefficients
        else:
            raise TypeError()
        if self.__Coefficients == 0:
            degree = [1 for _ in range(len(Variables))]
            self.__Coefficients = Tensor.zeros(degree)
    
    def __neg__(self):
        return Polynomial(-self.__Coefficients, self.__Variable)

    def __abs__(self):
        return self

    def __add__(self, other: object):
        if isinstance(other, Polynomial):
            map0 = []
            for i in range(len(self.__Variable)):
                map0.append(i)
            map1 = []
            nOverlapIdx = len(self.__Variable)
            for i in range(len(other.__Variable)):
                try:
                    idx = self.__Variable.index(other.__Variable[i])
                    map1.append(idx)
                except:
                    map1.append(nOverlapIdx)
                    nOverlapIdx += 1
            Coefficients = Tensor.PolyPlus(self.__Coefficients, other.__Coefficients, map0, map1)
            Variables = []
            Variables += self.__Variable
            for i in range(len(map1)):
                if map1[i] >= len(self.__Variable):
                    Variables.append(other.__Variable[i])
            return Polynomial(Coefficients, Variables)
        else:
            try:
                return self + Polynomial([other], [self.__Variable[0]])
            except:
                raise TypeError()
        
    def __radd__(self, other: object):
        return self + other
    
    def __sub__(self, other: object):
        return self + (-other)
        
    def __rsub__(self, other: object):
        return -self + other
    
    def __mul__(self, other: object):
        if isinstance(other, Polynomial):
            map0 = []
            for i in range(len(self.__Variable)):
                map0.append(i)
            map1 = []
            nOverlapIdx = len(self.__Variable)
            for i in range(len(other.__Variable)):
                try:
                    idx = self.__Variable.index(other.__Variable[i])
                    map1.append(idx)
                except:
                    map1.append(nOverlapIdx)
                    nOverlapIdx += 1
            Coefficients = Tensor.PolyConv(self.__Coefficients, other.__Coefficients, map0, map1)
            Variables = []
            Variables += self.__Variable
            for i in range(len(map1)):
                if map1[i] >= len(self.__Variable):
                    Variables.append(other.__Variable[i])
            return Polynomial(Coefficients, Variables)
        else:
            try:
                return self * Polynomial([other], [self.__Variable[0]])
            except:
                raise TypeError()
        
    def __rmul__(self, other: object):
        return self * other
    
    def __truediv__(self, other: object):
        if isinstance(other, Polynomial):
            return Fraction(self, other)
        else:
            raise TypeError()

    def coefficients(self) -> Tensor:
        return self.__Coefficients

    def variables(self) -> list[Variable]:
        return copy(self.__Variable)

    def isConst(self) -> bool:
        for item in self.__Coefficients._Tensor__elem[1:]:
            if item != 0:
                return False
        return True

    def sign(self) -> bool:
        if self.isConst() or ascending:
            for item in self.__Coefficients._Tensor__elem:
                if item == 0:
                    continue
                else:
                    return item.sign()
        if not ascending:
            for item in self.__Coefficients._Tensor__elem[::-1]:
                if item == 0:
                    continue
                else:
                    return item.sign()
        return True

    def __abs__(self):
        if self.sign():
            return self
        else:
            return -self

    def LaTeX(self):
        out = ''
        coef = self.__Coefficients
        elem: list[MathiumObject] = coef._Tensor__elem
        volume = coef._Tensor__volume
        if ascending:
            for i in range(volume):
                if elem[i] == 0:
                    continue
                if elem[i].sign():
                    out += '+'
                
                if elem[i] == 1 and i != 0:
                    pass
                elif elem[i] == -1 and i != 0:
                    out += '-'
                else:
                    out += elem[i].LaTeX()
                
                powers = coef.LineToCoor(i)
                for k in range(len(powers)):
                    if powers[k] == 0:
                        continue
                    out += self.__Variable[k].name
                    if powers[k] == 1:
                        continue
                    out += '^{' + str(powers[k]) + '}'
        else:
            for i in list(range(volume))[::-1]:
                if elem[i] == 0:
                    continue
                if elem[i].sign():
                    out += '+'
                
                if elem[i] == 1 and i != 0:
                    pass
                elif elem[i] == -1 and i != 0:
                    out += '-'
                else:
                    out += elem[i].LaTeX()
                
                powers = coef.LineToCoor(i)
                for k in range(len(powers)):
                    if powers[k] == 0:
                        continue
                    out += self.__Variable[k].name
                    if powers[k] == 1:
                        continue
                    out += '^{' + str(powers[k]) + '}'
        if out == '':
            out = '0'
        if out[0] == '+':
            out = out[1:]
        return out

def toMathium(obj) -> MathiumObject:
    if isinstance(obj, (int, float)):
        return Real(obj)
    elif isinstance(obj, MathiumObject):
        return obj
    else:
        NotImplemented

def toVariable(obj) -> Variable:
    if isinstance(obj, str):
        return Variable(obj)
    elif isinstance(obj, Variable):
        return obj
    else:
        try:
            return Variable(str(obj))
        except:
            NotImplemented

def display(obj: MathiumObject):
    ipy.display(ipy.Latex('$$\n' + obj.LaTeX() + '\n$$'))

i = Complex(0, 1)
j = Complex(0, 1)

def derivative(obj: MathiumObject, var: Union[str, Variable] = 'x'):
    obj = toMathium(obj)
    var = toVariable(var)
    if isinstance(obj, Real):
        return Real(0)
    elif isinstance(obj, Complex):
        Re = derivative(obj.Re(), var)
        Im = derivative(obj.Im(), var)
        return Complex(Re, Im)
    elif isinstance(obj, Fraction):
        num = obj.numerator()
        den = obj.denominator()
        num = derivative(num, var) * den - derivative(den, var) * num
        den = den * den
        return Fraction(num, den)
    elif isinstance(obj, Polynomial):
        if var not in obj.variables():
            return Real(0)
        else:
            degree = obj.coefficients().degree()
            variables = obj.variables()
            varIndex = variables.index(var)
            if degree[varIndex] == 1:
                return Real(0)
            if degree[varIndex] > 1:
                degree[varIndex] -= 1
            coef = Tensor.zeros(degree)
            for i in range(len(coef)):
                indices = coef.LineToCoor(i)
                indices[varIndex] += 1
                coef._Tensor__elem[i] = obj.coefficients()[indices] * indices[varIndex]
            return Polynomial(coef, variables)
    else:
        NotImplemented
        
D = derivative

