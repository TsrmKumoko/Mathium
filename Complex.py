from math import sqrt

class Complex:

    def __init__(self, Re, Im):
        self.Re = Re
        self.Im = Im

    def Re(self):
        return self.Re

    def Im(self):
        return self.Im

    def conj(self):
        return Complex(self.Re, -self.Im)

    def norm(self):
        return sqrt(self.Re ** 2 + self.Im ** 2)

    def __pos__(self):
        return self

    def __neg__(self):
        return Complex(-self.Re, -self.Im)

    def __add__(self, other):
        if isinstance(other, Complex):
            return Complex(self.Re + other.Re, self.Im + other.Im)
        elif isinstance(other, int) or isinstance(other, float):
            return Complex(self.Re + other, self.Im)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        Ans = self + other
        return Ans

    def __sub__(self, other):
        if isinstance(other, Complex):
            return Complex(self.Re - other.Re, self.Im - other.Im)
        elif isinstance(other, int) or isinstance(other, float):
            return Complex(self.Re - other, self.Im)
        else:
            return NotImplemented

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
        return self - other

    def __mul__(self, other):
        if isinstance(other, Complex):
            Re = self.Re * other.Re - self.Im * other.Im
            Im = self.Re * other.Im + self.Im * other.Re
            return Complex(Re, Im)
        elif isinstance(other, int) or isinstance(other, float):
            Re = self.Re * other
            Im = self.Im * other
            return Complex(Re, Im)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Complex):
            c = self * other.conj()
            r2 = other.Re ** 2 + other.Im ** 2
            return Complex(c.Re / r2, c.Im / r2)
        elif isinstance(other, int) or isinstance(other, float):
            return self * (1 / other)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, Complex):
            return other / self
        elif isinstance(other, int) or isinstance(other, float):
            return other * (Complex(1, 0) / self)
        else:
            return NotImplemented

def Re(c):
    return c.Re

def Im(c):
    return c.Im

def conj(c):
    return Complex(c.Re, -c.Im)

def norm(c):
    return sqrt(c.Re ** 2 + c.Im ** 2)

i = Complex(0, 1)
j = Complex(0, 1)
