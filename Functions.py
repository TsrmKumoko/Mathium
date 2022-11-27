from math import cos, sin, sqrt, atan
import IPython.display as ipy
from Complex import *
from Matrix import *
from Fraction import *
from Algebra import *

ascending = False

def LaTeX(obj) -> str:
    out = ''
    if isinstance(obj, Polynomial):
        coef = obj._Polynomial__Coefficients
        elem = coef._Tensor__elem
        volume = coef._Tensor__volume
        if ascending:
            for i in range(volume):
                if elem[i] == 0:
                    continue
                if elem[i] > 0:
                    out += '+'
                
                if elem[i] == 1 and i != 0:
                    pass
                elif elem[i] == -1 and i != 0:
                    out += '-'
                else:
                    out += LaTeX(elem[i])
                
                powers = coef.LineToCoor(i)
                for k in range(len(powers)):
                    if powers[k] == 0:
                        continue
                    out += obj._Polynomial__Variable[k].name
                    if powers[k] == 1:
                        continue
                    out += '^{' + LaTeX(powers[k]) + '}'
        else:
            for i in list(range(volume))[::-1]:
                if elem[i] == 0:
                    continue
                if elem[i] > 0:
                    out += '+'
                
                if elem[i] == 1 and i != 0:
                    pass
                elif elem[i] == -1 and i != 0:
                    out += '-'
                else:
                    out += LaTeX(elem[i])
                
                powers = coef.LineToCoor(i)
                for k in range(len(powers)):
                    if powers[k] == 0:
                        continue
                    out += obj._Polynomial__Variable[k].name
                    if powers[k] == 1:
                        continue
                    out += '^{' + LaTeX(powers[k]) + '}'
        if out[0] == '+':
            out = out[1:]
    elif isinstance(obj, Matrix):
        mtrx = obj.Elem
        out += '\\begin{bmatrix}\n'
        for i in range(len(mtrx)):
            for j in range(len(mtrx[0])):
                out += LaTeX(mtrx[i][j])
                if j < len(mtrx[0]) - 1:
                    out += ' & '
            if i < len(mtrx) - 1:
                out += ' \\\\\n'
            else:
                out += '\n\\end{bmatrix}'
    elif isinstance(obj, Complex):
        c = obj
        outRe = ''
        outIm = ''
        if c.Re != 0 or c.Im == 0:
            outRe += LaTeX(c.Re)
        if c.Im > 0:
            outIm += '+'
        if c.Im != 0:
            if c.Im == -1:
                outIm += '-'
            elif c.Im == 1:
                pass
            else:
                outIm += LaTeX(c.Im)
            outIm += '\mathrm{i}'
        out = outRe + outIm
    elif isinstance(obj, Fraction):
        frac = obj
        out += frac.sign
        num, den = frac.num, frac.den
        if isinstance(num, int):
            num = abs(num)
        if isinstance(den, int):
            den = abs(den)
        out += '\\frac{' + LaTeX(num) + '}{' + LaTeX(den) + '}'
    else:
        out = str(obj)
    return out
        
def display(obj):
    ipy.display(ipy.Latex('$$\n' + LaTeX(obj) + '\n$$'))

def Pol2Rec(r, theta):
    x = r * cos(theta)
    y = r * sin(theta)
    return x, y

def Rec2Pol(x, y):
    r = sqrt(x * x + y * y)
    theta = atan(y / x)
    return r, theta

def prll(*args):
    n = len(args)
    Ans = 0
    for i in range(n):
        Ans += 1 / args[i]
    return 1 / Ans

def fac(n: int) -> int:
    if n == 0:
        return 1
    elif n > 0:
        return n * fac(n - 1)
    else:
        return -1

def arr(m: int, n: int) -> int:
    if n >= 0 and m >= 0 and m >= n:
        return fac(m) / fac(m - n)
    else:
        return -1

def com(m: int, n: int) -> int:
    if n >= 0 and m >= 0 and m >= n:
        return fac(m) / fac(m - n) / fac(n)
    else:
        return -1
    
