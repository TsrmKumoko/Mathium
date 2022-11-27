class Fraction(object):

    def __init__(self, num, den) -> None:
        self.num = num
        self.den = den
        self.sign = ''
        if isinstance(num, int) and isinstance(den, int):
            if abs(num) / num * abs(den) / den < 0:
                self.sign = '-'
        elif isinstance(num, int):
            if abs(num) / num < 0:
                self.sign = '-'
        elif isinstance(den, int):
            if abs(den) / den < 0:
                self.sign = '-'

    def reduce(self):
        if isinstance(self.num, int) and isinstance(self.den, int):
            num, den = self.num, self.den
            while num % den != 0:
                mod = num % den
                num, den = den, mod
            self.num /= mod
            self.den /= mod
        else:
            TypeError(f'Cannot reduce {self.num.__class__} / {self.den.__class__} fraction.')
        return self
    
    def __neg__(self):
        self.num = -self.num
        num, den = self.num, self.den
        if isinstance(num, int) and isinstance(den, int):
            if abs(num) / num * abs(den) / den < 0:
                self.sign = '-'
        elif isinstance(num, int):
            if abs(num) / num < 0:
                self.sign = '-'
        return self

    def __add__(self, __o: object):
        if isinstance(__o, float):
            try:
                return self.num / self.den + float
            except:
                return float('NaN')
        elif isinstance(__o, Fraction):
            num = self.num * __o.den + self.den * __o.num
            den = self.den * __o.den
            return Fraction(num, den)
        else:
            try:
                num = self.num + self.den * __o
                den = self.den
                return Fraction(num, den)
            except:
                TypeError()

    def __radd__(self, __o: object):
        return self + __o
    
    def __sub__(self, __o: object):
        return self + (-__o)
    
    def __mul__(self, __o: object):
        if isinstance(__o, float):
            try:
                return self.num / self.den * float
            except:
                return float('NaN')
        elif isinstance(__o, Fraction):
            num = self.num * __o.num
            den = self.den * __o.den
            return Fraction(num, den)
        else:
            try:
                num = self.num * __o
                den = self.den
                return Fraction(num, den)
            except:
                TypeError()
    
    def __rmul__(self, __o: object):
        return self * __o
