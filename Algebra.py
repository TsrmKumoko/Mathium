from copy import copy, deepcopy
from Fraction import *

class Variable(object):

    def __init__(self, name) -> None:
        self.name = name

    def __eq__(self, __o: object) -> bool:
        if self.name == __o.name:
            return True
        else:
            return False
    
    def __ne__(self, __o: object) -> bool:
        if self.name != __o.name:
            return True
        else:
            return False
        
    def __add__(self, __o: object):
        if isinstance(__o, int) or isinstance(__o, float):
            return Polynomial([__o, 1], self)
        else:
            raise TypeError()
        
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
        self.__elem = []
        while len(stack) != 0:
            top = stack.pop(0)
            if isinstance(top, list):
                for item in top:
                    stack.append(item)
            else:
                self.__elem.append(top)

    def __getitem__(self, key):
        length = int(self.__volume / self.__deg[0])
        elem = self.__elem[key * length : (key + 1) * length]
        out = Tensor(elem, self.__deg[1:])
        if out.__volume == 1:
            return out.__elem[0]
        else:
            return out
        
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
        elem = [0 for _ in range(volume)]
        for i in range(ten0.__volume):
            for j in range(ten1.__volume):
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
                elem[idx] = ten0.__elem[i] * ten1.__elem[j]
        
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
        elem = [0 for _ in range(volume)]
        
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

class Polynomial(object):

    def __init__(self, Coefficients: list, Variables: list[str] = ['x']) -> None:
        if isinstance(Variables[0], str):
            self.__Variable = []
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
    
    def __neg__(self):
        return Polynomial(-self.__Coefficients, self.__Variable)

    def __add__(self, __o: object):
        if isinstance(__o, Polynomial):
            map0 = []
            for i in range(len(self.__Variable)):
                map0.append(i)
            map1 = []
            nOverlapIdx = len(self.__Variable)
            for i in range(len(__o.__Variable)):
                try:
                    idx = self.__Variable.index(__o.__Variable[i])
                    map1.append(idx)
                except:
                    map1.append(nOverlapIdx)
                    nOverlapIdx += 1
            Coefficients = Tensor.PolyPlus(self.__Coefficients, __o.__Coefficients, map0, map1)
            Variables = []
            Variables += self.__Variable
            for i in range(len(map1)):
                if map1[i] >= len(self.__Variable):
                    Variables.append(__o.__Variable[i])
            return Polynomial(Coefficients, Variables)
        else:
            try:
                return self + Polynomial([__o], [self.__Variable[0]])
            except:
                raise TypeError()
        
    def __radd__(self, __o: object):
        return self + __o
    
    def __sub__(self, __o: object):
        return self + (-__o)
        
    def __rsub__(self, __o: object):
        return -self + __o
    
    def __mul__(self, __o: object):
        if isinstance(__o, Polynomial):
            map0 = []
            for i in range(len(self.__Variable)):
                map0.append(i)
            map1 = []
            nOverlapIdx = len(self.__Variable)
            for i in range(len(__o.__Variable)):
                try:
                    idx = self.__Variable.index(__o.__Variable[i])
                    map1.append(idx)
                except:
                    map1.append(nOverlapIdx)
                    nOverlapIdx += 1
            Coefficients = Tensor.PolyConv(self.__Coefficients, __o.__Coefficients, map0, map1)
            Variables = []
            Variables += self.__Variable
            for i in range(len(map1)):
                if map1[i] >= len(self.__Variable):
                    Variables.append(__o.__Variable[i])
            return Polynomial(Coefficients, Variables)
        else:
            try:
                return self * Polynomial([__o], [self.__Variable[0]])
            except:
                raise TypeError()
        
    def __rmul__(self, __o: object):
        return self * __o
    
    def __truediv__(self, __o: object):
        if isinstance(__o, Polynomial):
            return Fraction(self, __o)
        else:
            raise TypeError()
