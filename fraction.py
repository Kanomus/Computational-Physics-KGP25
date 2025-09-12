# Class for fractions

class frac:
    def __init__(self, num, den):
        if den == 0: raise ZeroDivisionError("Denominator given is zero")
        if num*den < 0: self.sign = -1
        else: self.sign = 1
        self.num , self.den = _reduce_frac(abs(num), abs(den))
    
    def __str__(self):
        if self.den == 1: return str(self.num)
        out = ""
        if self.sign == -1: out += "-"
        out += str(self.num)+"/"+str(self.den)
        return out
    
    def __eq__(self, other):
        if type(other) is int: other = frac(other, 1)
        return type(self)==type(other) and self.num == other.num and self.den == other.den
    
    def __add__(self, other):
        if type(other) is int:
            other = frac(other, 1)
        if type(other) is frac:
            num = (self.sign*self.num*other.den) + (other.sign*other.num*self.den)
            den = (self.den * other.den)
            return frac(*_reduce_frac(num, den))
        raise TypeError(f"Cannot add variables of type {type(self)} and {type(other)}")
    
    def __sub__(self, other):
        if type(other) is int:
            other = frac(other, 1)
        if type(other) is frac:
            num = (self.sign*self.num*other.den) - (other.sign*other.num*self.den)
            den = (self.den * other.den)
            return frac(*_reduce_frac(num, den))
        raise TypeError(f"Cannot subtract variables of type {type(self)} and {type(other)}")
    
    def __mul__(self, other):
        if type(other) is int:
            other = frac(other, 1)
        if type(other) is frac:
            num = self.num * other.num
            den = self.den * other.den
            num, den = _reduce_frac(num, den)
            return frac(self.sign*other.sign*num, den)
        raise TypeError(f"Cannot multiply variables of type {type(self)} and {type(other)}")
    
    def __truediv__(self, other):
        if type(other) is int:
            other = frac(other, 1)
        if type(other) is frac:
            num = self.num * other.den
            den = self.den * other.num
            num, den = _reduce_frac(num, den)
            return frac(self.sign*other.sign*num, den)
        raise TypeError(f"Cannot divide variables of type {type(self)} and {type(other)}")


def _lcm(a : int, b : int) -> int:
    if type(a) is not int or type(b) is not int: raise TypeError(f"Expected integer arguments, got {type(a)} and {type(b)}")
    lcm = int(a) if a>b else int(b)
    while True:
        if lcm % a == 0 and lcm % b == 0: return int(lcm)
        lcm += 1

def _hcf(a : int, b : int) -> int:
    if type(a) is not int or type(b) is not int: raise TypeError(f"Expected integer arguments, got {type(a)} and {type(b)}")
    hcf = int(a) if a<b else int(b)
    while hcf > 0:
        if a % hcf == 0 and b % hcf == 0: return int(hcf)
        hcf -= 1

def _reduce_frac(a : int, b : int):
    hcf = _hcf(a,b)
    if hcf == None: return 0,1
    return int(a/hcf), int(b/hcf)