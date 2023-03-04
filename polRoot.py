#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument('-newt', action='store_true')
exgroup.add_argument('-sec', action='store_true')
exgroup.add_argument('-hybrid', action='store_true')
parser.add_argument('initP')
parser.add_argument('initP2', nargs='?')
parser.add_argument('-maxIter', default=10000)
parser.add_argument('polyFileName')
args = parser.parse_args()

initP = float(args.initP)
if args.initP2 is None and not args.newt:
    print("Missing P2")
    exit()
initP2 = None if args.initP2 is None else float(args.initP2)
maxIter = int(args.maxIter)
if not args.polyFileName.endswith(".pol"):
    print("File name error")
    exit()

class Polynomial:
    def __init__(self, n, L):
        self.n = n
        self.L = L
    def get(self, x: int):
        sum = 0
        for i in range(len(self.L)):
            sum = sum + self.L[i] * pow(x,self.n-i)
        return sum
    def getDerivative(self, x: int):
        sum = 0
        for i in range(self.n):
            sum = sum + self.L[i] * (self.n-i) * pow(x,self.n-i-1)
        return sum

class Output:
    def __init__(self, root, iterations, outcome):
        self.root = root
        self.iterations = iterations
        self.outcome = outcome
    def to_str(self) -> str:
        return str(self.root) + " " + str(self.iterations) + " " + str(self.outcome)

file = open(args.polyFileName,"r")
n = int(file.readline())
a = [float(x) for x in file.readline().split(" ")]
p = Polynomial(n, a)
file.close()

def bisection(f, a, b, maxIter, eps):
    fa,fb,c = f.get(a),f.get(b),0
    if fa * fb >= 0:
        return Output(a if fa == 0 else b, 0, "success" if fa * fb == 0 else "fail")
    error = b - a
    for it in range(maxIter):
        error = error / 2
        c = a + error
        fc = f.get(c)
        if abs(error) < eps or fc == 0:
            return Output(c, it + 1, "success")
        if fa * fc < 0:
            b,fb = c,fc
        else:
            a,fa = c,fc
    return Output(c, maxIter, "fail")

def newton(f, x, maxIter, eps, delta):
    fx,fd = f.get(x), None
    for it in range(maxIter):
        fd = f.getDerivative(x)
        if abs(fd) <= delta:
            return Output(x, it + 1, "fail")
        d = fx / fd
        x = x - d
        fx = f.get(x)
        if abs(d) <= eps or fx == 0:
            return Output(x, it + 1, "success")
    return Output(x, maxIter, "fail")

def secant(f, a, b, maxIter, eps):
    fa,fb = f.get(a),f.get(b)
    if abs(fa) > abs(fb):
        a,b,fa,fb = b,a,fb,fa
    for it in range(maxIter):
        d = (b - a) / (fb - fa)
        b,fb = a,fa
        d = d * fa
        if abs(d) <= eps or fa == 0:
            return Output(a, it + 1, "success")
        a = a - d
        fa = f.get(a)
    return Output(a, maxIter, "fail")

def hybrid(f, a, b, maxIter, eps):
    ans = bisection(f, a, b, int(maxIter/2), eps)
    bit = ans.iterations
    ans = newton(f, ans.root, int(maxIter/2), eps, 0)
    return Output(ans.root, bit + ans.iterations, ans.outcome)

ans = Output(None, None, None)
if args.newt:
    ans = newton(p, initP, maxIter, 0, 0)
elif args.sec:
    ans = secant(p, initP, initP2, maxIter, 0)
elif args.hybrid:
    ans = hybrid(p, initP, initP2, maxIter, 0)
else:
    ans = bisection(p, initP, initP2, maxIter, 0)

file_name = args.polyFileName.replace(".pol",".sol")
file = open(file_name, "w")
file.write(ans.to_str())
file.close()
