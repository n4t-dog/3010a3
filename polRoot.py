#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument('-newt', action='store_true')
exgroup.add_argument('-sec', action='store_true')
exgroup.add_argument('-hybrid', action='store_true')
parser.add_argument('initP')
exgroup.add_argument('initP2', nargs='?')
parser.add_argument('-maxIter', default=10000)
parser.add_argument('polyFileName')
args = parser.parse_args()

if args.newt: print('newt')
elif args.sec: print('sec')
elif args.hybrid: print('hybrid')
elif args.initP2 is None:
    print("OH FUCK")
    exit()
print(args.initP, args.initP2)
print("maxIter:", args.maxIter)
print("polyFileName:", args.polyFileName)
if not args.polyFileName.endswith(".pol"):
    print("OH FUCK")
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

file = open(args.polyFileName,"r")
n = int(file.readline())
a = [float(x) for x in file.readline().split(" ")]
f = Polynomial(n, a)
file.close()

def bisection(a, b, maxIter, eps):
    fa,fb = f.get(a),f.get(b)

    if fa * fb >= 0:
        print("Inadequate values for a and b.")
        return None
    error = b - a
    for it in range(maxIter):
        error = error / 2
        c = a + error
        fc = f.get(c)
        if abs(error) < eps or fc == 0:
            print("Algorithm has converged after",it,"iterations!")
            return c
        if fa * fc < 0:
            b,fb = c,fc
        else:
            a,fa = c,fc
    print("Max iterations reached without convergence...")
    return None

b = bisection(int(args.initP), int(args.initP2), args.maxIter, 0)
print("root:",b)
