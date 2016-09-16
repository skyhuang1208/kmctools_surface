#!/usr/bin/env python3

def func_fitted(x, a, b, c):        # fitting curves
    return x*(1-x)*(a+b*x+c*x**2)
def fitting(x, y):
    import scipy.optimize as opt
    return opt.curve_fit(func_fitted, x, y)[0]

import sys

if len(sys.argv) != 2:
    print("Fitting data to x(1-x)f(x) where we'll get a")
    exit("Usage: %s <input_file1> <input_file2>" % sys.argv[0])

x= []
y= []
with open(sys.argv[1], "r") as IFILE: # start reading data 
    for line in IFILE:
        (c, e)= line.strip().split()
        x.append(float(c))
        y.append(float(e))

    print("omega_s= : ", fitting(x, y))
