#!/usr/bin/env python3

def func_fitted(x, a, b):        # fitting curves
    return a+b*x
def fitting(x, y):
    import scipy.optimize as opt
    return opt.curve_fit(func_fitted, x, y)[0]

import sys

if len(sys.argv) != 2:
    print("calculate diffusivity with the input of out.msd")
    exit("Usage: %s <input_file>" % sys.argv[0])

x= []
y= []
with open(sys.argv[1], "r") as IFILE: # start reading data 
    for line in IFILE:
        (step, t, msd)= line.strip().split()
        x.append(float(t))
        y.append(float(msd))

    print("diffusivity: ", fitting(x, y)[1])
