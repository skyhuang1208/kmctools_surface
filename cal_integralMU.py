#!/usr/bin/env python3

import sys

if (len(sys.argv) != 2):
    print("\n integrate mu and output free energy\n")
    exit("use: %s [DATA INPUT FILE]\n" % sys.argv[0])

# PARAMETERS
# PARAMETERS

import scipy.integrate

x= []
y= []
i= [] # results of intergration
with open(sys.argv[1], "r") as IFILE: # start reading
    for line in IFILE:
        if line.strip()=='': continue 
        (temp, mu, avgC, errC, avgS, errS, *dump)= line.strip().split()
        x.append(float(avgC))
        y.append(float(mu))
        if len(x)>2 and 1==len(x)%2: i.append( (avgC, scipy.integrate.simps(y, x, 1, -1, 'last')) )

for a in range(len(i)): print(i[a][0], i[a][1])
