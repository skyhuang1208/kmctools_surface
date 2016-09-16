#!/usr/bin/env python3

from sys import argv
from os.path import isfile

if len(argv) != 3:
    print("add energy (out.energy) to chemical potential (mu*(nA-nB))")
    exit("usage: add.py [LOG file] [OUT.ENERGY file]")

# read and calculate
step= []
na= []
nb= []
energy= []
with open(argv[1], "r") as IFILE: # Reading his file
    for line in IFILE:
        (a, b, c, d, e, *f)= line.split()
        step.append(a)
        na.append(d)
        nb.append(e)

with open(argv[2], "r") as IFILE2:
    for line in IFILE2:
        (g, h, i)= line.split()
        energy.append(i)

for i in range(len(step)):
    print(int(step[i]), float(energy[i])+0.7*(float(na[i])-float(nb[i])))
