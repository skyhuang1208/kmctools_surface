#!/usr/bin/env python3

###### read kmc_par.h and change value of some parameters  #######

import sys

if (len(sys.argv) < 4) or (len(sys.argv)%2 != 0):
    print("\nread kmc_par.h and change some parameters")
    exit("use: auto_chpar.py kmc_par.h [variable] [value] [variable] [value] ...\n")


var= [] # targeted variables
val= [] # values of the vars
for n in range(1, int(len(sys.argv)/2)):
    var.append(sys.argv[n*2])
    val.append(sys.argv[n*2+1])

OFILE= open("TEMP_"+sys.argv[1], "w") # open output file

with open(sys.argv[1], "r") as IFILE: # start reading and change parameters
    for line in IFILE:
        isch= 0
        if "=" in line:
            (w0, w1)= line.strip().split("=")
            for n in range(len(var)):
                if var[n] in w0:
                    print(w0, "=", val[n], ";", file=OFILE)
                    isch= 1 # the value of the var is changed
                    break
        if 0==isch:
            print(line, file=OFILE, end="")

print("*** Calculation completed! ***")
