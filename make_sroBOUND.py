#!/usr/bin/env python3

import sys
from math import log

if len(sys.argv) != 3:
    print("\ninput bounadry and extract to files")
    exit("use: cal_bounadry.py data.txt [temperature]\n")

### variables ###
temp= sys.argv[2]

### reading boundary file ###
with open(sys.argv[1], "r") as IFILE: # start reading data 
    for line in IFILE:
        (conc, sro)= line.strip().split()

        with open("bound_"+sro, "a") as OFILE:
            print(conc, temp, file=OFILE)
