#!/usr/bin/env python3

import sys
import math
import decimal

if (len(sys.argv) != 2):
    print("\n duplicate the system to become 2x2x2\n")
    exit("use: make_duplicate.py [in_file]\n")

nx= 64
ny= 64
nz= 64

with open(sys.argv[1], "r") as IFILE: # start reading
    l= 0
    for line in IFILE:
        l +=1
        if l==1: 
            nline= int(line)
            print(nline*8)
        elif l==2:
            print(line, end='')
        else:
            (t, i, j, k, other)= line.strip().split(None, 4)
            for a in range(2):
                for b in range(2):
                    for c in range(2):
                        x= int(i)+a*nx
                        y= int(j)+b*ny
                        z= int(k)+c*nz

                        print(t, x, y, z, other)

