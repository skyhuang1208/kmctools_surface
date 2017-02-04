#!/usr/bin/env python3

import sys
import math
import decimal

if (len(sys.argv) != 5):
    print("\n Output the SRO values for every 0.01 in xB\n")
    exit("use: make_contour.py [in_file] #col(temp)# #col(xB)# #col(sro)#\n")

# PARAMETERS #
cT= int(sys.argv[2])-1
cX= int(sys.argv[3])-1
cS= int(sys.argv[4])-1

dx= 0.01 # period that output boundary
xSTOP= 0.4 # where last x is calculated
# PARAMETERS #

xOUT= 0.00 # the next time output values
with open(sys.argv[1], "r") as IFILE: # start reading
    i= 0
    for line in IFILE:
        if not line.strip(): continue
        i +=1

        [*inp]= line.strip().split()
        t= float(inp[cT])
        x= float(inp[cX])
        s= float(inp[cS])
        
        if i !=1:
            if t != tback: exit("temperature inconsistent")

            while i==2 and xOUT < xback:
                print(xOUT, t, s-(s-sback)*(x-xOUT)/(x-xback))
                xOUT += dx

            while xOUT > xback and xOUT < x and xOUT < (xSTOP+dx/2):
                print(xOUT, t, (s-sback)*(xOUT-xback)/(x-xback)+sback)
                xOUT += dx

        sback= s
        xback= x
        tback= t

print("") # add a blank line
