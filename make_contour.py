#!/usr/bin/env python3

import sys
import math
import decimal

if (len(sys.argv) != 5):
    print("\n Output the boundary at 0.01 to make contour map\n")
    exit("use: make_contour.py [in_file] #col(temp)# #col(xB)# #col(sro)#\n")

# PARAMETERS #
cT= int(sys.argv[2])-1
cX= int(sys.argv[3])-1
cS= int(sys.argv[4])-1

dx= 0.01 # period that output boundary
# PARAMETERS #

sback= -999
compdata= {}
with open(sys.argv[1], "r") as IFILE: # start reading
    for line in IFILE:
        [*inp]= line.strip().split()
        t= float(inp[cT])
        x= float(inp[cX])
        s= float(inp[cS])
        
        if sback != -999:
            if t != tback: exit("temperature inconsistent") 
            
            if s > sback:
                s0= sback
                s1= s
                x0= xback
                x1= x
                sout= math.ceil(sback/dx)*dx # output point
            else:
                s0= s
                s1= sback
                x0= x
                x1= xback
                sout= math.ceil(s/dx)*dx
            while(sout<s1):
                xout= x0 + (x1-x0)*(sout-s0)/(s1-s0)
                compdata.setdefault(sout, []).append(xout)
                sout += dx

        sback= s
        xback= x
        tback= t

for key, value in compdata.items():
    name= "bound_" + str("%.2f" % key)
    OFILE= open(name, "a") # open output file
    for v in value:
        print(v, t, key, file= OFILE)
print("*** Calculation's done. ***")
