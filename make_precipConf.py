#!/usr/bin/env python3

#####      make calibration sro files (sep, dimers, trimers)   #####
#####                                      By Sky on Mar. 2016 #####

import sys
import numpy

if len(sys.argv) != 5:
    print("Make a totally aggregate .sol file")
    exit("Usage: %s nx conc_B r0 dr" % sys.argv[0])

### parameters ###
nx= sys.argv[1]
nb= sys.argv[2]*nx*nx*nx
r0= sys.argv[3]
dr= sys.argv[4]
### parameters ###

def cal_dis(da, db, dc):
    if (da>(nx/2-1)) or (db>(nx/2-1)) or (dc>(nx/2-1)): exit("(cal_dis) the dis is larger than half box length")
    
    vbra=  [[-0.5,  0.5,  0.5],
            [ 0.5, -0.5,  0.5], 
            [ 0.5,  0.5, -0.5]]
    
    dx= da*vbra[0][0] + db*vbra[1][0] + dc*vbra[2][0]
    dy= da*vbra[0][1] + db*vbra[1][1] + dc*vbra[2][1]
    dz= da*vbra[0][2] + db*vbra[1][2] + dc*vbra[2][2]

    return numpy.sqrt(dx*dx + dy*dy + dz*dz)

i_range= 0
while n<nb:
    r= r0 + i_range*dr
    n=0
    states[i][j][k]= [[[0 for k in range(nx)] for j in range(nx)] for i in range(nx)]
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                if cal_dis(i-nx/2, j-nx/2, k-nx/2) <= r:
                    n += 1
                    states[i][j][k]= -1
    if i_range==0 and n>nb: exit("use larger r0: (n nb)", n, nb)
    print("time", i_range, ": n=", n, ", nb=", nb, ", pct=", float(n)/nb)
    i_range += 1

with open("temp.sol", "r") as OFILE:
    print(n, file= OFILE)
    print("T: 1 1", file= OFILE)
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                if states[i][j][k]== -1:
                    ltcp= i*nx*nx + j*nx + k
                    print(ltcp, file= OFILE)
