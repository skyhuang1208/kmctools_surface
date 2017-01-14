#!/usr/bin/env python3

import sys
import os.path

if len(sys.argv) != 2:
    print("Input history.dis and output time evolution of r^2 and r1*r2")
    exit("Usage: %s <history.dis>" % sys.argv[0])

# global variables
vbra= ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5));
iltcp= []
x= []
y= []
z= []
# global variables

def cal_dis2(i, j, k):
    x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
    y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
    z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]

    return x*x + y*y + z*z

def cal_dot(i1, j1, k1, i2, j2, k2):
    x1= i1*vbra[0][0] + j1*vbra[1][0] + k1*vbra[2][0]
    y1= i1*vbra[0][1] + j1*vbra[1][1] + k1*vbra[2][1]
    z1= i1*vbra[0][2] + j1*vbra[1][2] + k1*vbra[2][2]
    x2= i2*vbra[0][0] + j2*vbra[1][0] + k2*vbra[2][0]
    y2= i2*vbra[0][1] + j2*vbra[1][1] + k2*vbra[2][1]
    z2= i2*vbra[0][2] + j2*vbra[1][2] + k2*vbra[2][2]

    return x1*x2 + y1*y2 + z1*z2

vcc= []
sol= [0, 0, 0]
with open(sys.argv[1], 'r') as IN: # Read his.dis
    l= 0; N= 99999; 
    dis2D= 0; dis2S= 0; disDOT= 0 
    for line in IN:
        l +=1
        
        if   l==1: N= int(line)
        elif l==2: (T, step, time1_, time2_) = line.split()
        else:
            (tp, x, y, z) = line.split()
            if l==3 and int(tp)==-1: exit("Error: first atom is not a defect")
            if l >3 and int(tp)!=-1: exit("Error: a not-first atom is not a solute")

            x= int(x)
            y= int(y)
            z= int(z)
            if l==3:
                dis2D += cal_dis2(x, y, z)
                vcc= [x, y, z]
            else:
                dis2S  += cal_dis2(x, y, z)
                sol[0] += x; sol[1] += y; sol[2] += z

        if (l-2)==N:
            disSV= cal_dot(sol[0], sol[1], sol[2], vcc[0], vcc[1], vcc[2])
            disSS= cal_dot(sol[0], sol[1], sol[2], sol[0], sol[1], sol[2])
            if N==1: print(step, float(time1_), float(time2_), dis2D, "END")
            else:    print(step, float(time1_), float(time2_), dis2D, dis2S/(N-1), disSV, disSS)
            
            l= 0; dis2D= 0; dis2S= 0; disDOT= 0; vcc= []; sol= [0, 0, 0]

