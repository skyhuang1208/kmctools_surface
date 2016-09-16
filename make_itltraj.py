#!/usr/bin/env python3

import sys
import os.path

if len(sys.argv) != 4:
    print("Make trajactory movie")
    exit("Usage: %s <t0.xyz> <history.def> N_steps" % sys.argv[0])

# global variables
iltcp= []
x= []
y= []
z= []
# global variables

def output():
    print(len(iltcp))
    print("N_steps=", sys.argv[3])
    for l in iltcp:
        k= int(l[1]) # index
        print(l[0], x[k], y[k], z[k])

with open(sys.argv[1], 'r') as T0: # Reading t0 file
    i= 0
    for line in T0:
        i +=1
        if i==1: N_atoms= int(line)
        elif i==2: pass
        else:
            (tp, xi, yi, zi, *others) = line.split()
            x.append(xi)
            y.append(yi)
            z.append(zi)
            if 2==int(tp) or 3==int(tp):
                iltcp.append( (tp, i-3) )

if(N_atoms != len(x)): exit("number inconsistent (N_atoms, len_X) %d %d" % (N_atoms, len(x)))

with open(sys.argv[2], 'r') as IN: # Read his.def
    i= 0
    output()
    for line in IN:
        i +=1
        if 1==i%3:
            N= int(line)
            if(N != 1): exit("suppose to be 1 defect only but: %d" % N)
        if 2==i%3:
            (T, step, time_) = line.split()
            if(int(step) != (i+1)/3): exit("step inconsistent: %d %d" % (int(step), (i+1)/3))
        if 0==i%3:
            (tp, ltcp, *others) = line.split()
            iltcp.append( (tp, ltcp) )
            if 0==(i/3)%int(sys.argv[3]): output()

