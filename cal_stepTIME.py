#!/usr/bin/env python3

import sys
import math
import decimal

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print("\nCalculate steps of 0, 5, 10, 15, 20, 25, 30, 35 (s) after incubation")
    print("incubation defines by the time cltr reaches 50")
    exit("use: cal_stepTIME.py [out.csind] is_correct80(optional)\n")

goal=[5, 10, 15, 20, 25, 30, 35]
cltrINCUBATION= 50

i= 0
INCtime= 999999999999
with open(sys.argv[1], "r") as IFILE: # start reading
    for line in IFILE:
        (step, t, csize, n)= line.strip().split()
        step= int(step)
        t= float(t)
        csize= int(csize)
        n= int(n)

        if len(sys.argv)==3: t= t * (80*80*80) / (64*64*64) # apply time correction for 80 system box
        
        if csize > cltrINCUBATION and INCtime==999999999999: 
            print("INCUBATION STEP, TIME, CSIZE: ", step, t, csize)
            INCtime= t

        if csize>cltrINCUBATION and t > (INCtime+goal[i]):
            print(goal[i], "(s) after incubation STEP, TIME, CSIZE: ", step, t, csize)
            if i== (len(goal)-1): exit("DONE!")
            i +=1

print("calculation INcomplete: ", i, "goals are accomplished...")
