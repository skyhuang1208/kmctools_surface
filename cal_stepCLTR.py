#!/usr/bin/env python3

import sys
import math
import decimal

if (len(sys.argv) != 2):
    print("\n calculate the steps that when cluster size reaches 50, 100, 500, 1000, 2000\n")
    exit("use: cal_stepCLTR.py [out.csind]\n")

goal=[50, 100, 500, 1000, 2000]

i= 0
with open(sys.argv[1], "r") as IFILE: # start reading
    for line in IFILE:
        (step, t, csize, n)= line.strip().split()
        step= int(step); t= float(t); csize= int(csize); n= int(n)
        if csize > goal[i]:
            print("GOAL: ", goal[i], "; step, time, csize: ", step, t, csize)
            i +=1
            if i>=len(goal): exit("calculation complete!")

print("calculation INcomplete: ", i, "goals are accomplished...")
