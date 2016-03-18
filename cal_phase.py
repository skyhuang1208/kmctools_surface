#!/usr/bin/env python3

##### determine phase by short range order (sro) data 
##### criterion: 1. if sro is greater (less) than upper (lower) 
#####               -> separate (order)                 
#####            2. if final value, 1st derivative, 2nd derivative are all + (-) 
#####               -> separate (order)
#####            3. others are random solution
#####                               By Sky on March 2016 

import sys

if (len(sys.argv) != 7):
    print("\nread out.sro and determines phase (order, sep, or random)")
    exit("use: cal_phase.py [IN_SRO] [OUT_sep] [OUT_ord] [OUT_ran] sol% temp\n")

# PARAMETERS
upper=  0.15 
lower= -0.15
# PARAMETERS

x_last=0 # last value of SRO
with open(sys.argv[1], "r") as IFILE: # start reading and change parameters
    for line in IFILE:
        (step, t, sro)= line.strip().split()
        x_last= sro

name_out= str()
if   x_last >= upper:
    name_out= sys.argv[2]
elif x_last <= lower:
    name_out= sys.argv[3]
else:
    name_out= sys.argv[4]

OFILE= open("TEMP_"+sys.argv[1], "a") # open output file
print(sys.argv[5], sys.argv[6], file=OFILE)

print("*** Calculation completed! ***")
