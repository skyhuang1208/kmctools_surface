#!/usr/bin/env python3

import sys

if (len(sys.argv) != 8):
    print("\n calculate average value of concentration and SRO\n")
    exit("use: cal_avg.py [LOG_file] [SRO_file] [eSGC_file] #N_data_point# [OUT_file] #temperature# #mu#\n")

# PARAMETERS
N_pts= int(sys.argv[4])
N_col= 2 # column of data for out.sro
# PARAMETERS

data= []
iscal= 0
with open(sys.argv[1], "r") as IFILE: # start reading
    for line in IFILE:
        if line.strip()=='': continue 
        [*inp]= line.strip().split()
        
        if inp[0]=='0': iscal= 1
        if inp[0]=='****': iscal= 0
        
        if 1==iscal: data.append( (float(inp[3]), float(inp[4])) )

data2= []
with open(sys.argv[2], "r") as IFILE2: # start reading
    for line in IFILE2:
        if line.strip()=='': continue 
        [*inp]= line.strip().split()
        
        data2.append(float(inp[N_col]))

if   len(data2) <= N_pts: exit("error: Total data number is less than N of rows: %d %d" % (len(data2), N_pts)) # check data N
elif len(data2) <= 2*N_pts: print("Warning: Total data number is small: (N_data) (N_required)", len(data2), N_pts)

ISup= 0 # if energy goes up
with open(sys.argv[3], "r") as IFILE3: # start reading
    data3= []
    for line in IFILE3:
        if line.strip()=='': continue 
        [*inp]= line.strip().split()
        
        data3.append(float(inp[N_col]))
    if data3[0]<data3[-1]: ISup=  1
    else:                  ISup= -1

comp= []
sro= []
for i in range(-1, -N_pts-1, -1):
    comp.append( data[i][1]/(data[i][0]+data[i][1]) )
    sro.append( data2[i] )

import statistics as stat

OFILE= open(sys.argv[5], "a") # open output file
print(sys.argv[6], sys.argv[7], stat.mean(comp), stat.stdev(comp), stat.mean(sro), stat.stdev(sro), ISup, file= OFILE)
print("*** Calculation's done. ***")
