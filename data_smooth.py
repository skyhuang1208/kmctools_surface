#!/usr/bin/env python3

#********* data analysis functions *********#
# transform data to a smooth curve
# 1. fitting: use least-square curve fitting
# 2. smooth_SG: smooth curve by locally fitting
# 3. smooth_AA: sum[i-hw:i+hw] or running avg

def func_fitted(x, a, b, c):        # fitting curves
    return a + b * x**c
def fitting(x, y):
    import scipy.optimize as opt
    return opt.curve_fit(func_fitted, x, y, guess0)[0]

def smooth_SG(x, pwin=41, pord=2):  # smooth data by Savitzky-Golay
    from scipy.signal import savgol_filter
    return savgol_filter(x, pwin, pord)

def smooth_AA(x, pwin=41):          # smooth data by adjacent averaging
    import numpy
    hw= int(pwin/2.0) # half window
    smoothed= []
    for i in range(len(x)):
        smoothed.extend(numpy.mean(x[max(0, i-hw):min(len(x), i+hw+1)]))
    return smoothed
#####***** data analysis functions *****#####

from sys import argv
from os.path import isfile

if len(argv) != 2 and len(argv) != 3:
    print("\nSmooth curve by Savitzky-Golay or adjacent averaging")
    exit("Usage: %s <input_file> name(optional)\n" % argv[0])

# input parameters
N_col= input("Insert number of columns\n")
data= []
for i in range(int(N_col)): data.append([])

in_isc= input("Do you need chop vacuum data (rm zeros)? (yes/no)\n")
if in_isc=="yes": is_chop= True
else:             is_chop= False

in_pwin= input("Insert number of points of window (dafault:41)\n")
if in_pwin.isdigit(): N_pwin= int(in_pwin)
else:                 N_pwin= 41

# open output file
name_out= "out"
if len(argv)==3: name_out= name_out + argv[2]
OFILE= open(name_out, "w")

# read
if not isfile(argv[1]): exit("%s doesn\'t exist!" % argv[1])
with open(argv[1], "r") as IFILE: # Reading his file
    for line in IFILE:
        indata = line.strip().split()
        for i in range(len(data)):
            data[i].append(float(indata[i]))

#chop zero
TOL_BOUND= 0.05 # smaller than this considerd vacuum
left= -1; right= len(data[-1]) 
for i in range(len(data[-1])):
    if data[-1][i]>TOL_BOUND: left = i; break
for i in range(len(data[-1])-1, -1, -1):
    if data[-1][i]>TOL_BOUND: right= i; break

# calculate
if is_chop:
    dataSM= [0]*left
    dataSM.extend(smooth_SG(data[-1][left:right+1], N_pwin))
    dataSM.extend([0]*(len(data[-1])-right-1))
else:
    dataSM= smooth_SG(data[-1], N_pwin)

if len(dataSM) != len(data[0]): print("wrong dataSM: ", len(dataSM), len(data[0]))
# output
for i in range(len(data[0])):
    for j in range(len(data)-1):
        print(data[j][i], end=" ", file= OFILE)
    print(dataSM[i], file= OFILE)

