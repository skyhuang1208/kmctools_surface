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


#********* constants *********#
nx= 660
ny= 270 
nz=   4
nMl= 20
ntotal= nx*ny*nz
SEG_POSI= 150.0 # minimum position of segregation
SEG_FRAC=   0.1 # deviation fraction of segregation
TOL_BOUND= 0.05 # smaller than this considerd vacuum
#########* constants *#########

#********* variables *********#
step= 0
time= 0.0
dis= 0.5*(2**0.5)
pos= [(tempx-nx/2.0)*dis for tempx in range(0, nx)]
#########* variables *#########

#********* calculate degree of seg *********#
def find_edge(x, y):
    if len(x)!=nx or len(y)!=nx: exit("error(find_edge): x(%d) or y(%d) != nx" % (len(x), len(y)))
    
    y_mid= y[int(nx/2.0)]
    upper= y_mid*(1+SEG_FRAC)
#    lower= y_mid*(1-SEG_FRAC)
    lower= -1 
    ledge= -1; redge= nx
    for i in range(int(nx/2.0), -1, -1):
        if x[i]<-SEG_POSI and (y[i]>upper or y[i]<lower): ledge= i; break
    for i in range(int(nx/2.0), nx,  1):
        if x[i]> SEG_POSI and (y[i]>upper or y[i]<lower): redge= i; break
    return (ledge, redge)

def cal_seg(y, lxlo, lxhi, rxlo, rxhi, c_avg):
    lslen= lxhi-lxlo+1 # left  seg length
    rslen= rxhi-rxlo+1 # right seg length
    dseg= 0
    if lslen >0: dseg += sum(y[lxlo:lxhi+1]) - c_avg*lslen
    if rslen >0: dseg += sum(y[rxlo:rxhi+1]) - c_avg*rslen
    return dseg*dis
#####***** calculate degree of seg *****#####

from sys import argv
from os.path import isfile

if len(argv) != 2 and len(argv) != 3:
    print("calculate degree of segregation")
    exit("Usage: %s <history.sol> name(optional)" % argv[0])

# open output file
name_out= "out"
if len(argv)==3: name_out= name_out + argv[2]
OFILE= open(name_out, "w")

# read and calculate
if not isfile(argv[1]): exit("%s doesn\'t exist!" % argv[1])
with open(argv[1], "r") as IFILE: # Reading his file
    while True:
        # first line
        try: nlines= int(IFILE.readline())
        except: break # check EOF

        # second line
        (dump, step, time)= IFILE.readline().split()
    
        # reading the block
        den= [0]*nx
        for i in range(0, nlines):
            ltcp= int(IFILE.readline())
            x= int(ltcp/nz/ny)
            den[x]= den[x] + 1.0/ny/nz
        import numpy
        c_bar= numpy.mean(den)*nx/(nx-2.0*nMl)

        # computing segregation length
        lbound= 0; rbound= 0
        for x in range(0, nx): # 
            if den[x]>TOL_BOUND: lbound= x; break
        for x in range(nx-1, -1, -1):
            if den[x]>TOL_BOUND: rbound= x; break
        
        denSM= [0]*lbound
        denSM.extend(smooth_SG(den[lbound:rbound+1]))
        denSM.extend([0]*(nx-rbound-1))
        (lbegin, rbegin)= find_edge(pos, denSM)

        # computing degree of segregation
        print(step, time, cal_seg(den, lbound, lbegin, rbegin, rbound, c_bar), file= OFILE)
        print("step= ", step, ",bounds: ", lbound, lbegin, rbegin, rbound, ",c-avg: ", c_bar)

print("*** Calculation completed! ***")
