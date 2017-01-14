#!/usr/bin/env python3

import sys

def func_fitted(x, a):        # fitting curves
    return 6*a*x
def fitting(x, y):
    import scipy.optimize as opt
    return opt.curve_fit(func_fitted, x, y)[0]

if len(sys.argv) != 6:
    print("calculate diffusivities and transport coeffs with the input of r2_VS_t")
    exit("Usage: %s <input_file> <output> #TEMP# #SOL%%# #Vtotal#" % sys.argv[0])

kB= 8.617330350e-5
temp= float(sys.argv[3])
sol= float(sys.argv[4])
vtotal= float(sys.argv[5])

x= []   # solute
y= []
x1= []  # defect
y1= []
x2= []  # transport
y2= []
x3= []  # transport_ss
y3= []
with open(sys.argv[1], "r") as IFILE: # start reading data 
    for line in IFILE:
        (step, t1, t2, r2d, r2s, rds, rs2)= line.strip().split()
        # t1: un-rescaled time; t2: rescaled time
        #r2d: r^2 of defect; r2s: r^2 of solutes; rds: rd DOT rs
        x.append(float(t2))     # solute
        y.append(float(r2s))
        x1.append(float(t1))    # defect
        y1.append(float(r2d))
        x2.append(float(t2))    # transport
        y2.append(float(rds))
        x3.append(float(t2))    # transport_ss
        y3.append(float(rs2))

Ds= fitting(x, y)[0] # solute

import statistics

OFILE1= open(sys.argv[2], "a")
chop= [10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10, 5]

for c in chop:
    p= int(len(x1)/c) # period
    D= []
    L= []
    Lss= []
    for i in range(c):
        index0= p*i
        index1= p*(i+1)-1
        if x1[index1] != x1[index0]: 
            dcoeff= (y1[index1]-y1[index0])/6.0/(x1[index1]-x1[index0])
            D.append(dcoeff)
        if x2[index1] != x2[index0]: 
            transp= (y2[index1]-y2[index0])/(6.0*vtotal*kB*temp)/(x2[index1]-x2[index0])
            L.append(transp)
        if x3[index1] != x3[index0]: 
            tra_ss= (y3[index1]-y3[index0])/(6.0*vtotal*kB*temp)/(x3[index1]-x3[index0])
            Lss.append(tra_ss)

    if len(D) != 0: avgD= statistics.mean(D); stdD= statistics.stdev(D)
    else: avgD= 0; stdD= 0
    if len(L) != 0: avgL= statistics.mean(L); stdL= statistics.stdev(L)
    else: avgL= 0; stdL= 0
    if len(Lss) != 0: avgLss= statistics.mean(Lss); stdLss= statistics.stdev(Lss)
    else: avgLss= 1; stdLss= 1

    print(temp, sol, c, avgD, stdD, Ds, 0, avgL, stdL)
    if c==100: print(temp, sol, avgD, stdD, Ds, 0, avgL, stdL, avgLss, stdLss, "T Cs Dd Dd Ds Ds L L Ls Ls", file=OFILE1)

print(temp, sol, 1, y1[-1]/6.0/x1[-1], 0, y2[-1]/(6.0*vtotal*kB*temp)/x2[-1], 0, y3[-1]/(6.0*vtotal*kB*temp)/x3[-1])
