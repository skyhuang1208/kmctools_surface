#!/usr/bin/env python3

#********* constants *********#
bcc= ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5));
fcc= ((0.5, 0.5, 0),    (0.5, 0, 0.5),    (0, 0.5, 0.5));
#########* constants *#########

#********* parameters ********#
nx= 660
ny= 32
nz= 32
vbra= fcc 
l_nM= 20 # vacuum layer
l_out= 2 # output visual layer

lxlo=      l_nM - 5         #left boundaies
lxhi=      l_nM + l_out -1
rxlo= nx - l_nM - l_out     #right boundaies
rxhi= nx - l_nM + 4
#########  parameters #########

from sys import argv
from os.path import isfile

if len(argv) != 2:
    print("make layers near srf.")
    print("read in history.sol file. Only output 1 and -1")
    exit("Usage: %s <history.sol>" % argv[0])

# print parameters
print("Starting...")
print("nx=", nx, ", vacuum layer=", l_nM, ", output layer=", l_out)
print("boundaries:")
print(lxlo, lxhi)

# open output file
OFILE1= open("srf1.xyz", "w")
OFILE2= open("srf2.xyz", "w")
print(rxlo, rxhi)

# read history.sol
if not isfile(argv[1]): exit("%s doesnt exist!" % argv[3])
with open(argv[1], "r") as IFILE: # Reading his file
    Ncycle= 0
    while True:
        # first line
        try: nlines= int(IFILE.readline())
        except: break # check EOF

        # second line
        (dump, step, time_)= IFILE.readline().split()

        # reading the block and construct srf
        data1= {} # data is a dictionary
        data2= {}
        nB1= 0
        nB2= 0
        for a in range(0, nlines):
            (ltcp)= IFILE.readline().split()
            ltcp= int(ltcp)
            i=  int(ltcp/nz/ny)
            if i >= lxlo and i <= lxhi: 
                data1[ltcp]= -1
                nB1 += 1
            if i >= rxlo and i <= rxhi:
                data2[ltcp]= -1
                nB2 += 1

        # output file1
        print((lxhi-lxlo+1)*ny*nz-len(data1)+nB1, file= OFILE1)
        print(dump, step, time_, file= OFILE1)
        for i in range(lxlo, lxhi+1):
            for j in range(0, ny):
                for k in range(0, nz):
                    x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
                    y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
                    z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]
                    
                    ltcp= i*ny*nz + j*nz + k
                    tp= data1.get(ltcp)
                    if tp==None: print("1",  x, y, z, file=OFILE1)
                    elif -1==tp: print("-1", x, y, z, file=OFILE1)

        # output file2
        print((rxhi-rxlo+1)*ny*nz-len(data2)+nB2, file= OFILE2)
        print(dump, step, time_, file= OFILE2)
        for i in range(rxlo, rxhi+1):
            for j in range(0, ny):
                for k in range(0, nz):
                    x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
                    y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
                    z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]
                    
                    ltcp= i*ny*nz + j*nz + k
                    tp= data2.get(ltcp)
                    if tp==None: print("1",  x, y, z, file=OFILE2)
                    elif -1==tp: print("-1", x, y, z, file=OFILE2)
        
        # print to scream the progress
        Ncycle= Ncycle+1
        if 0==Ncycle%100: print(Ncycle, "snapshots")

print("*** Calculation completed! ***")
