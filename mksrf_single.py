#!/usr/bin/env python3

#********* data analysis functions *********#
#####***** data analysis functions *****#####


#********* constants *********#
bcc= ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5));
fcc= ((0.5, 0.5, 0),    (0.5, 0, 0.5),    (0, 0.5, 0.5));
#########* constants *#########

#********* parameters ********#
nx= 660
vbra= fcc 
l_nM= 20 # vacuum layer
l_out= 2 # output visual layer

lxlo=      l_nM - 5         #left boundaies
lxhi=      l_nM + l_out -1
rxlo= nx - l_nM - l_out     #right boundaies
rxhi= nx - l_nM + 4
#########  parameters #########

#********* calculate degree of seg *********#
#####***** calculate degree of seg *****#####

from sys import argv
from os.path import isfile

if len(argv) != 2:
    print("make plotable files for srf.\ninput .ltcp file!!")
    exit("Usage: %s <.ltcp>" % argv[0])

# print parameters
print("Starting...")
print("nx=", nx, ", vacuum layer=", l_nM, ", output layer=", l_out)
print("boundaries:")
print(lxlo, lxhi)
print(rxlo, rxhi)

# read and calculate
data1= []
data2= []
if not isfile(argv[1]): exit("%s doesn\'t exist!" % argv[1])
with open(argv[1], "r") as IFILE: # Reading his file
    # first line
    nlines= int(IFILE.readline())

    # second line
    (dump, step, time_)= IFILE.readline().split()
    
    # reading the block
    for a in range(0, nlines):
        (tp, i, j, k, *temp)= IFILE.readline().split()
        tp= int(tp)
        i= int(i)
        j= int(j)
        k= int(k)
        x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
        y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
        z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]
            
        if (i>=lxlo and i<=lxhi) and (1==tp or -1==tp): data1.append([tp, x, y, z])
        if (i>=rxlo and i<=rxhi) and (1==tp or -1==tp): data2.append([tp, x, y, z])

# open output file
OFILE1= open(step+"_srf1.xyz", "w")
OFILE2= open(step+"_srf2.xyz", "w")
        
# output to ofile1
print(len(data1), file= OFILE1)
print(dump, step, time_, file= OFILE1)
for data in data1:
    print(data[0], data[1], data[2], data[3], file= OFILE1)

# output to ofile2
print(len(data2), file= OFILE2)
print(dump, step, time_, file= OFILE2)
for data in data2:
    print(data[0], data[1], data[2], data[3], file= OFILE2)

print("*** Calculation completed! ***")
