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
#########  parameters #########

#********* calculate degree of seg *********#
#####***** calculate degree of seg *****#####

from sys import argv
from os.path import isfile

if len(argv) != 2:
    print("make plotable files for srf.\n split srf to 2 and convert to xyz coordinate")
    exit("Usage: %s <history.srf>" % argv[0])

# open output file
OFILE1= open("srf1.xyz", "w")
OFILE2= open("srf2.xyz", "w")

# read and calculate
if not isfile(argv[1]): exit("%s doesn\'t exist!" % argv[1])
with open(argv[1], "r") as IFILE: # Reading his file
    n= 0
    while True:
        data1= []
        data2= []
        
        # first line
        try: nlines= int(IFILE.readline())
        except: break # check EOF

        # second line
        line2= IFILE.readline().rstrip()
    
        # reading the block
        for a in range(0, nlines):
            (tp, i, j, k)= IFILE.readline().split()
            i= int(i)
            j= int(j)
            k= int(k)
            x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
            y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
            z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]
            
            if i<nx/2.0: data1.append([tp, x, y, z])
            else:        data2.append([tp, x, y, z])
         
        # output to ofile1
        print(len(data1), file= OFILE1)
        print(line2, file= OFILE1)
        for data in data1:
            print(data[0], data[1], data[2], data[3], file= OFILE1)

        # output to ofile2
        print(len(data2), file= OFILE2)
        print(line2, file= OFILE2)
        for data in data2:
            print(data[0], data[1], data[2], data[3], file= OFILE2)

        # print times calculated
        n= n+1
        if n%100==0: print("Snapshots: ", n)

print("*** Calculation completed! ***")
