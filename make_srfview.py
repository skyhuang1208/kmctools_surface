#!/usr/bin/env python3

import sys
import os.path

if len(sys.argv) !=2 and len(sys.argv) !=3 and len(sys.argv) !=4:
    print("Make history.srf to become a vitalization file")
    exit("Usage: %s <history.srf> [N_blocks_out] [is_projected]" % sys.argv[0])

nameIn= sys.argv[1]
if len(sys.argv)==4:
    N_blocks= int(sys.argv[2])
    is_pro= int(sys.argv[3])
elif len(sys.argv)==3:
    N_blocks= int(sys.argv[2])
    is_pro= 0
else:
    N_blocks= 1
    is_pro= 0

def output(line):
    if(is_pro):
        (tp, x, y, z, *dump)= line.split()
        print(tp, x, y, "0")
    else:
        sys.stdout.write("%s" % line)

N_max= 0
idb= 0 # id in block
with open(nameIn, "r") as IN: # Reading history files
    for line in IN:
        idb += 1

        if idb == 1:
            N_atoms= int(line)
            if N_atoms > N_max: N_max= N_atoms
        elif idb == N_atoms+2:
            idb= 0

idb= 0
i_blocks= 0
with open(nameIn, "r") as IN: # Reading history files
    for line in IN:
        idb += 1

        if idb == 1:
            N_atoms= int(line)
            if i_blocks%N_blocks == 0: 
                print(N_max)
        elif i_blocks%N_blocks == 0:
            if(idb==2): print(line)
            else:       output(line)

        if idb == N_atoms+2:
            if i_blocks%N_blocks == 0:
                for i in range(N_atoms, N_max):
                    output(line)
            idb= 0
            i_blocks += 1

