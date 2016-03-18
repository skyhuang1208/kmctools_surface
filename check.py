#!/usr/bin/env python3

###### read kmc_par.h and change value of some parameters  #######

import sys

val_vcc= ""
val_sol= ""
val_temp= ""
with open(sys.argv[1], "r") as IFILE: # start reading and change parameters
    for line in IFILE:
        if "=" in line:
            (w0, w1)= line.strip().split("=")
            (val, *dump)= w1.split(";")
            val= val.strip()
            if "par_compV" in w0: val_vcc=  val
            if "par_compA" in w0: val_sol=  val
            if "par_temp"  in w0: val_temp= val

print("vcc= ", val_vcc, ", AAA= ", val_sol, ", temp= ", val_temp)
