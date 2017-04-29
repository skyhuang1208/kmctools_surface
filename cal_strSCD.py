#!/usr/bin/env python3

from sys import argv
from os.path import isfile
import math

if len(argv) != 6:
    print("Calculate N of clusters.")
    exit("Usage: %s <species.out> [sys_volume] [dpa_rate] [sys_name] [temperature]" % argv[0])

#********* parameters ********#
sys_vol=  float(argv[2])
dpa_rate= float(argv[3])
sys_name=	argv[4]
T=          int(argv[5])
targeted_dpas= (0.5, 1.0, 2.0, 99999999999999999999999999999)
dpa_names=  ("05", "10", "20")
atomVol= 0.5
alatt= 3.165e-10    # lattice constant (unit: m)
limit_size= 1.0e-9  # the limit size to be seen in TEM (unit: m)
burVec= math.sqrt(3)/2
### strength cal pars
alpha= (0.15, (0.25, 0.30, 0.35, 0.40), 0.60) # SIA, void, prec
M= 3.06
mu= 161000.0 # MPa
#########  parameters #########

def convert_cltr(key):
    islt0= (key>0); key= abs(key)
    nRe= key%1000; key //= 1000
    if islt0: nSIA= key; nV= 0
    else:     nSIA= 0; nV= key
    return( [nSIA, nV, nRe] )

SIAlength= lambda n: 2*alatt*(n*atomVol/math.pi/burVec)**(1/2)
VCClength= lambda n: 2*alatt*(3*n*atomVol/4/math.pi)**(1/3)
SOLlength= lambda n: 2*alatt*(3*n*atomVol/4/math.pi)**(1/3)

def cala2Nd(ctype, d, n):
    if d<limit_size: return 0 
    if ctype==1:
        if   d<4:   a= alpha[1][int(d)-1]
        else:       a= alpha[1][3]
    else:           a= alpha[ctype]
    return (a * a * d * n)

### read species.out ###
print("Start reading species.out...")
ncltrs= {}; nlines= 0; phytime= 0; iprint= 0
OFILE= open("hardening.out", "a")
with open(argv[1], "r") as IFILE: # Reading his file
    for line in IFILE:
        nlines +=1
        data= line.strip().split()
      
        if data[0]=='i_step': nlines= 1

        if nlines==1:
            dpa= phytime*dpa_rate
            if dpa > targeted_dpas[iprint]:
                ## PRINT OUT ##
                print("Output strength at %f dpa." % dpa)
                a2Nd= [0, 0, 0] # \alpha^2*d*n, !! n is number density !!
                for c, n in ncltrs.items():
                    if c[0]!=0 and c[1]!=0: exit("Error: contains both SIA and V: %d %d %d" % (c[0], c[1], c[2]))
                    elif c[0]!=0: # SIA-Re
                        if c[0]>=c[2]:  a2Nd[0] += cala2Nd( 0, SIAlength(c[0]+c[2]), n/sys_vol )
                        else:           a2Nd[2] += cala2Nd( 2, SOLlength(c[0]+c[2]), n/sys_vol )
                    elif c[1]!=0: # V-Re
                        if c[1]>=c[2]:  a2Nd[1] += cala2Nd( 1, VCClength(c[1]+c[2]), n/sys_vol )
                        else:           a2Nd[2] += cala2Nd( 2, SOLlength(c[1]+c[2]), n/sys_vol )
                    else:               a2Nd[2] += cala2Nd( 2, SOLlength(c[2]), n/sys_vol )
                hardi= 3.2*M*mu*(burVec*alatt)*math.sqrt(a2Nd[0])
                hardv= 3.2*M*mu*(burVec*alatt)*math.sqrt(a2Nd[1])
                hards= 3.2*M*mu*(burVec*alatt)*math.sqrt(a2Nd[2])
                print(sys_name, T, targeted_dpas[iprint], hardi, hardv, hards, file=OFILE)

                iprint +=1
            ncltrs={}
            
        elif nlines==2: pass
        elif nlines==3: phytime= float(data[3])
        else:
            cltr= convert_cltr(int( data[1].replace(',', '') ))
            n= int(data[3])
            key= ( cltr[0], cltr[1], cltr[2] )
            ncltrs[key]= ncltrs.get(key, 0) + n

if iprint !=3: print("Warning: not getting distributions of all targeted dpa - %d/3" % iprint)
print("**** Calculation's done. ****")
