#!/usr/bin/env python3

from sys import argv
from os.path import isfile
import math

if len(argv) != 6 and len(argv) != 7:
    print("Calculate N of clusters.")
    exit("Usage: %s <species.out> [sys_volume] [dpa_rate] [sys_name] [temperature] isPRINTall(1:on)(optional)" % argv[0])

#********* parameters ********#
sys_vol=  float(argv[2])
dpa_rate= float(argv[3])
sys_name=	argv[4]
T=          int(argv[5])
if len(argv)==7 and argv[6]=="1":   isALL= True
else:                               isALL= False
targeted_dpas= (0.5, 1.0, 2.0, 99999999999999999999999999999)
dpa_names=  ("05", "10", "20")
atomVol= 0.5
alatt= 0.3165 # lattice constant (unit: nm)
lbin= 0.1
limit_size= 1.0 # the limit size to be seen in TEM (unit: nm)
burVec= math.sqrt(3)/2
#########  parameters #########

def convert_cltr(key):
    islt0= (key>0); key= abs(key)
    nRe= key%1000; key //= 1000
    if islt0: nSIA= key; nV= 0
    else:     nSIA= 0; nV= key
    return( [nSIA, nV, nRe] )

SIAlength= lambda n: (n*atomVol/math.pi/burVec)**(1/2)
VCClength= lambda n: (3*n*atomVol/4/math.pi)**(1/3)
SOLlength= lambda n: (3*n*atomVol/4/math.pi)**(1/3)

### read species.out ###
if not isfile(argv[1]): exit("Input file %s doesnt exist!" % argv[1])
OFILE1= open(sys_name+str(T)+"k.nden", "a")

print("Start reading species.out...")
if isALL: print("isALL turned on. Only output density with all.")
ncltrs= {}; nlines= 0; phytime= 0; iprint= 0
with open(argv[1], "r") as IFILE: # Reading his file
    for line in IFILE:
        nlines +=1
        data= line.strip().split()
      
        if data[0]=='i_step': nlines= 1

        if nlines==1:
            ldist= [ {}, {}, {} ] # length distribution: SIA, VCC, SOL
            cltrDEN= [0, 0, 0]      # cluster density of SIA, VCC, SOL
            for c, n in ncltrs.items():
                if c[0]!=0 and c[1]!=0: exit("Error: contains both SIA and V: %d %d %d" % (c[0], c[1], c[2]))
                elif c[0]!=0: # SIA-Re
                    if c[0]>=c[2]:  ctype= 0; r= SIAlength(c[0]+c[2])
                    else:           ctype= 2; r= SOLlength(c[0]+c[2])
                elif c[1]!=0: # V-Re
                    if c[1]>=c[2]:  ctype= 1; r= VCClength(c[1]+c[2])
                    else:           ctype= 2; r= SOLlength(c[1]+c[2])
                else:               ctype= 2; r= SOLlength(c[2])
                
                l= 2*r*alatt # get diameter in nm
                key= l//lbin # sum in bin index
                ldist[ctype][key]= ldist[ctype].get(key,0) + n
                if isALL:
                    if (c[0]+c[1]+c[2])!=1: cltrDEN[ctype] += n
                else:
                    if l>limit_size: cltrDEN[ctype] += n
            ncltrs= {}
            
            ## PRINT OUT ##
            dpa= phytime*dpa_rate
            print(T, dpa, cltrDEN[0]/sys_vol, cltrDEN[1]/sys_vol, cltrDEN[2]/sys_vol, file=OFILE1)
            if (not isALL) and dpa > targeted_dpas[iprint]:
                print("Output distribution at %f dpa." % dpa)
                OFILE2= [0, 0, 0]
                OFILE2[0]= open(sys_name+str(T)+"k"+dpa_names[iprint]+"dpa.disi", "w")
                OFILE2[1]= open(sys_name+str(T)+"k"+dpa_names[iprint]+"dpa.disv", "w")
                OFILE2[2]= open(sys_name+str(T)+"k"+dpa_names[iprint]+"dpa.diss", "w")

                ibinSINGLE= ( (2*SIAlength(1)*alatt)//lbin, (2*VCClength(1)*alatt)//lbin, (2*SOLlength(1)*alatt)//lbin )
                for ctype, dist in enumerate(ldist):
                    if dist.get(ibinSINGLE[ctype])==None: 
                        ibin= 0 # index of previous bin in dist
                    else:
                        ibin= ibinSINGLE[ctype]
                        print(0, dist[ibinSINGLE[ctype]]/sys_vol, file=OFILE2[ctype])

                    dsorted= sorted(dist.items(), key= lambda x: x[0]) # sort by length
                    for i, n in dsorted: # index of bin, numbers
                        l= i*lbin
                        if (i-ibin) > 1: # if prev and current not adjacent, cut them with 1s
                            print((ibin+1)*lbin, 1, file=OFILE2[ctype])
                            print(l, 1, file= OFILE2[ctype])
                        print(l, n/sys_vol, file=OFILE2[ctype])
                        print(l+lbin, n/sys_vol, file=OFILE2[ctype])
                        ibin= i
                    print(l+lbin, 1, file=OFILE2[ctype])
                    
                iprint +=1
            
        elif nlines==2: pass
        elif nlines==3: phytime= float(data[3])
        else:
            cltr= convert_cltr(int( data[1].replace(',', '') ))
            n= int(data[3])
            key= ( cltr[0], cltr[1], cltr[2] )
            ncltrs[key]= ncltrs.get(key, 0) + n

if iprint !=3: print("Warning: not getting distributions of all targeted dpa - %d/3" % iprint)
print("**** Calculation's done. ****")
