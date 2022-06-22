#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt


SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'
infile=SCENARIO_DIR+'SSP2-2.6_IMAGE_dq_aerosol.txt'
outfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_dq_aerosol_smooth.txt'

trans_year=2005

lines=open(infile,'r').readlines()
year,dq_aero=[],[]
for line in lines:
    split=line.replace('\n','').split()
    year.append(int(split[0]))
    dq_aero.append(float(split[1]))

dq_aero_smooth=np.copy(dq_aero)

mean_range=1
n_smooths=20
for j in range(n_smooths):
    for i in range(len(dq_aero_smooth)):
        if i<mean_range:
            dq_aero_smooth[i]=np.mean(dq_aero_smooth[:i+mean_range+1])
        elif (i>len(dq_aero_smooth)-mean_range):
            dq_aero_smooth[i]=np.mean(dq_aero_smooth[i-mean_range:])
        else:
            dq_aero_smooth[i]=np.mean(dq_aero_smooth[i-mean_range:i+mean_range+1])
            #print(i)
            #print(i-mean_range,i+mean_range+1)


plt.plot(year,dq_aero,label='original')
plt.plot(year,dq_aero_smooth,label='smooth',lw=1.5)
plt.xlabel('year')
plt.ylabel('dQ - Aerosol')
plt.legend(loc=2)
plt.show()

outf=open(outfile,'w')
for yr,AERO in zip(year,dq_aero_smooth):
    outf.write('%4i  %10.5f\n'%(yr,AERO))

outf.close()

