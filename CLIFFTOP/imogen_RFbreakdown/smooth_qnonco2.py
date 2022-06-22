#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt


infile='SSP2-2.6_IMAGE_qnonco2.txt'
outfile='SSP2-2.6_IMAGE_qnonco2_smooth.txt'

lines=open(infile,'r').readlines()

year,qnonco2=[],[]
for line in lines:
    split=line.split()
    year.append(int(split[0]))
    qnonco2.append(float(split[1].replace('\n','')))

qnonco2_smooth=qnonco2.copy()
mean_range=1
n_smooths=10
for j in range(n_smooths):
    for i in range(len(qnonco2)):
        if i<mean_range:
            qnonco2_smooth[i]=np.mean(qnonco2_smooth[:i+mean_range+1])
        elif (i>len(qnonco2)-mean_range):
            qnonco2_smooth[i]=np.mean(qnonco2_smooth[i-mean_range:])
        else:
            qnonco2_smooth[i]=np.mean(qnonco2_smooth[i-mean_range:i+mean_range+1])
            #print(i)
            #print(i-mean_range,i+mean_range+1)


plt.plot(year,qnonco2,label='SSP2')
plt.plot(year,qnonco2_smooth,label='SSP2-smooth',lw=1.5)
plt.xlabel('$Q_{non-co2}$ (W m$^{-2}$)')
plt.legend(loc=2)
plt.show()

outf=open(outfile,'w')
for yr,Q in zip(year,qnonco2_smooth):
    outf.write('%4i  %10.5f\n'%(yr,Q))

outf.close()

