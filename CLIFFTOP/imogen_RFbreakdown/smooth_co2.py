#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt


SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'
infile1=SCENARIO_DIR+'rcp2.6_concs_co2_vn1p1.txt'
infile2=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2.txt'
outfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2_vn1p1.txt'

trans_year=2005

# Read in RCP data (for 1850-2005)
lines=open(infile1,'r').readlines()
year1,rcp_co2=[],[]
for line in lines:
    split=line.replace('\n','').split()
    year1.append(int(split[0]))
    rcp_co2.append(float(split[1]))

# Read in SSP data (for 2005-2100)
lines=open(infile2,'r').readlines()
year2,ssp_co2=[],[]
for line in lines:
    split=line.replace('\n','').split()
    year2.append(int(split[0]))
    ssp_co2.append(float(split[1]))

trans_index=year1.index(trans_year)

co2_smooth=rcp_co2.copy()
co2_smooth[trans_index:]=ssp_co2[trans_index:]

mean_range=1
n_smooths=10
start_point=trans_index-2
for j in range(n_smooths):
    for i in range(len(co2_smooth)):
        if i<mean_range:
            co2_smooth[i]=np.mean(co2_smooth[:i+mean_range+1])
        elif (i>len(co2_smooth)-mean_range):
            co2_smooth[i]=np.mean(co2_smooth[i-mean_range:])
        else:
            co2_smooth[i]=np.mean(co2_smooth[i-mean_range:i+mean_range+1])
            #print(i)
            #print(i-mean_range,i+mean_range+1)


plt.plot(year2,ssp_co2,label='SSP2')
plt.plot(year1,rcp_co2,label='RCP2.6')
plt.plot(year1,co2_smooth,label='SSP2/RCP2.6-smooth',lw=1.5)
plt.xlabel('CO2 (ppmv)')
plt.legend(loc=2)
plt.show()

outf=open(outfile,'w')
for yr,CO2 in zip(year1,co2_smooth):
    outf.write('%4i  %10.5f\n'%(yr,CO2))

outf.close()

