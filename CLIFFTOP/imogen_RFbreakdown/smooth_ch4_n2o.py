#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt


SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'
infile1=SCENARIO_DIR+'rcp2.6_concs_ch4_n2o_vn1p1.txt'
infile2=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o.txt'
outfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o_vn1p1.txt'

trans_year=2005

# Read in RCP data (for 1850-2005)
lines=open(infile1,'r').readlines()
year1,rcp_ch4,rcp_n2o=[],[],[]
for line in lines:
    split=line.replace('\n','').split()
    year1.append(int(split[0]))
    rcp_ch4.append(float(split[1]))
    rcp_n2o.append(float(split[2]))

# Read in SSP data (for 2005-2100)
lines=open(infile2,'r').readlines()
year2,ssp_ch4,ssp_n2o=[],[],[]
for line in lines:
    split=line.replace('\n','').split()
    year2.append(int(split[0]))
    ssp_ch4.append(float(split[1]))
    ssp_n2o.append(float(split[2]))

trans_index=year1.index(trans_year)

ch4_smooth=rcp_ch4.copy()
ch4_smooth[trans_index:]=ssp_ch4[trans_index:]

n2o_smooth=rcp_n2o.copy()
n2o_smooth[trans_index:]=ssp_n2o[trans_index:]

mean_range=1
n_smooths=10
start_point=trans_index-2
for j in range(n_smooths):
    for i in range(start_point,len(ch4_smooth)):
        if i<mean_range:
            ch4_smooth[i]=np.mean(ch4_smooth[:i+mean_range+1])
            n2o_smooth[i]=np.mean(n2o_smooth[:i+mean_range+1])
        elif (i>len(ch4_smooth)-mean_range):
            ch4_smooth[i]=np.mean(ch4_smooth[i-mean_range:])
            n2o_smooth[i]=np.mean(n2o_smooth[i-mean_range:])
        else:
            ch4_smooth[i]=np.mean(ch4_smooth[i-mean_range:i+mean_range+1])
            n2o_smooth[i]=np.mean(n2o_smooth[i-mean_range:i+mean_range+1])


plt.plot(year2,ssp_ch4,label='SSP2')
plt.plot(year1,rcp_ch4,label='RCP2.6')
plt.plot(year1,ch4_smooth,label='SSP2/RCP2.6-smooth',lw=1.5)
plt.xlabel('CH4 (ppbv)')
plt.legend(loc=2)
plt.show()

plt.plot(year2,ssp_n2o,label='SSP2')
plt.plot(year1,rcp_n2o,label='RCP2.6')
plt.plot(year1,n2o_smooth,label='SSP2/RCP2.6-smooth',lw=1.5)
plt.xlabel('N2O (ppbv)')
plt.legend(loc=2)
plt.show()

outf=open(outfile,'w')
for yr,CH4,N2O in zip(year1,ch4_smooth,n2o_smooth):
    outf.write('%4i  %10.5f %10.5f\n'%(yr,CH4,N2O))

outf.close()

