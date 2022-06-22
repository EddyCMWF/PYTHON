#!/usr/bin/env python

import numpy as np


def func_MN(M,N):
    #func_MN = 0.47 * np.log( 1 +  (0.6356* ( (M*N*1e-6)**0.75 ) )       \
    #                        + (0.007 * (M*1e-3) * ( (M*N*1e-6)**1.52) ) )
    func_MN = 0.47 * np.log( 1 +  (2.01e-5*((M*N)**0.75))       \
                            + (5.31e-15*M*((M*N)**1.52) ) )
    return func_MN



#Define constants:
alpha_ch4=0.036


#Read in CH4 and N2O data
infile='/users/eow/edwcom/CLIFFTOP/Tprofiles/ch4_n2o_conc.csv'
inf=open(infile,'r')
lines=inf.readlines()
inf.close()
years,ch4,n2o,rcpRF=[],[],[],[]
for line in lines:
    split=line.split(',')
    years.append(int(split[0]))
    ch4.append(float(split[1]))
    n2o.append(float(split[2]))
    rcpRF.append(float(split[3]))

years=np.array(years)

ch4=np.array(ch4)
ch4_0=np.zeros_like(ch4)+ch4[0]

n2o=np.array(n2o)
n2o_0=np.zeros_like(n2o)+n2o[0]

rcpRF=np.array(rcpRF)

#f_ch4_n2o0  = func_MN(ch4,n2o[0])
#f_ch40_n2o0 = func_MN(ch4[0],n2o[0])
f_ch4_n2o0  = func_MN(ch4,n2o_0)
f_ch40_n2o0 = func_MN(ch4_0,n2o_0)


delRF_ch4_n2o =  ( alpha_ch4*(np.sqrt(ch4)-np.sqrt(ch4_0)) )  \
               - ( f_ch4_n2o0 - f_ch40_n2o0 )

delRF_ch4_h2o = alpha_ch4*(np.sqrt(ch4)-np.sqrt(ch4_0))*0#.15

delRF_ch4=delRF_ch4_n2o+delRF_ch4_h2o

for i in range(270):
    print(years[i],ch4[i],n2o[i],rcpRF[i],delRF_ch4[i],rcpRF[i]-delRF_ch4[i])
 


import matplotlib.pyplot as plt
year_index=np.where(years==1765)[0]
fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(8,4))

ax.plot(years[:300],rcpRF[:300],label='rcp3PD')
#ax.plot(years[:300],delRF_ch4[:300],label='Calculated')
ax.plot(years[:300],delRF_ch4[:300]-delRF_ch4[year_index],label='Calculated')
#ax.plot(years[:300],rcpRF[:300]-delRF_ch4[:300],label='difference')

ax.legend(loc=2)

plt.show()





