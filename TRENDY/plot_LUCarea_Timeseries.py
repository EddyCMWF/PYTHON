#!/bin/env python3.5

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


TRENDY_DIR = '/prj/CLIFFTOP/TRENDY/'
TRENDY_RUN='S2'
Cveg_file=TRENDY_DIR+TRENDY_RUN+'/CRU-NCEP_1p1.update.vn4.6.'+TRENDY_RUN+'.Annual.cVeg.nc'
Csoil_file=TRENDY_DIR+TRENDY_RUN+'/CRU-NCEP_1p1.update.vn4.6.'+TRENDY_RUN+'.Annual.cSoil.nc'
area_file=TRENDY_DIR+'Area.nc'

area=nc.Dataset(area_file,'r').variables['area'][:]

CVinf=nc.Dataset(Cveg_file,'r')
timevar=CVinf.variables['time']
time=nc.num2date(timevar[:],units=timevar.units,calendar=timevar.calendar)
date=np.array([ tim.year for tim in time ]) 
CV = CVinf.variables['cVeg'][:]*area*1e-12
CV = np.sum( CV.reshape( CV.shape[0],-1 ), axis=1 )
CVinf.close()

CSinf=nc.Dataset(Csoil_file,'r')
CS = CSinf.variables['cSoil'][:]*area*1e-12
CS = np.sum( CS.reshape( CS.shape[0],-1 ), axis=1 )
CSinf.close()

print(time)
fig,axes=plt.subplots( ncols=1,nrows=3,figsize=(15,20) )

ax=axes[0]
ax.plot(date,CV+CS)
ax.set_ylabel('Total Land Carbon Stock (GtC)')

ax=axes[1]
ax.plot(date,CV)
ax.set_ylabel('Veg Carbon Stock (GtC)')


ax=axes[2]
ax.plot(date,CS)
ax.set_ylabel('Soil Carbon Stock (GtC)')

fig.suptitle('TRENDY Carbon Stocks ('+TRENDY_RUN+')',fontsize=20)
fig.savefig(TRENDY_DIR+'TRENDY_CarbonStocks_'+TRENDY_RUN+'.png')

plt.show()
