#!/bin/env python2.7

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

plotdir='/group_workspaces/jasmin2/clifftop/CLIFFTOP/ECP_output/plots/'
data_dir='/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/'
gcm='CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0'
scenario='1p5deg'
year1='1902'
year2='1903'
year3='1904'
file1=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Annual_carbon.'+year1+'.nc'
file2=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Annual_carbon.'+year2+'.nc'
file3=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Annual_carbon.'+year3+'.nc'

file1=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_carbon.'+year1+'.nc'
file2=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_carbon.'+year2+'.nc'
file3=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_carbon.'+year3+'.nc'

file1=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_h2o_t.'+year1+'.nc'
file2=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_h2o_t.'+year2+'.nc'
file3=data_dir+gcm+'/vn4.8_imogen_'+gcm+'_'+scenario+'.Monthly_h2o_t.'+year3+'.nc'

grid_file='/group_workspaces/jasmin2/clifftop/COMMON_DATA/ANCILS/grid_info.nc'


inf1=nc.Dataset(file1,'r')
inf2=nc.Dataset(file2,'r')
inf3=nc.Dataset(file3,'r')
grinf=nc.Dataset(grid_file,'r')

grindex=grinf.variables['land_index'][:]
Area=grinf.variables['Area'][:]

var='cv'
option=0
level=0

if option==0:
    gpp1=inf1.variables[var][0,:].squeeze()
    gpp2=inf2.variables[var][0,:].squeeze()
    gpp3=inf3.variables[var][0,:].squeeze()
elif option==1:
    gpp1=inf1.variables[var][0,level,:].squeeze()
    gpp2=inf2.variables[var][0,level,:].squeeze()
    gpp3=inf3.variables[var][0,level,:].squeeze()
elif option==2:
    gpp1=np.max(inf1.variables[var][0,:].squeeze(),axis=0)
    gpp2=np.max(inf2.variables[var][0,:].squeeze(),axis=0)
    gpp3=np.max(inf3.variables[var][0,:].squeeze(),axis=0)

    
gpp1_2D=np.ma.masked_array(gpp1[grindex],mask=grindex.mask)
gpp2_2D=np.ma.masked_array(gpp2[grindex],mask=grindex.mask)
gpp3_2D=np.ma.masked_array(gpp3[grindex],mask=grindex.mask)


diff21=gpp2-gpp1
diff21_2D=np.ma.masked_array(diff21[grindex],mask=grindex.mask)

diff31=gpp3-gpp1
diff31_2D=np.ma.masked_array(diff31[grindex],mask=grindex.mask)

diff32=gpp3-gpp2
diff32_2D=np.ma.masked_array(diff32[grindex],mask=grindex.mask)


plt.figure(figsize=[20,10])
plt.subplot(231)
plt.imshow(diff21_2D,origin='bottom')
plt.title('Diff 2-1')
plt.colorbar()

plt.subplot(232)
plt.imshow(diff32_2D,origin='bottom')
plt.title('Diff 3-2')
plt.colorbar()

plt.subplot(233)
plt.imshow(diff31_2D,origin='bottom')
plt.title('Diff 3-1')
plt.colorbar()

plt.subplot(234)
plt.imshow(gpp1_2D,origin='bottom')
plt.title(var+' '+year1)
plt.colorbar()

plt.subplot(235)
plt.imshow(gpp2_2D,origin='bottom')
plt.title(var+' '+year2)
plt.colorbar(pad=0.04)

plt.subplot(236)
plt.imshow(gpp3_2D,origin='bottom')
plt.title(var+' '+year3)
plt.colorbar()

plt.savefig(plotdir+var+'_comparison_maps.png')
plt.show()




inf1.close()
inf2.close()
inf3.close()




