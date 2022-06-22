#!/bin/env python2.7

import netCDF4 as nc
import numpy as np
import sys
import matplotlib.pyplot as plt



def plot_absolute(data1,data2,data3,      \
                  grindex,var,plotdir, \
                  year1,year2,year3):

    plt.figure(figsize=[20,5])

    plt.subplot(1,3,1)
    plt.imshow(np.ma.masked_array(data1[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year1)

    plt.subplot(1,3,2)
    plt.imshow(np.ma.masked_array(data2[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year2)

    plt.subplot(1,3,3)
    plt.imshow(np.ma.masked_array(data3[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year3)

    plt.savefig(plotdir+var+'_dump_comparison.png')
    plt.show()


def plot_diff_and_absolute(data1,data2,data3,      \
                            grindex,var,plotdir, \
                            year1,year2,year3):
    diff21=data2-data1
    diff32=data3-data2
    diff31=data3-data1

    plt.figure(figsize=[20,9])

    plt.subplot(2,3,1)
    plt.imshow(np.ma.masked_array(data1[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year1)

    plt.subplot(2,3,2)
    plt.imshow(np.ma.masked_array(data2[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year2)

    plt.subplot(2,3,3)
    plt.imshow(np.ma.masked_array(data3[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year3)

    plt.subplot(2,3,4)
    plt.imshow(np.ma.masked_array(diff21[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year2+'-'+year1)

    plt.subplot(2,3,5)
    plt.imshow(np.ma.masked_array(diff32[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year3+'-'+year2)

    plt.subplot(2,3,6)
    plt.imshow(np.ma.masked_array(diff31[grindex],mask=grindex.mask),origin='bottom')
    plt.colorbar(pad=0.04,shrink=0.8)
    plt.title(var+' '+year3+'-'+year1)

    plt.savefig(plotdir+var+'_dump_comparison.png')
    plt.show()


run=sys.argv[1]
scenario=sys.argv[2]
year1=sys.argv[3]
year2=sys.argv[4]

run='CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0'
scenario='1p5deg'
year1='1901'
year2='1902'
year3='1902'

plotdir='/group_workspaces/jasmin2/clifftop/CLIFFTOP/ECP_output/plots/'
data_dir='/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/'+run+'/'

dump_file1=data_dir+'vn4.8_imogen_'+run+'_'+scenario+'.dump.'+year1+'0101.0.nc'
print(dump_file1)
dump_file2=data_dir+'vn4.8_imogen_'+run+'_'+scenario+'.dump.'+year2+'0101.0.backup.nc'
print(dump_file2)
dump_file3=data_dir+'vn4.8_imogen_'+run+'_'+scenario+'.dump.'+year3+'0101.0.nc'
print(dump_file3)

year2+='-orig'
year3+='-new'

grid_file='/group_workspaces/jasmin2/clifftop/COMMON_DATA/ANCILS/grid_info.nc'
grinf=nc.Dataset(grid_file,'r')
grindex=grinf.variables['land_index'][:]
Area=grinf.variables['Area'][:]

var_list=[]
DATA_DICT_1={}
inf1=nc.Dataset(dump_file1,'r')
for var in inf1.variables:
    var_list.append(str(var))
    DATA_DICT_1[str(var)]=inf1.variables[var][:]
inf1.close()

DATA_DICT_2={}
inf2=nc.Dataset(dump_file2,'r')
for var in inf2.variables:
    DATA_DICT_2[str(var)]=inf2.variables[var][:]
inf2.close()

DATA_DICT_3={}
inf3=nc.Dataset(dump_file3,'r')
for var in inf3.variables:
    DATA_DICT_3[str(var)]=inf3.variables[var][:]
inf3.close()


for var in var_list:
    diff=DATA_DICT_2[var]-DATA_DICT_1[var]
    print(var, np.min(diff),np.mean(diff),np.max(diff), np.std(diff) )
    
for var in var_list:
    diff=DATA_DICT_3[var]-DATA_DICT_1[var]
    print(var, np.min(diff),np.mean(diff),np.max(diff), np.std(diff) )

for var in var_list:
    diff=DATA_DICT_3[var]-DATA_DICT_2[var]
    print(var, np.min(diff),np.mean(diff),np.max(diff), np.std(diff) )

map_var_list=['canopy','cs','gs','snow_tile','t_soil','tstar_tile', 'sthuf','lai','canht','frac_agr_prev','wood_prod_fast','frac_past_prev','sthzw','zw','rgrain','tsoil_deep','seed_rain','cv','tsnow','rho_snow','snow_depth', 'snow_grnd','nsnow','frac',]
for var in map_var_list:
    print(var,DATA_DICT_1[var].shape)



soil_map_vars=['t_soil', 'sthuf','tsoil_deep']
tile_map_vars=['canopy','snow_tile','tstar_tile','rgrain','rho_snow','snow_depth', 'snow_grnd','nsnow','frac']
pft_map_vars=['lai','canht']
gb_map_vars=['sthzw','zw','frac_agr_prev','wood_prod_fast','frac_past_prev','cv']

var='frac'
max_frac_1=np.argmax(DATA_DICT_1[var],axis=0)
max_frac_2=np.argmax(DATA_DICT_2[var],axis=0)
max_frac_3=np.argmax(DATA_DICT_3[var],axis=0)

max_pft_1=np.argmax(DATA_DICT_1[var][:13,:],axis=0)
max_pft_2=np.argmax(DATA_DICT_2[var][:13,:],axis=0)
max_pft_3=np.argmax(DATA_DICT_3[var][:13,:],axis=0)


data1=np.ma.masked_array(max_frac_1[grindex],mask=grindex.mask)
data2=np.ma.masked_array(max_frac_2[grindex],mask=grindex.mask)
data3=np.ma.masked_array(max_frac_3[grindex],mask=grindex.mask)


plt.figure(figsize=[20,4])

plt.subplot(1,3,1)
plt.imshow(data1,origin='bottom')
plt.colorbar(pad=0.04,shrink=0.8)
plt.title(var+' '+year1)

plt.subplot(1,3,2)
plt.imshow(data2,origin='bottom')
plt.colorbar(pad=0.04,shrink=0.8)
plt.title(var+' '+year2)

plt.subplot(1,3,3)
plt.imshow(data3,origin='bottom')
plt.colorbar(pad=0.04,shrink=0.8)
plt.title(var+' '+year3)

plt.savefig(plotdir+var+'_dump_comparison.png')
plt.show()

for var in gb_map_vars:
    data1=DATA_DICT_1[var]
    data2=DATA_DICT_2[var]
    data3=DATA_DICT_3[var]
    
    plot_diff_and_absolute(data1,data2,data3,grindex,var,plotdir,year1,year2,year3)


for var in soil_map_vars:
    data1=np.mean(DATA_DICT_1[var][:3,:],axis=0)
    data2=np.mean(DATA_DICT_2[var][:3,:],axis=0)
    data3=np.mean(DATA_DICT_3[var][:3,:],axis=0)
    
    plot_diff_and_absolute(data1,data2,data3,grindex,var,plotdir,year1,year2,year3)


npts=len(max_frac_1)

for var in tile_map_vars:
    data1=np.array( [ DATA_DICT_1[var][max_frac_1[i],i] for i in range(npts) ] )
    data2=np.array( [ DATA_DICT_2[var][max_frac_2[i],i] for i in range(npts) ] )
    data3=np.array( [ DATA_DICT_3[var][max_frac_3[i],i] for i in range(npts) ] )

    plot_diff_and_absolute(data1,data2,data3,grindex,var,plotdir,year1,year2,year3)


for var in pft_map_vars:
    data1=np.array( [ DATA_DICT_1[var][max_pft_1[i],i] for i in range(npts) ] )
    data2=np.array( [ DATA_DICT_2[var][max_pft_2[i],i] for i in range(npts) ] )
    data3=np.array( [ DATA_DICT_3[var][max_pft_3[i],i] for i in range(npts) ] )

    plot_diff_and_absolute(data1,data2,data3,grindex,var,plotdir,year1,year2,year3)










