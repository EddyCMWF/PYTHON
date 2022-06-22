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


spin1='07'
spin2='08'
start_year=1850
end_year=1860

plotdir='/group_workspaces/jasmin2/clifftop/CLIFFTOP/ECP_output/plots/'
data_dir='/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/'

dump_file1=data_dir+'vn4.8_imogen_spinup'+spin1+'.dump.YYYY0101.0.nc'
print(dump_file1)
dump_file2=data_dir+'vn4.8_imogen_spinup'+spin2+'.dump.YYYY0101.0.nc'
print(dump_file2)


grid_file='/group_workspaces/jasmin2/clifftop/COMMON_DATA/ANCILS/grid_info.nc'
grinf=nc.Dataset(grid_file,'r')
grindex=grinf.variables['land_index'][:]
Area=grinf.variables['Area'][:]

var_list=['sthuf','cv','cs','t_soil']
DATA_DICT_1={var:[] for var in var_list}
DATA_DICT_2={var:[] for var in var_list}

for year in range(1850,1861):
    inf1=nc.Dataset(dump_file1.replace('YYYY',str(year)),'r')
    inf2=nc.Dataset(dump_file2.replace('YYYY',str(year)),'r')
    for var in var_list:
        DATA_DICT_1[var].append(inf1.variables[var][:])
        DATA_DICT_2[var].append(inf2.variables[var][:])

    inf1.close()
    inf2.close()

DATA_DICT={}
for var in var_list:
    DATA_DICT_1[var]=np.array(DATA_DICT_1[var])
    DATA_DICT_2[var]=np.array(DATA_DICT_1[var])
    
    DATA_DICT[var]=np.append(DATA_DICT_1[var],DATA_DICT_2[var],axis=0)


#plot time series of sthuf for each layer
nlayers=14
var='sthuf'
for i in range(14):
    TS=np.mean(DATA_DICT[var][:,i,:],axis=1)
    plt.plot(TS,label=i)
plt.legend(loc=8,ncol=14)
plt.show()

#plot time series of sthuf for each layer
nlayers=14
var='t_soil'
for i in range(14):
    TS=np.mean(DATA_DICT[var][:,i,:],axis=1)
    plt.plot(TS,label=i)
plt.legend(loc=8,ncol=14)
plt.show()

#plot time series of sthuf for each layer
nlayers=14
var='cv'
TS=np.mean(DATA_DICT[var][:],axis=1)
plt.plot(TS,label=var)
plt.legend(loc=8,ncol=14)
plt.show()



