#!/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats


ancil_dir='/prj/CLIFFTOP/COMMON_DATA/ANCILS/'
area_file_1d=ancil_dir+'Area_in_iris_format.nc'
area_file_2d=ancil_dir+'grid_info.nc'

IMAGE_dir='/prj/CLIFFTOP/CLUES/'
IMAGE_file=IMAGE_dir+'GREGION_IMAGE_5arcmin.asc'

# In[67]:

#Read IMAGE data:
IMinfile = open(IMAGE_file,'r')
IMlines = IMinfile.readlines()
IMinfile.close()

ncols  = int(IMlines.pop(0).split()[-1])
nrows  = int(IMlines.pop(0).split()[-1])
minlon = float(IMlines.pop(0).split()[-1])
minlat = float(IMlines.pop(0).split()[-1])
res    = float(IMlines.pop(0).split()[-1])
fillval= int(IMlines.pop(0).split()[-1])

lats_1d=(np.arange(nrows)*res)+minlat
lons_1d=(np.arange(ncols)*res)+minlon

lons_2d, lats_2d = np.meshgrid(lons_1d, lats_1d)
IMAGE_data = np.zeros_like(lats_2d)

for iline in range(nrows):
    split = IMlines[iline].split()
    data_line = np.array([ int(number) for number in split ])
    IMAGE_data[iline,:]=data_line

IMAGE_data = np.ma.masked_equal(IMAGE_data,fillval)
IMAGE_data = IMAGE_data[::-1,:]

plt.figure(figsize=[15,11])
plt.imshow(IMAGE_data, origin='bottom')
plt.colorbar()
plt.show()


#read in the 1d imogen data:
inf=nc.Dataset(area_file_1d,'a')
im_lat=inf.variables['latitude'][:].squeeze()
im_lon=inf.variables['longitude'][:].squeeze()



lat_res=2.5
lon_res=3.75
n_impoints=len(im_lat)
fill_value=-999.
im_IMAGE_1d = np.zeros_like(im_lat)

for ipt in range(n_impoints):
    lat,lon = im_lat[ipt],im_lon[ipt]
    
    cell_mask = (  (lats_2d>=lat)&(lats_2d<lat+lat_res)
                 & (lons_2d>=lon)&(lons_2d<lon+lon_res) 
                 & (IMAGE_data.mask==False)
                 )
    
    
    if len(IMAGE_data[cell_mask])>=1:
        print(ipt,lat,lon,stats.mode(IMAGE_data[cell_mask]))
        im_IMAGE_1d[ipt] = stats.mode(IMAGE_data[cell_mask])[0]
    else:
        print(ipt,lat,lon,fillval)
        im_IMAGE_1d[ipt] = fillval
    
    
badex=np.where(im_IMAGE_1d==fillval)[0]


print(badex)
print(im_lat[badex])
print(im_lon[badex])

# fix gap fills manually:
fill_1=[0,1]
im_IMAGE_1d[badex[fill_1]]=4
fill_2=[2]
im_IMAGE_1d[badex[fill_2]]=3
fill_3=[3]
im_IMAGE_1d[badex[fill_3]]=16

IMAGE_region_names=['Canada','USA','Mexico','Central America','Brazil','Rest of South America',
                    'Northern Africa','Western Africa','Eastern Africa','South Africa',
                    'Western Europe','Central Europe','Turkey','Ukraine region','Central Asia',
                    'Russia region','Middle East','India','Korea region','China region',
                    'Southeastern Asia','Indonesia region','Japan','Oceania','Rest of South Asia'
                    'Rest of Southern Africa']

#add transcomm to Area file
outvar=inf.createVariable('IMAGE_region','int',('time','y','x'))
outvar=inf.variables['IMAGE_region']
outvar.long_name='IMAGE index'
outvar.region_names=IMAGE_region_names
outvar.index_values=list(range(1,27))


outvar[:] = im_IMAGE_1d


# In[124]:

inf.close()


# In[125]:

# Read in 2d grid info to check and store in 2d
inf_2d=nc.Dataset(area_file_2d,'a')
lat_2d = inf_2d.variables['latitude'][:]
lon_2d = inf_2d.variables['longitude'][:]
land_index = inf_2d.variables['land_index'][:]


# In[130]:

im_IMAGE_2d=np.ma.masked_array(im_IMAGE_1d[land_index],mask=land_index.mask,fill_value=fillval)
im_IMAGE_2d.data[im_IMAGE_2d.mask==True]=fillval
plt.figure(figsize=[20,15])
plt.imshow(im_IMAGE_2d[:,:],origin='bottom',vmin=-5)
plt.colorbar()
plt.show()

# In[133]:

#add transcomm to Area file
outvar=inf_2d.createVariable('IMAGE_region','float32',('y','x'),fill_value=-999)
outvar.long_name='IMAGE index'
outvar[:] = im_IMAGE_2d


# In[134]:

inf_2d.close()


# In[ ]:



