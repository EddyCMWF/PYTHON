#!/bin/env python

import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats
import pdb

ancil_dir='/prj/CLIFFTOP/COMMON_DATA/ANCILS/'
area_file_1d=ancil_dir+'Area_in_iris_format.nc'
area_file_2d=ancil_dir+'grid_info.nc'

SWAMPS_dir='/prj/CLIFFTOP/SWAMPS/'
SWAMPS_infile=SWAMPS_dir+'fw_swamps-glwd_2000-2012.nc'
SWAMPS_outfile=SWAMPS_dir+'fw_swamps-glwd_2000-2012_IMOGENgrid.nc'
SWAMPS_outfile_2d=SWAMPS_dir+'fw_swamps-glwd_2000-2012_IMOGENgrid_2d.nc'

# In[67]:

#Read IMAGE data:
SWAMPSinf = nc.Dataset(SWAMPS_infile,'r')
SWAMPSlons = SWAMPSinf.variables['lon'][:]
SWAMPSlats = SWAMPSinf.variables['lat'][:]
SWAMPSdata = SWAMPSinf.variables['Fw'][:]
SWAMPStime = SWAMPSinf.variables['time'][:]
SWAMPStimeunits = SWAMPSinf.variables['time'].units
SWAMPStimecalendar= SWAMPSinf.variables['time'].calendar
SWAMPSinf.close()

nlats=len(SWAMPSlats)
nlons=len(SWAMPSlons)
ntsteps=len(SWAMPStime)
lons_2d, lats_2d = np.meshgrid(SWAMPSlons, SWAMPSlats)

SWAMPS_data = np.ma.masked_equal(SWAMPSdata,0)
print(SWAMPS_data.shape)

SWAMPS_annmean = np.mean(SWAMPS_data.reshape(-1,12,nlats,nlons),axis=1)
SWAMPS_annmax = np.max(SWAMPS_data.reshape(-1,12,nlats,nlons),axis=1)

#plt.figure(figsize=[15,11])
#plt.imshow(SWAMPS_data[0,:]) # , origin='bottom')
#plt.colorbar()
#plt.show()

#read in the 1d imogen data:
inf=nc.Dataset(area_file_1d,'r')
im_lat=inf.variables['latitude'][:].squeeze()
im_lon=inf.variables['longitude'][:].squeeze()
im_area=inf.variables['area'][:].squeeze()
inf.close()

lat_radias=lats_2d*np.pi/180.
# calculate the area of the gid
rearth=6371000.
latdel=0.5
londel=0.5
SWAMPSarea=rearth*rearth*np.radians(londel)*(np.sin(np.radians(lats_2d+latdel))-np.sin(np.radians(lats_2d)))

#pdb.set_trace()
SWAMPS_wetlarea = SWAMPSdata*SWAMPSarea

lat_res=2.5
lon_res=3.75
n_impoints=len(im_lat)
fill_value=-999.
im_SWAMPS_1d = np.zeros([ntsteps,n_impoints])

for ipt in range(n_impoints):
    lat,lon = im_lat[ipt],im_lon[ipt]
    
    #Sum the Fwetl of fully enclosed pixels:
    ur_y = np.where( SWAMPSlats-0.25 == lat )[0][0]
    ll_y = np.where( SWAMPSlats+0.25 == lat+lat_res )[0][0]
    n_y = ur_y-ll_y+1
    ll_x = np.min( np.where( SWAMPSlons-0.25 >=lon )[0] )
    ur_x = np.max( np.where( SWAMPSlons+0.25 <=lon+lon_res )[0] )
    n_x = ur_x-ll_x+1
     
    #print(ll_y,ll_x,ur_y,ur_x)
    GC_wetl_area = np.sum(SWAMPS_wetlarea[:,ll_y:ur_y+1,ll_x:ur_x+1].reshape(-1,n_y*n_x),axis=1)
    GC_area = np.sum(SWAMPSarea[ll_y:ur_y+1,ll_x:ur_x+1])

    #print(ipt,lat,lon,GC_wetl_area,GC_area)
    # Add the half values for the border column,
    #   check to see if column is to the East or west
    if np.round(lon*2)==lon*2:
        # Column to East
        #print('East')
        c_x = np.where(SWAMPSlons-0.25 == lon)[0][0]
    else:
        # Column to West
        #print('West')
        c_x = np.where(SWAMPSlons+0.25 == lon+lon_res)[0][0]

    #print(c_x)
    # Add the extra half values:
    GC_wetl_area += np.sum(SWAMPS_wetlarea[:,ll_y:ur_y+1,c_x],axis=1)*0.5
    GC_area += np.sum(SWAMPSarea[ll_y:ur_y+1,c_x])*0.5

    # Calculate the Fwetl at imogen resolution
    GC_Fwetl = GC_wetl_area/GC_area
    
    #print(ipt,lat,lon,GC_wetl_area.max(),GC_area,GC_Fwetl.max())

    im_SWAMPS_1d[:,ipt] = GC_Fwetl 

#read in the 2d imogen data:
# Read in 2d grid info to check and store in 2d
inf_2d=nc.Dataset(area_file_2d,'r')
im_lat_2d = inf_2d.variables['latitude'][:]
im_lon_2d = inf_2d.variables['longitude'][:]
land_index = inf_2d.variables['land_index'][:]
inf_2d.close()

im_SWAMPS_wetlarea = im_SWAMPS_1d*im_area

mask_t2d = np.array([land_index.mask for ts in range(ntsteps)])
#
im_SWAMPS_2d = np.ma.masked_array(im_SWAMPS_1d[:,land_index],mask=mask_t2d)
im_SWAMPS_2d[im_SWAMPS_2d.mask==True]=fill_value
#
im_SWAMPS_wetlarea_2d = np.ma.masked_array(im_SWAMPS_wetlarea[:,land_index],mask=mask_t2d)
im_SWAMPS_wetlarea_2d[im_SWAMPS_wetlarea_2d.mask==True]=fill_value
#
fig,axes = plt.subplots(nrows=5,ncols=1,figsize=(8,15))
#
raw_extent = [-180.,180,-90,90]
raw_image=axes[0].imshow(SWAMPS_wetlarea.mean(axis=0),vmin=0,vmax=0.4e9,cmap='YlGn',extent=raw_extent)
#raw_image=axes[0].imshow(SWAMPSdata.mean(axis=0),vmin=0,vmax=0.4,cmap='YlGn',extent=raw_extent)
#raw_image=axes[0].imshow(SWAMPSdata[6,:],vmin=0,vmax=0.4,cmap='YlGn')
plt.colorbar(raw_image,ax=axes[0])
#
im_extent = [-180.,180,-55,85]
im_image=axes[1].imshow(im_SWAMPS_wetlarea_2d.mean(axis=0),vmin=0,vmax=0.4e10,cmap='YlGn',origin='bottom',extent=im_extent)
#im_image=axes[1].imshow(im_SWAMPS_1d[6,land_index],vmin=0,vmax=0.4,cmap='YlGn',origin='bottom')
plt.colorbar(im_image,ax=axes[1])
#
axes[2].plot(np.sum(SWAMPS_wetlarea.reshape(ntsteps,-1),axis=1),label='Raw')
axes[2].plot(np.sum(im_SWAMPS_wetlarea,axis=1),label='IMOGEN')
#axes[2].legend(ncol=2)#,loc='lower right')
#
axes[3].plot(np.mean(SWAMPS_wetlarea,axis=0).sum(axis=1)/latdel,SWAMPSlats,label='Raw')
axes[3].plot(np.mean(im_SWAMPS_wetlarea_2d,axis=0).sum(axis=1)/lat_res,im_lat_2d[:,0],label='IMOGEN',lw=2)
axes[3].legend(loc='lower right')
#
axes[4].plot(SWAMPSlons,np.mean(SWAMPS_wetlarea,axis=0).sum(axis=0)/londel,label='Raw')
axes[4].plot(im_lon_2d[0,:],np.mean(im_SWAMPS_wetlarea_2d,axis=0).sum(axis=0)/lon_res,label='IMOGEN',lw=2)
axes[4].legend()#loc='lower right')
#
plt.show()

print(im_SWAMPS_wetlarea_2d.shape)
#pdb.set_trace()
print(SWAMPS_outfile_2d)
outf=nc.Dataset(SWAMPS_outfile_2d,'w')
outf.createDimension('time',ntsteps)
outf.createDimension('latitude',im_lat_2d.shape[0])
outf.createDimension('longitude',im_lat_2d.shape[1])
outvar=outf.createVariable('time','float32',('time'))
outvar.units=SWAMPStimeunits
outvar.calendar=SWAMPStimecalendar
outvar[:]=SWAMPStime
outvar=outf.createVariable('latitude','float32',('latitude'))
outvar[:]=im_lat_2d[:,0]
outvar=outf.createVariable('longitude','float32',('longitude'))
outvar[:]=im_lon_2d[0,:]
outvar=outf.createVariable('Fw','float32',('time','latitude','longitude'),fill_value=fill_value)
outvar.units='1'
outvar.long_name='Fraction of wetland in gridcell'
outvar[:]=im_SWAMPS_2d
outf.Title='SWAMPS Fw data regridded onto the IMOGEN grid by Edward Comyn-Platt'
outf.note='A fraction of the global total is lost due to the coarse resolution of the IMOGEN grid'
outf.close()



print(SWAMPS_outfile)
outf=nc.Dataset(SWAMPS_outfile,'w')
outf.createDimension('time',ntsteps)
outf.createDimension('x',n_impoints)
outf.createDimension('y',1)
outvar=outf.createVariable('time','float32',('time'))
outvar.units=SWAMPStimeunits
outvar.calendar=SWAMPStimecalendar
outvar[:]=SWAMPStime
outvar=outf.createVariable('latitude','float32',('y','x'))
outvar[:]=im_lat
outvar=outf.createVariable('longitude','float32',('y','x'))
outvar[:]=im_lon
outvar=outf.createVariable('Fw','float32',('time','y','x'),fill_value=fill_value)
outvar.units='1'
outvar.long_name='Fraction of wetland in gridcell'
outvar[:]=im_SWAMPS_1d
outf.Title='SWAMPS Fw data regridded onto the IMOGEN grid by Edward Comyn-Platt'
outf.note='A fraction of the global total is lost due to the coarse resolution of the IMOGEN grid'
outf.close()







