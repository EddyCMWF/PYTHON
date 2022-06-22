#!/usr/bin/env python
################################################################################
#
# Create a combined land mask between CEH-GEAR, CEH land cover (CHESS version), 
# IHDTM elevation on 1km grid
#
# ELR 2/10/2014
#
################################################################################

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors 
import datetime as dt
import mpl_toolkits.basemap.pyproj as pyproj 


osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

# http://spatialreference.org/ref/epsg/4258/
etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')

def bng_to_ll(x,y):
	lonlat = [ pyproj.transform(osgb36, etrs89, x[i], y[i] ) for i in range(len(x)) ]

	lon = np.array([l[0] for l in lonlat])
	lat = np.array([l[1] for l in lonlat])

	return lon, lat



def read_pdef(pdefname,nx,ny,bswap=False):

	# read the index (4byte ints)
	indx=np.fromfile(pdefname,dtype='int32',count=nx*ny).byteswap(bswap)
	indx=indx.reshape([ny,nx])

	# read the weights (floats) 
	wt=np.fromfile(pdefname,dtype='float32',count=-1).byteswap(bswap)
	wt=wt[nx*ny:(2*nx*ny)]
	wt=wt.reshape([ny,nx])

	return indx,wt

mv=-99999.0

################################################################################
# Read the land cover
f_lc='/prj/chess/data/1km/ancil/data/gb_lcov2000.gra'
np_lc=227964
nz_lc=8
lc_vector=np.fromfile(f_lc,dtype='float32').reshape([nz_lc,np_lc]).byteswap()

ny_lc=1044
nx_lc=600

# read the pdef
f_pdef_lc='/prj/chess/data/1km/ancil/data/gb_pdef.gra'
pdef_ind_lc, pdef_wt_lc = read_pdef(f_pdef_lc,nx_lc,ny_lc,True)
pdef_ind_lc=pdef_ind_lc.flatten()
where_pts = np.where(pdef_ind_lc!=-999)

# convert to python indices
pdef_ind_lc[where_pts]=pdef_ind_lc[where_pts]-1

# Put on full grid
lc=np.ones([nz_lc,ny_lc*nx_lc])*mv
for z in range(nz_lc):
	lc[z,where_pts]=lc_vector[z,pdef_ind_lc[where_pts]]

lc=lc.reshape([nz_lc,ny_lc,nx_lc])
lc=np.ma.masked_equal(lc,mv)

# put the land cover on the full 1300*700 grid
lc_full=np.ones([nz_lc,1300,700])*mv
lc_full[:,12:12+ny_lc,55:55+nx_lc]=lc[:]
lc_full=np.ma.masked_equal(lc_full,mv)

# Get the mask
lc_mask=lc_full[0,:,:].mask

################################################################################
# Read the original land cover
lc2_f=open('/prj/chess/emrobi/data/landcover2000/LCM2000GBSUBCLASSDOMINANT.asc','r')
ncols_lc2=int(lc2_f.readline().strip().split()[1])
nrows_lc2=int(lc2_f.readline().strip().split()[1])
xllcorner_lc2=float(lc2_f.readline().strip().split()[1])
yllcorner_lc2=float(lc2_f.readline().strip().split()[1])
cellsize_lc2=float(lc2_f.readline().strip().split()[1])
mv_lc2=float(lc2_f.readline().strip().split()[1])

lc2=[]
for lin in lc2_f.readlines():
	lc2.append([ float(l) for l in lin.strip().split() ])

lc2=np.ma.masked_equal(np.array(lc2)[::-1,:],mv_lc2)
seaval=1
litrocval=3
litsedval=4
lc2=np.ma.masked_equal(lc2,seaval)
lc2=np.ma.masked_equal(lc2,litrocval)
lc2=np.ma.masked_equal(lc2,litsedval)

# Get the mask
lc2_mask=lc2[:,:].mask

################################################################################
# Read the land cover from Jon
f_lc3='/users/rec/jon/JULESdrv/ancil/lcov2000_01.bin'

ny_lc3=1300
nx_lc3=700
mv_lc3=-99999

lc3=np.fromfile(f_lc3,dtype='float32').reshape([ny_lc3,nx_lc3]).byteswap()

lc3=np.ma.masked_equal(lc3,mv_lc3)


# Get the mask
lc3_mask=lc3[:,:].mask


################################################################################
# Read the elevation

f_elev='/local/localscratch/emrobi/CHESS/data/1km/v3/input_data/uk_ihdtm_elev_1km.bin'

ny_elev=1300
nx_elev=700

elev=np.ma.masked_equal(np.fromfile(f_elev,dtype='float32').reshape([ny_elev,nx_elev]).byteswap(),-99999.0)
elev=elev[::-1,:]

elev_mask=elev.mask

################################################################################
# read a GEAR file
f_gear='/prj/uk_rain/grids_generator/OUTPUT/permanent/GB/daily/netCDF/compressed/CEH_GEAR_daily_GB_1983.nc'

fg=nc.Dataset(f_gear,'r')

gear=fg.variables['rainfall_amount'][0,::-1,:]

gear_full=np.ones([1300,700])*gear.fill_value
gear_full[:1251,:]=gear[:,:700]
gear_full=np.ma.masked_equal(gear_full,gear.fill_value)

gear_mask=gear_full.mask

################################################################################
# read average wind speeds
f_ws='/prj/chess/code/metprod/input_data/wspeed_2m_avg.flt'
ny_ws=1300
nx_ws=700
ws=np.ma.masked_equal(np.fromfile(f_ws,dtype='float32').reshape([ny_ws,nx_ws]).byteswap(),-9999)
ws=ws[::-1,:]

ws_mask=ws.mask

################################################################################
# Extra masking for reasons of distance
# The scilly isles are mentioned in Jon's doc, but are excluded by the land
# cover map anyway, so no extra masking for them here
# but we do need to exclude the isle of man
x,y=np.meshgrid(np.arange(0,nx_lc3*1000,1000),np.arange(0,ny_lc3*1000,1000))

distmask1=x<50000
distmask2=y>1069000

#isleofman=np.logical_and(np.logical_and(x>210000,x<260000),np.logical_and(y>461000,y<510000))


################################################################################
# combined mask
# We use 'lc3', ie Jon's initial regridding to 1km, before masking any points

full_mask=~np.logical_and(np.logical_and(np.logical_and(~lc3_mask,~elev_mask),np.logical_and(~gear_mask,~ws_mask)),np.logical_and(~distmask1,~distmask2))

lc3_only = np.logical_and(~lc3_mask,full_mask)
n_lc3_only=len(np.where(lc3_only)[0])

elev_only = np.logical_and(~elev_mask,full_mask)
n_elev_only=len(np.where(elev_only)[0])

gear_only = np.logical_and(~gear_mask,full_mask)
n_gear_only=len(np.where(gear_only)[0])

ws_only = np.logical_and(~ws_mask,full_mask)
n_ws_only=len(np.where(ws_only)[0])

print "Number of points excluded from land cover map: %d"%n_lc3_only
print "Number of points excluded from elevation: %d"%n_elev_only
print "Number of points excluded from GEAR: %d"%n_gear_only
print "Number of points excluded from average wind speed: %d"%n_ws_only

fnum=0


fnum+=1
plt.figure(fnum)
plt.imshow(~full_mask[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='bone_r',vmin=0,vmax=5)
plt.imshow(np.ma.masked_where(~lc3_only,lc3_only)[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='hsv')
plt.title('Land cover map\n%d points excluded'%n_lc3_only)

fnum+=1
plt.figure(fnum)
plt.imshow(~full_mask[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='bone_r',vmin=0,vmax=5)
plt.imshow(np.ma.masked_where(~elev_only,elev_only)[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='hsv')
plt.title('Elevation\n%d points excluded'%n_elev_only)

fnum+=1
plt.figure(fnum)
plt.imshow(~full_mask[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='bone_r',vmin=0,vmax=5)
plt.imshow(np.ma.masked_where(~gear_only,gear_only)[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='hsv')
plt.title('GEAR\n%d points excluded'%n_gear_only)

fnum+=1
plt.figure(fnum)
plt.imshow(~full_mask[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='bone_r',vmin=0,vmax=5)
plt.imshow(np.ma.masked_where(~ws_only,ws_only)[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000],cmap='hsv')
plt.title('Wind speeds\n%d points excluded'%n_ws_only)

fnum+=1
plt.figure(fnum)
plt.imshow(np.ma.masked_where(full_mask,full_mask)[::-1,:],interpolation='nearest',extent=[0,700000,0,1300000])
plt.title('Full land mask')

fnum+=1
plt.figure(fnum)
lc_both=np.zeros([1300,700])
lc_both[np.where(np.logical_and(lc_mask,lc2_mask))]=mv
lc_both[np.where(np.logical_and(~lc_mask,~lc2_mask))]=3 
lc_both[np.where(np.logical_and(lc_mask,~lc2_mask))]=2 
lc_both[np.where(np.logical_and(~lc_mask,lc2_mask))]=1 
lc_both=np.ma.masked_equal(lc_both,mv)

newcmap=colors.ListedColormap(['r','b',[0.4,0.4,0.4]])
plt.imshow(lc_both[::-1,:],cmap=newcmap,interpolation='nearest',vmin=0.5,vmax=3.5)
cb=plt.colorbar(ticks=[1.0,2.0,3.0])
cb.set_ticklabels(['JF','orig','both'])

fnum+=1
plt.figure(fnum)
lc2_gear=np.zeros([1300,700])
lc2_gear[np.where(np.logical_and(gear_mask,lc2_mask))]=mv
lc2_gear[np.where(np.logical_and(~gear_mask,~lc2_mask))]=3 
lc2_gear[np.where(np.logical_and(gear_mask,~lc2_mask))]=2 
lc2_gear[np.where(np.logical_and(~gear_mask,lc2_mask))]=1 
lc2_gear=np.ma.masked_equal(lc2_gear,mv)

newcmap=colors.ListedColormap(['r','b',[0.4,0.4,0.4]])
plt.imshow(lc2_gear[::-1,:],cmap=newcmap,interpolation='nearest',vmin=0.5,vmax=3.5)
cb=plt.colorbar(ticks=[1.0,2.0,3.0])
cb.set_ticklabels(['GEAR','orig LCM2000','both'])

fnum+=1
plt.figure(fnum)
lc_both3=np.zeros([1300,700])
lc_both3[np.where(np.logical_and(lc_mask,lc3_mask))]=mv
lc_both3[np.where(np.logical_and(~lc_mask,~lc3_mask))]=3 
lc_both3[np.where(np.logical_and(lc_mask,~lc3_mask))]=2 
lc_both3[np.where(np.logical_and(~lc_mask,lc3_mask))]=1 
lc_both3=np.ma.masked_equal(lc_both3,mv)

newcmap=colors.ListedColormap(['r','b',[0.4,0.4,0.4]])
plt.imshow(lc_both3[::-1,:],cmap=newcmap,interpolation='nearest',vmin=0.5,vmax=3.5)
cb=plt.colorbar(ticks=[1.0,2.0,3.0])
cb.set_ticklabels(['JF','JF orig','both'])

fnum+=1
plt.figure(fnum)
lc3_gear=np.zeros([1300,700])
lc3_gear[np.where(np.logical_and(gear_mask,lc3_mask))]=mv
lc3_gear[np.where(np.logical_and(~gear_mask,~lc3_mask))]=3 
lc3_gear[np.where(np.logical_and(gear_mask,~lc3_mask))]=2 
lc3_gear[np.where(np.logical_and(~gear_mask,lc3_mask))]=1 
lc3_gear=np.ma.masked_equal(lc3_gear,mv)

newcmap=colors.ListedColormap(['r','b',[0.4,0.4,0.4]])
plt.imshow(lc3_gear[::-1,:],cmap=newcmap,interpolation='nearest',vmin=0.5,vmax=3.5)
cb=plt.colorbar(ticks=[1.0,2.0,3.0])
cb.set_ticklabels(['GEAR','JF orig LCM2000','both'])


################################################################################
# create output file
landmask=np.zeros([1300,700])
landmask[np.where(~full_mask)]=1.0

landmask.astype('float32').byteswap().tofile('/local/localscratch/emrobi/CHESS/data/1km/v3/input_data/gb_landmask_1km_v3.gra')

################################################################################
# Also create the land mask netCDF while we're at it
outf=nc.Dataset('/local/localscratch/emrobi/CHESS/data/1km/v3/ancil/chess_landfrac.nc','w')

nync=1057
nxnc=656

outf.createDimension('y',nync)
outf.createDimension('x',nxnc)

outf.createVariable('y','f',('y',),fill_value=mv)
outf.variables['y'].units='m'
outf.variables['y'].long_name = "northing - OSGB36 grid reference" 
outf.variables['y'].standard_name = "projection_y_coordinate" 
outf.variables['y'].units='m'
outf.variables['y'].point_spacing = "even" 
outf.variables['y'][:]=np.arange(0,nync*1000,1000)+500.0

outf.createVariable('x','f',('x',),fill_value=mv)
outf.variables['x'].units='m'
outf.variables['x'].long_name = "easting - OSGB36 grid reference" 
outf.variables['x'].standard_name = "projection_x_coordinate" 
outf.variables['x'].units='m'
outf.variables['x'].point_spacing = "even" 
outf.variables['x'][:]=np.arange(0,nxnc*1000,1000)+500.0

outf.createVariable('lat','f',('y','x'),fill_value=mv)
outf.variables['lat'].long_name = "latitude of grid box centre" 
outf.variables['lat'].standard_name = "latitude" 
outf.variables['lat'].units = "degrees_north" 

outf.createVariable('lon','f',('y','x'),fill_value=mv)
outf.variables['lon'].long_name = "longitude of grid box centre" 
outf.variables['lon'].standard_name = "longitude" 
outf.variables['lon'].units = "degrees_east" 

xg,yg=np.meshgrid(np.arange(0,nxnc*1000,1000),np.arange(0,nync*1000,1000))
xg+=500
yg+=500

glon,glat = bng_to_ll(xg.flatten(),yg.flatten())
glon=glon.reshape([nync,nxnc])
glat=glat.reshape([nync,nxnc])

outf.variables['lat'][:]=glat[:]
outf.variables['lon'][:]=glon[:]

# Coord variable
outf.createVariable('crs','i',[])
outf.variables['crs'].long_name = "coordinate_reference_system" 
outf.variables['crs'].grid_mapping_name = "transverse_mercator" 
outf.variables['crs'].EPSG_code = "EPSG:27700" 
outf.variables['crs'].semi_major_axis = 6377563.396 
outf.variables['crs'].semi_minor_axis = 6356256.91 
outf.variables['crs'].inverse_flattening = 299.3249646 
outf.variables['crs'].latitude_of_projection_origin = 49. 
outf.variables['crs'].longitude_of_projection_origin = -2. 
outf.variables['crs'].false_easting = 400000. 
outf.variables['crs'].false_northing = -100000. 
outf.variables['crs'].scale_factor_at_projection_origin = 0.9996012717 

# data variable
outf.createVariable('landfrac','f',('y','x'),fill_value=mv,zlib=True,complevel=6,shuffle=True)
outf.variables['landfrac'].standard_name='land_fraction'
outf.variables['landfrac'].long_name='Land fraction'
outf.variables['landfrac'].units='none'
outf.variables['landfrac'].grid_mapping='crs'
outf.variables['landfrac'][:]=landmask[:nync,:nxnc]

outf.title = "Land fraction" 
outf.description = "Land fraction at 1km resolution over Great Britain" 
outf.institution = "CEH Wallingford - NERC" 
now=dt.datetime.now()
outf.history = "Created "+now.strftime('%Y-%M-%d %h:%m:%s')
outf.date_created = now.strftime('%Y-%M-%d')


outf.geospatial_lat_min = glat.min()
outf.geospatial_lat_max = glat.max()
outf.geospatial_lon_min = glon.min()
outf.geospatial_lon_max = glon.max()
outf.version = "beta_version" 
outf.Conventions = "CF-1.6" 
outf.creator_name = "Emma L. Robinson" 
outf.creator_email = "emrobi@ceh.ac.uk" 

outf.close()
#plt.show()
