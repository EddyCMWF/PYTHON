#!/usr/bin/python

import gdal
import numpy as np
import netCDF4 as nc

topo_dir   = '/prj/wetlands_africa/topo/ti_marthews/'

CRU_topo_outfile='/users/eow/edwcom/CRUNCEP/topoidx_CRUNCEP_0p5_lp_global.nc'

CRUNCEP_gfile = '/users/eow/edwcom/CRUNCEP/cru_ncep_land.nc'


# read in lons from topo long file
img = gdal.Open(topo_dir+'long_05x05_weighted.dat')
band= img.GetRasterBand(1)
tlons = band.ReadAsArray() + 0.25  #add on 0.25 for center point
del img
del band
# read in lons from topo lat file
img = gdal.Open(topo_dir+'lat_05x05_weighted.dat')
band= img.GetRasterBand(1)
tlats = band.ReadAsArray() - 0.25  # subtract 0.25 for center point
del img
del band
# read in lons from topo long file
img = gdal.Open(topo_dir+'topidx_05x05_weighted.dat')
band= img.GetRasterBand(1)
topidx = band.ReadAsArray() #add on 0.25 for center point
del img
del band
# read in lons from topoidx_SD file
img = gdal.Open(topo_dir+'topidx_SD_05x05_weighted.dat')
band= img.GetRasterBand(1)
topidx_SD = band.ReadAsArray() #add on 0.25 for center point
del img
del band

#zip topo lat and lons into list of tuples
topo_latlon = zip(tlats.flat,tlons.flat)



# Read in lat and lons from CRUNCEP land file
Cinf = nc.Dataset(CRUNCEP_gfile,'r')
clats = Cinf.variables['Latitude'][:]
clons = Cinf.variables['Longitude'][:]
cland = Cinf.variables['land'][:]
Cinf.close()

C_latlon = zip(clats,clons)

index = [ topo_latlon.index(pt) for pt in C_latlon ]

C_topidx    = topidx.flat[index]
C_topidx_SD = topidx_SD.flat[index]

C_fexp = np.ones_like(clats)

# swap nan in C_topidx with closest neighbour
goodex=np.where(np.isfinite(C_topidx))[0]
badex=np.where(np.isnan(C_topidx))[0]
for bad_pt in badex:
    swapdex = np.argmin( ((clats[goodex]-clats[bad_pt])**2) + \
                         ((clons[goodex]-clons[bad_pt])**2)  )
    C_topidx[bad_pt]=C_topidx[swapdex]

# swap nan in C_topidx_SD with closest neighbour
goodex=np.where(np.isfinite(C_topidx_SD))[0]
badex=np.where(np.isnan(C_topidx_SD))[0]
for bad_pt in badex:
    swapdex = np.argmin( ((clats[goodex]-clats[bad_pt])**2) + \
                         ((clons[goodex]-clons[bad_pt])**2)  )
    C_topidx_SD[bad_pt]=C_topidx_SD[swapdex]

Coutf=nc.Dataset(CRU_topo_outfile,'w')

Coutf.createDimension('land',len(clats))

#land index
outvar       = Coutf.createVariable('land','int',('land'))
outvar.units = "-"
outvar.note  = "Index for grid, row-wise from NW corner" 
outvar[:]    = cland

#longitude
outvar                = Coutf.createVariable('longitude','float',('land'))
outvar.units          = "degrees east"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre longitude"
outvar[:]             = clons

#latitude
outvar                = Coutf.createVariable('latitude','float',('land'))
outvar.units          = "degrees north"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre latitude"
outvar[:]             = clats

#ti_mean
outvar                = Coutf.createVariable('ti_mean','float',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index mean"
outvar[:]             = np.array(C_topidx)

#ti_std
outvar                = Coutf.createVariable('ti_std','float',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index standard deviation"
outvar[:]             = np.array(C_topidx_SD)

Coutf.setncattr('title','Aggregate topographic index data for CRU-NCEP')
Coutf.setncattr('institution','CEH - Wallingford')
Coutf.setncattr('source','Simon Dadson and Toby Marthews, U. Oxford')
Coutf.setncattr('contact', 'E. Comyn-Platt (edwcom@ceh.ac.uk)')
Coutf.setncattr('note','Extracted 0.5 degree data provided by T. Marthews')

Coutf.close()



