#!/usr/bin/env python

import gdal
import numpy as np
import netCDF4 as nc

topo_file  = '/prj/wetlands_africa/topo/ti_marthews/topidx.tif'

WFDEI025_topo_outfile='/users/eow/edwcom/WFD_EI/topoidx_WFDEI_0p25_global.nc'

WFDEI_gfile = '/users/eow/edwcom/WFD_EI/start_dumps/J4.6_WFDEI_0.25_PRESCRIBED.dump.spin1.20000101.0.nc'

# Read in lat and lons from WFDEI dumpfile
Winf = nc.Dataset(WFDEI_gfile,'r')
Wlats = Winf.variables['latitude'][:]
Wlons = Winf.variables['longitude'][:]
Winf.close()

#W_latlon = zip(Wlats,Wlons)
nPTs = len(Wlats)

topidx_mean=np.zeros_like(Wlats)
topidx_std=np.zeros_like(Wlats)
fexp=np.zeros_like(Wlats)+3.0

inf=gdal.Open(topo_file, gdal.GA_ReadOnly)
orig_res=15./3600.   # 15 arcseconds

fill_value=-9999.

for iPT in range(nPTs):
    if (Wlats[iPT]<-56.5)|(Wlats[iPT]>86):
        topidx_mean[iPT]=fill_value
        topidx_std[iPT]=fill_value
    else:
        latmin=(Wlats[iPT]-0.125)
        latmax=(Wlats[iPT]+0.125)
        lonmin=(Wlons[iPT]-0.125)
        lonmax=(Wlons[iPT]+0.125)
        y_start=int((176-(latmax+90.))/orig_res)
        y_end  =int((176-(latmin+90.))/orig_res)
        x_start=int((lonmin+180.)/orig_res)
        x_end  =int((lonmax+180.)/orig_res)
        
        indata=inf.ReadAsArray( xoff=x_start,yoff=y_start, xsize=(x_end-x_start), ysize=(y_end-y_start))
        indata=np.ma.masked_array(indata,indata<-1,fill_value=fill_value)

        topidx_mean[iPT]=np.mean(indata)
        topidx_std[iPT]=np.std(indata)

    
topidx_mean=np.ma.masked_equal(topidx_mean,fill_value)    
topidx_std=np.ma.masked_equal(topidx_std,fill_value)    

# swap nan in C_topidx with closest neighbour
goodex=np.where(topidx_mean.mask==False)[0]
badex=np.where(topidx_mean.mask==True)[0]
print('len(goodex),len(badex)')
print(len(goodex),len(badex))

for bad_pt in badex:
    swapdex = np.argmin( ((Wlats[goodex]-Wlats[bad_pt])**2) + \
                         ((Wlons[goodex]-Wlons[bad_pt])**2)  )
    topidx_mean[bad_pt]=topidx_mean[swapdex]
    topidx_std[bad_pt]=topidx_std[swapdex]

# Output data to file
outf=nc.Dataset(WFDEI025_topo_outfile,'w')

outf.createDimension('land',len(Wlats))

#longitude
outvar                = outf.createVariable('longitude','float',('land'))
outvar.units          = "degrees east"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre longitude"
outvar[:]             = Wlons

#latitude
outvar                = outf.createVariable('latitude','float',('land'))
outvar.units          = "degrees north"
outvar.missing_values = -9999.0
outvar.longname       = "grid box centre latitude"
outvar[:]             = Wlats

#ti_mean
outvar                = outf.createVariable('ti_mean','float',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index mean"
outvar[:]             = topidx_mean

#ti_std
outvar                = outf.createVariable('ti_std','float',('land'))
outvar.units          = "-"
outvar.missing_values = -9999.0
outvar.longname       = "topogrgaphic index standard deviation"
outvar[:]             = topidx_std

outf.setncattr('title','Aggregated topographic index data for WFDEI-0.25degree')
outf.setncattr('institution','CEH - Wallingford')
outf.setncattr('source','Simon Dadson and Toby Marthews, U. Oxford')
outf.setncattr('owner', 'E. Comyn-Platt (edwcom@ceh.ac.uk)')
outf.setncattr('note','Extracted 15 arcsecond data provided by T. Marthews')

outf.close()



