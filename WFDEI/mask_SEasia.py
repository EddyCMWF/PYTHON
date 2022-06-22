#!/bin/env python
#
import netCDF4 as nc
import numpy as np

base_dir='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/MAY_2016_results/'

runs=['Jv4.5_WFDEI_nti_NG-HWSD_gridded_monthly_ch4',\
      'Jv4.5_WFDEI_nti_NG-HWSD_nounf_gridded_monthly_ch4',\
      'Jv4.5_WFDEI_nti_TRIFFID_gridded_monthly_ch4',\
      'Jv4.5_WFDEI_nti_TRIFFID_nounf_gridded_monthly_ch4',\
      ]

lat_lims=[0,45]
lon_lims=[60,150]
start_year,end_year=1980,2014
lat_name,lon_name='latitude','longitude'

for run in runs:
    for year in range(start_year,end_year+1):
        inf_name=base_dir+run+'.'+str(year)+'.nc'
        outf_name=base_dir+run+'_SEasiaMASKED'+str(year)+'.nc'
        inf=nc.Dataset(inf_name,'r')
        outf=nc.Dataset(outf_name,'w')

        lons_2d, lats_2d = np.meshgrid(inf.variables[lon_name][:], \
                                       inf.variables[lat_name][:]  )
        
        mask_index = (lats_2d>lat_lims[0])&(lats_2d<lat_lims[1]) & \
                     (lons_2d>lon_lims[0])&(lons_2d<lon_lims[1]) 
        
        for dim in inf.dimensions:
            outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))
        
        for var in inf.variables:
            invar=inf.variables[var]
            outvar=outf.createVariable(str(var),invar.dtype,invar.dimensions)
            for att in invar.ncattrs():
                outvar.setncattr(str(att),invar.getncattr(str(att)))
            if str(var) in ['fwetl','fch4_wetl']:
                outdata=invar[:]
                outdata[:,mask_index==True].data=0.
                outvar[:]=outdata
            else:
                outvar[:]=invar[:]
        
        for att in outf.ncattrs():
            outf.setncattr(str(att),inf.getncattr(str(att)))
        
        outf.note='SE Asia (60-150E, 0-45N) has been set to zero to prevent double counting with rice paddies'
        
        outf.close()
        inf.close()
        

