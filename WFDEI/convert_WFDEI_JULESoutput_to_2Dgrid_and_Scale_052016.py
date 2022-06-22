#!/usr/bin/python
#
# Python script to convert CRUNCEP JULES fch4_wetl output onto 2D grid with annual ncdf files
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# April 2015
#
# Contains
import numpy as np
import netCDF4 as nc
import netcdftime as nctime

gridfile='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
area_file='/users/eow/edwcom/WFD_EI/wfdei-land-area.nc'
frac_file='/users/eow/edwcom/WFD_EI/wfdei-land-frac.nc'

data_dir='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/MAY_2016_results/'
run='JULES_WFDEI_nti_TRIFFID'
profile='monthly_soil'

out_dir=data_dir+'gridded/'
out_run='J4.5_WFDEI-GPCC_nti_TRIFFID'
out_profile='monthly_gridded'

syear=1980
eyear=2014
scale_year=2000
scale_value=180e9 # 180 Tg = 180e9 kg
JULES_to_TgperYr_fact= 1e-9*86400*365*(16.01/12.01)

#params=['fch4_wetl','fwetl']
fch4_conv_factor=(16.01/12.01)*1e-9
time_name='time'
fill_value=-999.9
#

metdata_source='WFDEI at 0.5 degrees'
soil_carbon_source='TRIFFID Dynamic'
topography_source='Marthews (2014)'
soilproperties_source='Van Genuchten'


grinf=nc.Dataset(gridfile,'r')
index=grinf.variables['land_index'][:]
lats=grinf.variables['latitude'][:]
lons=grinf.variables['longitude'][:]
grinf.close()
Annual_Mask=np.array([index.mask for i in range(12)])

area_inf=nc.Dataset(area_file,'r')
AREA=area_inf.variables['area'][:]
area_inf.close()

frac_inf=nc.Dataset(frac_file,'r')
frac=frac_inf.variables['land_frac'][:]
frac_inf.close()

#Scale data
scale_file=run+'.'+profile+'.'+str(scale_year)+'.nc'
scale_inf=nc.Dataset(data_dir+scale_file,'r')
fch4_sc_y = scale_inf.variables['fch4_wetl'][:].squeeze()
scale_inf.close()

fch4_sc_y_abs=np.mean(fch4_sc_y * AREA * JULES_to_TgperYr_fact,axis=0)
fch4_sc_y_abs_2D=np.ma.masked_array(fch4_sc_y_abs[index-1],mask=index.mask)
fch4_sc_y_tot=np.sum(fch4_sc_y_abs)
fch4_scale_factor=float(np.round(scale_value/fch4_sc_y_tot,3))
#fch4_scale_factor= 180./300.    # Based on scaling the 2000 total global emissions to 180 Tg y-1


for year in range(syear,eyear+1):
    # Repeat for NTI data
    infile=run+'.'+profile+'.'+str(year)+'.nc'
    outfile=out_run+'.'+out_profile+'.'+str(year)+'.nc'
    inf=nc.Dataset(data_dir+infile,'r')
    outf=nc.Dataset(out_dir+outfile,'w')

    # Create Dimensions
    outf.createDimension('latitude',len(lats))
    outf.createDimension('longitude',len(lons))
    outf.createDimension('time',len(inf.dimensions['time']))

    # create dimension variables
    outvar=outf.createVariable('latitude', 'float32', ('latitude') )
    outvar.units='Degrees North'
    outvar[:]=lats

    outvar=outf.createVariable('longitude', 'float32', ('longitude') )
    outvar.units='Degrees East'
    outvar[:]=lons

    outvar=outf.createVariable('time', 'float32', ('time') )
    for attr in inf.variables['time'].ncattrs():
        outvar.setncattr( str(attr),inf.variables['time'].getncattr(str(attr)) )
    outvar[:]=inf.variables['time'][:]

    #Copy 'fwetl' straight over
    param = 'fwetl'
    outvar=outf.createVariable(param, inf.variables[param].dtype, ('time','latitude','longitude') )
    for attr in inf.variables[param].ncattrs():
        outvar.setncattr( str(attr),inf.variables[param].getncattr(str(attr)) )
    outdata=inf.variables[param][:].squeeze()
    outdata=np.ma.masked_array(outdata[:,index-1],mask=Annual_Mask,\
                                fill_value=inf.variables[param]._FillValue  )
    outdata.data[outdata.mask==True]=outdata.fill_value
    outvar[:]=outdata

    #Scale 'fch4_wetl' before copying
    param = 'fch4_wetl'
    outvar= outf.createVariable(param, inf.variables[param].dtype, ('time','latitude','longitude') )
    for attr in  inf.variables[param].ncattrs():
        if str(attr)=='units':
            outvar.setncattr( str(attr),'kg CH4 m^-2 s^-1')
        else:
            outvar.setncattr( str(attr),inf.variables[param].getncattr(str(attr)) )
    # Add extra note stating that the JULES output has been rescaled
    outvar.setncattr('Correction_Scale_Factor',str(fch4_scale_factor) )
    
    outdata=inf.variables[param][:].squeeze()
    outdata=np.ma.masked_array(outdata[:,index-1],mask=Annual_Mask,\
                                fill_value=inf.variables[param]._FillValue  )
    outdata.data[outdata.mask==True]=outdata.fill_value
    outvar[:]=outdata*fch4_conv_factor*fch4_scale_factor

    outf.title='JULES output: Methane emissions from wetlands'
    outf.owner='Edward Comyn-Platt, edwcom@ceh.ac.uk'
    outf.metdata=metdata_source
    outf.SCdata=soil_carbon_source
    outf.topography=topography_source
    outf.soilproperties=soilproperties_source
    outf.note='fch4_wetl has been scaled to give total annual emissions for 2000 ~180 Tg'
    outf.close()
    inf.close()


