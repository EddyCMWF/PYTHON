#!/usr/bin/env python

# 
import numpy as np
import netCDF4 as nc
import maths_tools.MetTools as MT
import maths_tools.FluxSiteTools as FST
import netcdftime as nctime

CRICH_dir='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Rob_Data/'
infile=CRICH_dir+'Crichton_preliminary.nc'

out_Drivefile=CRICH_dir+'Crichton_JULES_Drive_ECP.nc'
out_Fluxfile=CRICH_dir+'Crichton_Fluxes_ECP.nc'

site_lon=-3.583439
site_lat=55.042767

in_drive_vars= ['mean_incoming_flux_sortwave_radiation',\
                'mean_net_flux_total_radiation',\
                'mean_speed_horizontal_wind',\
                'standard_deviation_speed_horizontal_wind',\
                'mean_direction_horizontal_wind_from',\
                'mean_pressure_air',\
                'mean_temperature_air',\
                'mean_relative_humidity_air',\
                'mean_millimoles_per_cubic_meter_CO2', ]

nVARs=len(in_drive_vars)

out_drive_vars=['sw_down', \
                'rad_net', \
                'wind',    \
                'wind_std',\
                'wind_dir',\
                'pstar',   \
                't',       \
                'q',       \
                'co2',     \
                ]

out_drive_units=['W m^-2', \
                 'W m^-2', \
                 'm s^-1', \
                 'm s^-1', \
                 'degrees north', \
                 'Pa',\
                 'K',\
                 'kg kg^-1',\
                 'kg kg^-1',\
                 ]

out_drive_con_fact = [ lambda x: x, \
                       lambda x: x, \
                       lambda x: x, \
                       lambda x: x, \
                       lambda x: x, \
                       lambda x: x*1e3, \
                       lambda x: x+273.15, \
                       lambda x: MT.rh2sh(x,p,T), \
                       lambda x: MT.conc2mixr(x,p=p,T=T,Msp=44.01), \
                       ]
                        
in_flux_vars=['mean_incoming_flux_soil_heat', \
              'standard_deviation_incoming_flux_soil_heat', \
              'mean_outgoing_flux_sensible_heat', \
              'standard_deviation_outgoing_flux_sensible_heat',\
              'mean_outgoing_flux_latent_heat',\
              'standard_deviation_outgoing_flux_latent_heat',\
              'mean_outgoing_flux_CO2',\
              'standard_deviation_outgoing_flux_CO2']

nFVARs = len(in_flux_vars)

out_flux_vars=['soil_heat', \
               'soil_heat_std',\
               'sensible_heat',\
               'sensible_heat_std',\
               'latent_heat',\
               'latent_heat_std',\
               'NEE',\
               'NEE_std']

out_flux_units=['W m^-2', \
                'W m^-2', \
                'W m^-2', \
                'W m^-2', \
                'W m^-2', \
                'W m^-2', \
                'umol m^-2 s-1', \
                'umol m^-2 s-1', \
                 ]




inf=nc.Dataset(infile,'r')
infCrich = inf.groups['Crichton'].groups['Preliminary']
p=infCrich.variables['mean_pressure_air'][:]*1e3
T=infCrich.variables['mean_temperature_air'][:]+273.15
SWd=infCrich.variables['mean_incoming_flux_sortwave_radiation'][:]
NEE=infCrich.variables['mean_outgoing_flux_CO2'][:]

# Calculate the GPP
GPP, TER, GPPFitParams = \
        FST.GPP_from_NEE_SW_T(NEE,SWd,T,\
                               spike_filter='median', FIT_maxfev=1000, \
                               ReturnResp=True, \
                               ReturnFitParams=True)


outf_drive = nc.Dataset(out_Drivefile,'w')
# create dimensions
outf_drive.createDimension('time',len(infCrich.dimensions['time']))
outf_drive.createDimension('land',1)

outvar=outf_drive.createVariable('lat','float32',('land'))
outvar[:]=site_lat
outvar=outf_drive.createVariable('lon','float32',('land'))
outvar[:]=site_lon

outvar=outf_drive.createVariable('time','float64',('time'))
outvar.units='seconds since 1970-01-01 00:00:00'
outvar[:]= infCrich.variables['time'][:]

for iVAR in range(nVARs):
    print out_drive_vars[iVAR]
    outvar=outf_drive.createVariable( out_drive_vars[iVAR], 'float32', \
                                ('time','land') )
   
    outvar.units=out_drive_units[iVAR]
    outvar[:]=out_drive_con_fact[iVAR]( infCrich.variables[in_drive_vars[iVAR]][:] )

outf_drive.close()


outf_flux = nc.Dataset(out_Fluxfile,'w')
# create dimensions
outf_flux.createDimension('time',len(infCrich.dimensions['time']))
outf_flux.createDimension('land',1)

outvar=outf_flux.createVariable('lat','float32',('land'))
outvar[:]=site_lat
outvar=outf_flux.createVariable('lon','float32',('land'))
outvar[:]=site_lon

outvar=outf_flux.createVariable('time','float64',('time'))
outvar.units='seconds since 1970-01-01 00:00:00'
outvar[:]= infCrich.variables['time'][:]

for iVAR in range(nFVARs):
    print out_flux_vars[iVAR]
    outvar=outf_flux.createVariable( out_flux_vars[iVAR], 'float32', \
                                     ('time','land') )
   
    outvar.units=out_flux_units[iVAR]
    outvar[:]=infCrich.variables[in_flux_vars[iVAR]][:]

# output the GPP and TER
outvar=outf_flux.createVariable('GPP','float32',('time','land') )
outvar.units='umol m^-2 s-1'
outvar[:]=GPP
outvar=outf_flux.createVariable('TER','float32',('time','land') )
outvar.units='umol m^-2 s-1'
outvar[:]=TER


outf_flux.close()

inf.close()




