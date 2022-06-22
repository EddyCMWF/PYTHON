#!/usr/bin/python
#
# Script to fetch the metadata for each of the PALS sites
#    Such as start and end data and geographical coordinates
#
#
###########################################
import netCDF4 as nc
import netcdftime as nctime
import csv


sitename_file = '/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/PALS_comparison/sites.txt'
site_metadata_infile  = '/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/PALS_comparison/PALS_sites.csv'

site_metadata_outfile = '/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/PALS_comparison/sites_meta.csv'

DATA_DIR = '/data/grp/fluxdata/PALS_sites/'

sites = open(sitename_file,'r').readlines()

albmar_metadata = list(csv.reader(open(site_metadata_infile,'r'),delimiter=','))
metdat_hdrs=albmar_metadata.pop(0)

outf=open(site_metadata_outfile,'w')
outf.write( '#%19s, %20s, %20s, %15s, %15s, %20s, %20s, %10s, %10s, \n' %
            ('Site', \
             'Country', \
             'Vegetation Type', \
             'Canopy Hgt', \
             'Meas. Hgt', \
             'Start Date', \
             'End Date',   \
             'lat', 'lon')   )

for i in range(len(sites)):
    site    = albmar_metadata[i][0] 
    country = albmar_metadata[i][1]
    veg_type= albmar_metadata[i][2]
    canht   = albmar_metadata[i][5]
    measht  = albmar_metadata[i][6]
    lat     = albmar_metadata[i][7]
    lon     = albmar_metadata[i][8]

    file=site+'Fluxnet.1.4_flux_threehourly.nc'
    print DATA_DIR+file
    #print DATA_DIR+file
    inf=nc.Dataset(DATA_DIR+file,'r')
    time_obj=nctime.num2date( inf.variables['time'][:],                  \
                              units=inf.variables['time'].units,      \
                              calendar='standard' ) 
    
    outf.write( '%20s, %20s, %20s, %15s, %15s, %20s, %20s, %10s, %10s, \n' % \
                (site,               \
                 country,            \
                 veg_type,           \
                 str(canht),         \
                 str(measht),        \
                 str(time_obj[0]),   \
                 str(time_obj[-1]),  \
                 str(lat), str(lon), )  )
    
outf.close()


