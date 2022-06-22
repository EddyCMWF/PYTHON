#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import csv
import netcdftime as nctime

Ndep_file='/scratch/emrobi/Ndeposition/accmip_ndep_acchist_interp_wfd_1850_2013.nc'

# read N deposition
Ndinf=nc.Dataset(Ndep_file,'r')
Ndep = Ndinf.variables['ndep'][:]
land = Ndinf.variables['land'][:]
Ndep_time=Ndinf.variables['time'][:]
Ndinf.close()

Ndep_2D = np.zeros([len(Ndep_time),360,720])-999.
for iTIME in range(len(Ndep_time)):
    Ndep_2D[iTIME,:].flat[land]=Ndep[iTIME,:]

Ndep_2D = np.ma.masked_equal(Ndep_2D,-999.)



PALS_dir='/users/eow/edwcom/GREENHOUSE/PALS_comparison/'

site_metadata_file=PALS_dir+'sites_meta.csv'
site_md_list = list(csv.reader(open(site_metadata_file,'r'),delimiter=','))
site_md_headers=site_md_list.pop(0)

for site_md in site_md_list:
    # strip relevant meta data of blank spaces
    site       = site_md[0].strip()
    Start_year = int(site_md[5].strip()[:4])
    End_year   = int(site_md[6].strip()[:4])
    latitude   = float(site_md[7].strip())
    longitude  = float(site_md[8].strip())
    if longitude<0:
        longitude=longitude+360.
    
    # convert Start_year and End_year into sp and ep
    sp=Start_year-1850
    ep=End_year-1850
    
    #convert lat and lon into i j points
    i = int(np.round(longitude*2))
    j = int(np.round((latitude+90)*2))
    
    print latitude, longitude
    print j,i
    
    pt_Ndep_ts = Ndep_2D[sp:ep+1,j,i]
        
    outf=open(PALS_dir+'ancil_data/site_depostion_files/'+site+'_Ndeposition_mean.dat','w')
    outf.write('# mean N deposition for '+site+'\n')
    outf.write(' %15.5e \n' % (np.mean(pt_Ndep_ts)) )
    outf.close()

