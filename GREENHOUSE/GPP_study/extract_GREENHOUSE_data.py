#!/usr/bin/python2.7
import os
import netCDF4 as nc
import numpy as np
import netcdftime as nctime
import csv
#
#####################################################################################
#
GH_dir     = '/group_workspaces/jasmin/greenhouse/GREENHOUSE/models/datasets/'
pheno_file = GH_dir+'downscaled_pheno_drivers.csv'
met_file   = GH_dir+'downscaled_met_drivers.csv'
#
outdir     = '/group_workspaces/jasmin/greenhouse/GREENHOUSE/ECP_datasets/'
#
# read in pheno data
pheno_csv=open(pheno_file,'r')
data=list(csv.reader(pheno_csv,delimiter=','))
pheno_hdrs=data.pop(0)
pheno_dict={}
for hdr in pheno_hdrs:
    pheno_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(pheno_hdrs,dat_lin):
        pheno_dict[hdr].append(float(dat))
#
pheno_csv.close()
del data
#
# read in the met data    
met_csv=open(met_file,'r')
data=list(csv.reader(met_csv,delimiter=','))
met_hdrs=data.pop(0)
met_dict={}
for hdr in met_hdrs:
    met_dict[hdr]=[]
#
for dat_lin in data:
    for hdr,dat in zip(met_hdrs,dat_lin):
        met_dict[hdr].append(float(dat))
#
met_csv.close()
del data
#
# Create list of FULL_hdrs, i.e. all uniques params
FULL_hdrs=list(pheno_hdrs)
for hdr in met_hdrs:
    if not hdr in FULL_hdrs:
        FULL_hdrs.append(hdr)
#
# get list of lat,lon pairs (tuples)#
# get a unique 'set' of lat_lon pairs
unique_latlons=list(set(lat_lons))
#
# 
# Initiate full, site by site dictionary to store data.
FULL_DICT = {}
#
# for each site, initiate a new dictionary within FULL_DICT
for site_num in range(len(unique_latlons)):
    FULL_DICT['site_'+str(site_num)]={}
    #within each sub dictionary initiate list for each parameter
    for hdr in FULL_hdrs:
        FULL_DICT['site_'+str(site_num)][hdr]=[]
    # store lat and long as scaler vals
    FULL_DICT['site_'+str(site_num)]['lat']=[unique_latlons[site_num][0]]
    FULL_DICT['site_'+str(site_num)]['long']=[unique_latlons[site_num][1]]
#
#
# Now loop through each latlon pair, find site number and 
# store data in appropriate dictionary
for pt_num in range(len(lat_lons)):
    site = lat_lons[pt_num]
    site_num = unique_latlons.index(site)
    #
    # Loop through pheno data and extract to FULL_DICT
    for hdr in pheno_hdrs:
        # do not need to store lat long each time.
        if not hdr in ['lat','long']:
            FULL_DICT['site_'+str(site_num)][hdr].append( \
                    pheno_dict[hdr][pt_num] )
    #
    # Loop through met fields and extract hourly data 
    #   to appropriate dictionary
    # calculate met start and end points
    met_sp = long(pt_num)*24
    met_ep = met_sp+24
    for hdr in met_hdrs:
        if not hdr in pheno_hdrs:
            FULL_DICT['site_'+str(site_num)][hdr].append( \
                    met_dict[hdr][met_sp:met_ep])
    #
# Now extract soil params from WFD-EI ancils
#
# First, calculate indices of data based on 2D 0.5 degree grid
Uniq_LL_np = np.array(unique_latlons)
Uniq_lats  = Uniq_LL_np[:,0]
Uniq_lons  = Uniq_LL_np[:,1]
#
x_index    = np.floor((Uniq_lons+180.)*2).astype('int64')
y_index    = np.floor((Uniq_lats+90.)*2).astype('int64')
#
# open soil file and read in data
vg_soil_file_in='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.soil_HWSD_class3_van_genuchten2d.nc'
#
# open file
soil_inf=nc.Dataset(vg_soil_file_in,'r')
# read params using tags from CAP, replace with sensible names in my output # # 
soil_b       = soil_inf.variables['field1381'][y_index,x_index]
soil_sathh   = soil_inf.variables['field342'][y_index,x_index]
soil_satcon  = soil_inf.variables['field333'][y_index,x_index]
soil_sm_sat  = soil_inf.variables['field332'][y_index,x_index]
soil_sm_crit = soil_inf.variables['field330'][y_index,x_index]
soil_sm_wilt = soil_inf.variables['field329'][y_index,x_index]
soil_hcap    = soil_inf.variables['field335'][y_index,x_index]
soil_hcon    = soil_inf.variables['field336'][y_index,x_index]
soil_albsoil = soil_inf.variables['field1395'][y_index,x_index]
# close file
soil_inf.close()
#
# add soil data to Dictionary
for site_num in range(len(unique_latlons)):
    FULL_DICT['site_'+str(site_num)]['soil_b']      = [soil_b[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_sathh']  = [soil_sathh[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_satcon'] = [soil_satcon[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_sm_sat'] = [soil_sm_sat[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_sm_crit']= [soil_sm_crit[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_sm_wilt']= [soil_sm_wilt[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_hcap']   = [soil_hcap[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_hcon']   = [soil_hcon[site_num]]
    FULL_DICT['site_'+str(site_num)]['soil_albsoil']= [soil_albsoil[site_num]]


# Repeat for topo data
topo_file_in='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/topoidx_WFDEI_0p5_2D_global.nc'
# open file
topo_inf=nc.Dataset(topo_file_in,'r')
# extract relevant points
ti_mean = topo_inf.variables['ti_mean'][y_index,x_index]
ti_std  = topo_inf.variables['ti_std'][y_index,x_index]
fexp    = topo_inf.variables['fexp'][y_index,x_index]
# close file
topo_inf.close()
#
# add topo to Dictionary
for site_num in range(len(unique_latlons)):
    FULL_DICT['site_'+str(site_num)]['ti_mean'] = [ti_mean[site_num]]
    FULL_DICT['site_'+str(site_num)]['ti_std']  = [ti_std[site_num]]
    FULL_DICT['site_'+str(site_num)]['fexp']    = [fexp[site_num]]

# repeat for land cover data
land_cover_file = '/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.veg.frac2d.nc'
# open file
LC_inf = nc.Dataset(land_cover_file,'r')
Land_Frac = LC_inf.variables['field1391'][:,y_index,x_index]
#close file
LC_inf.close()
# Create array of binary land cover, i.e. 1 for max tile 0 for rest
MAX_land_frac = np.zeros_like(Land_Frac)
MAX_land_frac[np.argmax(Land_Frac,axis=0),np.arange(len(unique_latlons))] = 1.
# Create array of binary PFT, i.e. 1 for max PFT tile 0 for all other tiles
MAX_PFT       = np.zeros_like(Land_Frac)
MAX_PFT[np.argmax(Land_Frac[:5,:],axis=0),np.arange(len(unique_latlons))] = 1.
# surface types for reference:
surface_types= ['BF, NT, C3, C4, shrub, urban, lake, soil, ice']
#
# add data to dictionary
for site_num in range(len(unique_latlons)):
    FULL_DICT['site_'+str(site_num)]['Land_Frac']    = Land_Frac[:,site_num]
    FULL_DICT['site_'+str(site_num)]['Max_surftype'] = [MAX_land_frac[:,site_num]]
    FULL_DICT['site_'+str(site_num)]['Max_PFT']      = [MAX_PFT[:,site_num]]
    #FULL_DICT['site_'+str(site_num)]['surf_types']   = surface_types
#
#
# Repeat for LAI data
lai_file='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.veg.func2d.nc'
# open file
LAI_inf = nc.Dataset(lai_file,'r')
# extract data
LAI_Jules = LAI_inf.variables['field1392'][:,y_index,x_index]
CanH_Jules= LAI_inf.variables['field1393'][:,y_index,x_index]
# close file
LAI_inf.close()
#
# add to DICT
for site_num in range(len(unique_latlons)):
    FULL_DICT['site_'+str(site_num)]['LAI_CAP_Jules']  = LAI_Jules[:,site_num]
    FULL_DICT['site_'+str(site_num)]['CanH_Jules'] = CanH_Jules[:,site_num]
#
