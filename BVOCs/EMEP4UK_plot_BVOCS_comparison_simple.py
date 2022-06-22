#!/usr/bin/python
#
# Edward Comyn-Platt (edwcom@ceh.ac.uk)
# 
###################################################
#

import os,sys
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import data_info_EMEP4UK as data_info
import plot_tools as PT

INTERACTIVE='Y' ##sys.argv[1] #'Y'
OUT_DIR='/users/eow/edwcom/EMEP/EMEP4UK/plots/JULES_EMEP4UK_comparison/'
iso_mean_range=[0,2]
iso_diff_range=[-1,1]
terp_mean_range=[0,4]
terp_diff_range=[-1,1]
lonrange=[-13,11]
latrange=[51.5,56.5]


if INTERACTIVE=='Y':
    iDISPLAY=raw_input('Display images? (Y/N) ')
else:
    iDISPLAY=sys.argv[2]


BVOC_SOURCES = data_info.BVOC_sources()
nSOURCES = len(BVOC_SOURCES)

# Read in data
# SOURCE[0] = JULES_BC
for year in range(2001,2014):


print 'Available sources: '
for iSOURCE in range(nSOURCES):
    print iSOURCE, BVOC_SOURCES[iSOURCE]['NAME']
if INTERACTIVE=='Y':
    SOURCE_NUMS = raw_input('Select sources seperated by commas, for all sources enter ALL: ')
    if SOURCE_NUMS=='ALL':
        SOURCE_NUMS=range(nSOURCES)
    else:
        SOURCE_NUMS=SOURCE_NUMS.split(',')
        SOURCE_NUMS = [ int(SOURCE_N) for SOURCE_N in SOURCE_NUMS ]
else:
    SOURCE_NUMS = sys.argv[3]
print SOURCE_NUMS



print 'Available time resolutions:'
Tresolutions=['hourly','daily','monthly']
TRES_list   = [3600,86400,-1]
for iTRES in range(len(Tresolutions)):
    print iTRES, Tresolutions[iTRES]

if INTERACTIVE=='Y':
    iTRES = input('Select a time resolutionfor comparison: ')
else:
    iTRES = 1
TRES = TRES_list[iTRES]
print TRES

iso_data=[]
terp_data=[]
time_data=[]

for iSOURCE in SOURCE_NUMS:
    # Read isoprene data onto 2D grid and convert to same units using special funciton:
    iso_temp,terp_temp,time_vector = \
                    data_info.read_netCDF_BVOC_data_to_TRES( BVOC_SOURCES[iSOURCE],\
                                                             TRES=TRES, \
                                                             return_time_vector=True )

    iso_data.append(iso_temp)
    terp_data.append(terp_data)
    time_data.append(time_vector)

# read lat/lon data from index file
grid_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc'
grinf=nc.Dataset(grid_file,'r')
lats_2D=grinf.variables['lats'][:]
lons_2D=grinf.variables['lons'][:]
grinf.close()

# Plot map of mean emission rate
if PLOTS[0]=='Y':
    for jSOURCE in range(len(SOURCE_NUMS)):
        iSOURCE=SOURCE_NUMS[jSOURCE]
        print 'Producing maps of mean emission rates for: '+BVOC_SOURCES[iSOURCE]['NAME']
        start_year=BVOC_SOURCES[iSOURCE]['start_year']
        end_year=BVOC_SOURCES[iSOURCE]['end_year']
        plot_map_data = np.mean(iso_data[jSOURCE],axis=0)
        data_range=iso_mean_range
        plot_title='Mean Isoprene Emission Rate - '   \
                    +  BVOC_SOURCES[iSOURCE]['NAME']  \
                    + ' - '+str(start_year)+','+str(end_year)
        plot_title=plot_title.replace('_','-')
        file_name=OUT_DIR+'Isoprene_Mean_Emission_Map_'      \
                   +BVOC_SOURCES[iSOURCE]['NAME']+'_'        \
                   +str(start_year)+'_'+str(end_year)+'.png'
        cbar_title= '$g.m^{-2}.yr^{-1}$'
        
        PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                    DATA_RANGE=data_range,                          \
                    COLOURS=['beige','greenyellow','darkgreen'],    \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,            \
                    CBAR_ORIENTATION='vertical',                    \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                    PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,  \
                    iDISPLAY='Y', FILE_PLOT=file_name,                           \
                    LATDEL=2., LONDEL=2., RESOLUTION='i',                       \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)

