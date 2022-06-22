#!/bin/env python2.7

import numpy as np
import netCDF4 as nc
import sys,os
import pandas as pd
#from imogen import data_info
#from PlotTools import plot_tools as PTs
import data_info
import plot_tools as PTs
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as col

import pdb

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

GtC_to_ppm=0.471
ppm_to_kgC = 1e12/GtC_to_ppm

PRESENT_DAY_YEAR=2015

Tile_names=data_info.TILE_short_names()
Tile_colours = data_info.TILE_colours()
nTiles=len(Tile_names)

# Configuration to compare, CONFIG 1 is plotted 1st, so ideally should be larger 
CONFIG1 = optional_argparse('-config1', 'BASELINE_CONFIG_NoLULUC') 
RUNID1  = optional_argparse('-runid1',  'BL' )
#CONFIG2 = optional_argparse('-config2', 'METHANE_FEEDBACK_NoLULUC')
#RUNID2  = optional_argparse('-runid2',  'CH4_FB')
#CONFIGS = [CONFIG1,CONFIG2]
CONFIGS = [CONFIG1] 
RUNIDS  = [RUNID1]
#RUNIDS  = [RUNID1,RUNID2]
nCONFIGS = len(CONFIGS)
print('Configs: ',CONFIGS)

# Directories containing JULES output and plot output directories:i
DATA_DIR = optional_argparse('-datadir','/prj/CLIFFTOP/ECP_output/')
print('DATA_DIR: '+DATA_DIR)
#DATA_DIR = '/group_workspaces/jasmin2/clifftop/CLIFFTOP/ECP_output/EQUILIBRIUM_OUTPUT/'+CONFIG+'/'

PLOT_TAG = optional_argparse('-plottag', 'BL_NoLULUC_Permafrost')
PLOT_DIR = optional_argparse('-plotdir',DATA_DIR+'plots/'+PLOT_TAG+'/')
print('PLOT_DIR: '+PLOT_DIR)
os.system('mkdir -p '+PLOT_DIR)

# Directory of Ocean Uptake data:
OCEAN_UPTAKE_DIR = '/prj/CLIFFTOP/toy_jules/Ocean_Draw_Down/'
OCEAN_DICT = { config: { '1p5deg':np.load(OCEAN_UPTAKE_DIR+config+'/ocean_uptake_accum.npy')[:,0,:]*-1e12,
                         '1p81p5deg':np.load(OCEAN_UPTAKE_DIR+config+'/ocean_uptake_accum.npy')[:,1,:]*-1e12,
                         '2deg':np.load(OCEAN_UPTAKE_DIR+config+'/ocean_uptake_accum.npy')[:,2,:]*-1e12,
                        }
                    for config in CONFIGS           }
OCEAN_START_YEAR=1850
print("Ocean Uptake data from: " + OCEAN_UPTAKE_DIR)

# Directory of ancillary data:
ANCILS_DIR='/prj/CLIFFTOP/COMMON_DATA/ANCILS/'
#ANCILS_DIR='./'
# Grid File (My index for converting the 1D jules output to a 2D grid)
GRID_file= ANCILS_DIR+'grid_info.nc'
grinf=nc.Dataset(GRID_file,'r')
lats_2din = grinf.variables['latitude'][:]
lons_2din = grinf.variables['longitude'][:]
#Area_2d = grinf.variables['Area'][:]       # I don't actually use this but it's here
land_index = grinf.variables['land_index'][:]
grinf.close()

lats_2d = np.zeros( [lats_2din.shape[0]+1,lats_2din.shape[1]+1 ] )
lats_2d[:lats_2din.shape[0],:lats_2din.shape[1]]=lats_2din
lats_2d[:-1,-1]=lats_2din[:,-1]
lats_2d[-1,:]=85.

lons_2d = np.zeros( [lons_2din.shape[0]+1,lons_2din.shape[1]+1 ] )
lons_2d[:lons_2din.shape[0],:lons_2din.shape[1]]=lons_2din
lons_2d[-1,:-1]=lons_2din[-1,:]
lons_2d[:,-1]=180.
#pdb.set_trace()

# 1Dimension grid cell area data for calculating totals etc.
AREA_file=ANCILS_DIR+'Area_in_iris_format.nc'
Ainf=nc.Dataset(AREA_file,'r')
AREA_1D = Ainf.variables['area'][:].squeeze()
lats_1D = Ainf.variables['latitude'][:].squeeze()
lons_1D = Ainf.variables['longitude'][:].squeeze()
Ainf.close()
#print(AREA_file)

# Select GCMs to plot:
GCMs=data_info.GCMs()    # This returns a list of all GCM names
#GCMs.pop(GCMs.index('CEN_IPSL_MOD_IPSL-CM5A-LR'))
#GCMs = ['CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0','CEN_MOHC_MOD_HadGEM2-ES','CEN_NOAA-GFDL_MOD_GFDL-ESM2G']
good_GCMs=[ 0,1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20,25,26,27,28,29,30 ]
GCMs = [ GCMs[i] for i in good_GCMs ] 
nGCMs = len(GCMs)
for i in range(nGCMs):
    print(i,GCMs[i])
#quit()
# Select Scenarios to plot:
SCENARIOs = ['presentday','1p5deg', '1p81p5deg' , '2deg' ]  # tag in the JULES output file directory 
SCENARIOyears = [PRESENT_DAY_YEAR,2099,2099,2099]  # tag in the JULES output file directory 
                                              #  (see filename construction later)
SCENARIO_names=['Present Day ('+str(PRESENT_DAY_YEAR)+')','1.5$^o$C (2100)','1.5$^o$C Overshoot (2100)','2.0$^o$C (2100)']  # Name to appear on plots etc.
nSCENARIOs=len(SCENARIOs)

# File containing the pre-industrial conditions to use as a baseline:
PI_year=1850
PI_scen='1p5deg'
PI_gcm='CEN_MOHC_MOD_HadGEM2-ES'
PI_config=CONFIG1
PI_runid=RUNID1

###################################################################################################
# Soil layer thicknesses - For calculation of permafrost depths 
dz_soil= np.array([0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293,
                   0.21517585,0.23784142,0.25980762,0.28117066,0.30200527,
                   0.32237098,0.34231625,0.36188121])
z_soil = np.cumsum(dz_soil)
#dweight_soil = dz_soil/np.sum(dz_soil)
###################################################################################################
# Names of variables to store global total stock in dictionary:
stock_vars = ['CV','CS','AtmCO2_ppm','AtmCO2_kg','frac','Woody_Products','OceanCO2',
              'PermaCS','OldCS','PermaArea0m','PermaArea1m','PermaArea2m','PermaArea3m']

# Variables to store map data:
map_vars = ['CV','CS','Max_Frac','PermaCS','OldCS']+[Tname for Tname in Tile_names ] +\
           ['Permafrost0m','Permafrost1m','Permafrost2m','Permafrost3m']

###################################################################################################
# Read in the Scenario data

# Dictionaries for storing data:
DATA_DICT= { config: { scenario:{var: [] for var in stock_vars} for scenario in SCENARIOs } for config in CONFIGS }
MAPDATA_DICT= { config: { scenario:{var: [] for var in map_vars} for scenario in SCENARIOs } for config in CONFIGS }

for iconfig in range(nCONFIGS):
  config=CONFIGS[iconfig]
  runid=RUNIDS[iconfig]
  for iscenario in range(nSCENARIOs):
    scenario=SCENARIOs[iscenario]
    scenyear=SCENARIOyears[iscenario]
    if scenario=='presentday':
        scentag='1p5deg'
    else:
        scentag=scenario

    #print('Scenario: ',scenario)
    for igcm in range(nGCMs):
        gcm=GCMs[igcm]
        #print('GCM: ',gcm)
        DUMP_FILE=DATA_DIR+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+scentag+'.dump.'+str(scenyear)+'0101.0.nc'
        #print(DUMP_FILE)
        DINF = nc.Dataset(DUMP_FILE,'r')
        Ann_File=DATA_DIR+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+scentag+'.Annual_carbon.'+str(scenyear)+'.nc'
        #print(Ann_File)
        Ainf=nc.Dataset(Ann_File,'r')
        TMon_File=DATA_DIR+config+'/'+gcm+'/'+runid+'_'+gcm+'_'+scentag+'.Monthly_h2o_t.'+str(scenyear)+'.nc'
        #print(TMon_File)
        TMinf=nc.Dataset(TMon_File,'r')
        # Read in Vegetation Carbon which  is just on land points:
        #CV = DINF.variables['cv'][:]
        CV = Ainf.variables['cv'][:].squeeze()
        MAPDATA_DICT[config][scenario]['CV'].append(CV)
        # Sum CV over land points
        CV = np.sum(CV*AREA_1D)
        #print('CV = ',CV)
        DATA_DICT[config][scenario]['CV'].append(CV)
        CS = Ainf.variables['cs_gb'][:].squeeze()
        # Store map data:
        MAPDATA_DICT[config][scenario]['CS'].append(CS)
        CS = np.sum( CS*AREA_1D )
        #print('CS = ', CS)
        DATA_DICT[config][scenario]['CS'].append(CS)
        
        # Permafrost stuff:
        FRAC_old_cs = DINF.variables['frac_cs_old_pool_gb'][:]
        CS_pool_layer = Ainf.variables['cs'][:].squeeze()
        T_soil = TMinf.variables['t_soil'][:].squeeze().max(axis=0)
        
        OldC = np.sum( np.sum(CS_pool_layer*FRAC_old_cs, axis=0), axis=0)
        MAPDATA_DICT[config][scenario]['OldCS'].append( np.copy(OldC) )
        DATA_DICT[config][scenario]['OldCS'].append( np.sum(OldC*AREA_1D) )

        PermaMask = np.zeros_like(T_soil)
        PermaMask[T_soil<273.15] = 1.0
        OldCS_layer = np.sum(CS_pool_layer*FRAC_old_cs,axis=0)
        #print(OldCS_layer.shape)
        #print(PermaMask.shape)
        PermaC = np.sum( OldCS_layer*PermaMask, axis=0 )
        #print(PermaC.shape)
        #pdb.set_trace()
        MAPDATA_DICT[config][scenario]['PermaCS'].append(np.copy(PermaC))
        DATA_DICT[config][scenario]['PermaCS'].append( np.sum(PermaC*AREA_1D) )
        
        Tmax_0m = T_soil.max(axis=0)
        permafrost_0m = np.zeros_like(Tmax_0m)
        permafrost_0m[Tmax_0m<273.15]=1.
        MAPDATA_DICT[config][scenario]['Permafrost0m'].append(np.copy(permafrost_0m))
        DATA_DICT[config][scenario]['PermaArea0m'].append( np.sum(permafrost_0m*AREA_1D) )

        Tmax_1m = T_soil[z_soil>1.,:].max(axis=0)
        permafrost_1m = np.zeros_like(Tmax_1m)
        permafrost_1m[Tmax_1m<273.15]=1.
        MAPDATA_DICT[config][scenario]['Permafrost1m'].append(np.copy(permafrost_1m))
        DATA_DICT[config][scenario]['PermaArea1m'].append( np.sum(permafrost_1m*AREA_1D) )

        Tmax_2m = T_soil[z_soil>2.,:].max(axis=0)
        permafrost_2m = np.zeros_like(Tmax_2m)
        permafrost_2m[Tmax_2m<273.15]=1.
        MAPDATA_DICT[config][scenario]['Permafrost2m'].append(np.copy(permafrost_2m))
        DATA_DICT[config][scenario]['PermaArea2m'].append( np.sum(permafrost_2m*AREA_1D) )

        Tmax_3m = T_soil[-1,:]
        permafrost_3m = np.zeros_like(Tmax_3m)
        permafrost_3m[Tmax_3m<273.15]=1.
        MAPDATA_DICT[config][scenario]['Permafrost3m'].append(np.copy(permafrost_3m))
        DATA_DICT[config][scenario]['PermaArea3m'].append( np.sum(permafrost_3m*AREA_1D) )

        # The Wood Products Pool:
        WP = ( DINF.variables['wood_prod_fast'][:]+        \
               DINF.variables['wood_prod_med'][:] +        \
               DINF.variables['wood_prod_slow'][:]  ) * AREA_1D
        WP = np.sum(WP)
        DATA_DICT[config][scenario]['Woody_Products'].append(WP)
        #Atmospheric CO2 
        AtmCO2_ppm = DINF.variables['co2_ppmv'][0]
        AtmCO2_kg = AtmCO2_ppm*ppm_to_kgC
        
        DATA_DICT[config][scenario]['AtmCO2_kg'].append(AtmCO2_kg)
        DATA_DICT[config][scenario]['AtmCO2_ppm'].append(AtmCO2_ppm)
        # Ocean CO2 from dtemp_o
        OCEAN_CO2 = OCEAN_DICT[config][scentag][igcm,scenyear-OCEAN_START_YEAR]
        DATA_DICT[config][scenario]['OceanCO2'].append(OCEAN_CO2)
        
        #Read in Frac Data
        FRAC = Ainf.variables['frac'][:].squeeze()
        for iTile in range(nTiles):
            MAPDATA_DICT[config][scenario][Tile_names[iTile]].append(FRAC[iTile,:])
        MAX_FRAC = np.argmax(FRAC,axis=0)
        MAPDATA_DICT[config][scenario]['Max_Frac']=MAX_FRAC
        
        FRAC = np.sum(FRAC*AREA_1D.squeeze()*1e-10,axis=1)  # m^2 to Mha
        DATA_DICT[config][scenario]['frac'].append(FRAC)
        
        TMinf.close()
        Ainf.close()
        DINF.close()
    #print('=====================================================')
            
        
###################################################################################################
# Read in the Pre Industrial Data into it's own dictionaries:
PreI_DICT={}
MapPreI_DICT={}
#print('GCM: ',gcm)
DUMP_FILE=DATA_DIR+PI_config+'/'+PI_gcm+'/'+PI_runid+'_'+PI_gcm+'_'+PI_scen+'.dump.'+str(PI_year)+'0101.0.nc'
#print(DUMP_FILE)
DINF = nc.Dataset(DUMP_FILE,'r')
Ann_File=DATA_DIR+PI_config+'/'+PI_gcm+'/'+PI_runid+'_'+PI_gcm+'_'+PI_scen+'.Annual_carbon.'+str(PI_year)+'.nc'
#print(Ann_File)
Ainf=nc.Dataset(Ann_File,'r')
TMon_File=DATA_DIR+PI_config+'/'+PI_gcm+'/'+PI_runid+'_'+PI_gcm+'_'+PI_scen+'.Monthly_h2o_t.'+str(PI_year)+'.nc'
#print(Ann_File)
TMinf=nc.Dataset(TMon_File,'r')

# Read in CV, land points only:
#CV = DINF.variables['cv'][:]
CV = Ainf.variables['cv'][:].squeeze()
# store map data
MapPreI_DICT['CV']=CV
# total:
CV = np.sum(CV*AREA_1D)
#print('CV = ',CV)
PreI_DICT['CV']=CV

# Read in Soil Carbon, and sum over pools:
#CS = np.sum(DINF.variables['cs'][:],axis=0)
# sum over layers:
#CS = np.sum( CS,axis=0 )
CS = Ainf.variables['cs_gb'][:].squeeze()
#store map data
MapPreI_DICT['CS']=CS
# global total:
CS = np.sum( CS*AREA_1D )
#print('CS = ', CS)
PreI_DICT['CS']=CS

# Permafrost stuff:
FRAC_old_cs = DINF.variables['frac_cs_old_pool_gb'][:]
CS_pool_layer = Ainf.variables['cs'][:].squeeze()
T_soil = TMinf.variables['t_soil'][:].squeeze().max(axis=0)

OldC = np.sum( np.sum(CS_pool_layer*FRAC_old_cs, axis=0), axis=0)
MapPreI_DICT['OldCS']= np.copy(OldC) 
PreI_DICT['OldCS']=  np.sum(OldC*AREA_1D) 

PermaMask = np.zeros_like(T_soil)
PermaMask[T_soil<273.15] = 1.0
OldCS_layer = np.sum(CS_pool_layer*FRAC_old_cs,axis=0)
PermaC = np.sum( OldCS_layer*PermaMask, axis=0 )
MapPreI_DICT['PermaCS']=np.copy(PermaC)
PreI_DICT['PermaCS']=np.sum(PermaC*AREA_1D)

Tmax_0m = T_soil.max(axis=0)
permafrost_0m = np.zeros_like(Tmax_0m)
permafrost_0m[Tmax_0m<273.15]=1.
MapPreI_DICT['Permafrost0m']=np.copy(permafrost_0m)
PreI_DICT['PermaArea0m']=np.sum(permafrost_0m*AREA_1D)

Tmax_1m = T_soil[z_soil>1.,:].max(axis=0)
permafrost_1m = np.zeros_like(Tmax_1m)
permafrost_1m[Tmax_1m<273.15]=1.
MapPreI_DICT['Permafrost1m']= np.copy(permafrost_1m)
PreI_DICT['PermaArea1m']= np.sum(permafrost_1m*AREA_1D)

Tmax_2m = T_soil[z_soil>2.,:].max(axis=0)
permafrost_2m = np.zeros_like(Tmax_2m)
permafrost_2m[Tmax_2m<273.15]=1.
MapPreI_DICT['Permafrost2m']= np.copy(permafrost_2m)
PreI_DICT['PermaArea2m']= np.sum(permafrost_2m*AREA_1D) 

Tmax_3m = T_soil[-1,:]
permafrost_3m = np.zeros_like(Tmax_3m)
permafrost_3m[Tmax_3m<273.15]=1.
MapPreI_DICT['Permafrost3m']=np.copy(permafrost_3m)
PreI_DICT['PermaArea3m']=np.sum(permafrost_3m*AREA_1D)

# Fill Woody Products with a dummy zero
PreI_DICT['Woody_Products']=0.0
        
#Atmospheric CO2 
AtmCO2_ppm = DINF.variables['co2_ppmv'][0]
AtmCO2_kg = AtmCO2_ppm*ppm_to_kgC
#print('AtmCO2 = ',AtmCO2_kg)
PreI_DICT['AtmCO2_kg']=AtmCO2_kg
PreI_DICT['AtmCO2_ppm']=AtmCO2_ppm
        
# Ocean CO2 from dtemp_o
OCEAN_CO2 = 0
PreI_DICT['OceanCO2']=OCEAN_CO2
#print('OCEAN_CO2 = ',OCEAN_CO2)

#FRAC = DINF.variables['frac'][:]
FRAC = Ainf.variables['frac'][:].squeeze()
for iTile in range(nTiles):
    MapPreI_DICT[Tile_names[iTile]]=FRAC[iTile,:]

MAX_FRAC = np.argmax(FRAC,axis=0)
MapPreI_DICT['Max_Frac']=MAX_FRAC
#print(MAX_FRAC.shape)
#print(MAX_FRAC)
FRAC = np.sum(FRAC*AREA_1D.squeeze()*1e-10,axis=1)  # m^2 to Mha
#print(FRAC.shape)
PreI_DICT['frac']=FRAC
#print(FRAC)
        
DINF.close()
Ainf.close()
TMinf.close()

#pdb.set_trace()

# Output Tabular data:
#####################################################################################################
# Permafrost stuff:
# From Pre-Industrial:
output_values=('Csoil','OldCS','PermaCS','delCsoil','delOldC','delPermaC')
novals=len(output_values)
outf=open(PLOT_DIR+'Permafrost_PI.txt','w')
outf.write( 'Carbon Stores: \n')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals*'%15s ' % output_values+'\n' )
outf.write( (42+novals*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = (PreI_DICT['CS']*1e-12,PreI_DICT['OldCS']*1e-12,PreI_DICT['PermaCS']*1e-12,
                 0,0,0)
    outf.write( '%25s %15s'%(config,scenario)+novals*'%15.3f '%(outdata)+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        outdata = ( np.mean(DATA_DICT[config][scenario]['CS'])*1e-12,
                    np.mean(DATA_DICT[config][scenario]['OldCS'])*1e-12,
                    np.mean(DATA_DICT[config][scenario]['PermaCS'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['CS'])-PreI_DICT['CS'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['OldCS'])-PreI_DICT['OldCS'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['PermaCS'])-PreI_DICT['PermaCS'])*1e-12,
                    )
        #print(novals,len(outdata))
        outf.write( '%25s %15s'%(config,scenario)+novals*'%15.3f '%outdata+'\n' )
    outf.write( (42+novals*16)*'_'+'\n' )
outf.write( (42+novals*16)*'#'+'\n\n\n' )
outf.write( (42+novals*16)*'#'+'\n' )

output_values2=('PermaArea0m','PermaArea1m','PermaArea2m','PermaArea3m')
novals2=len(output_values2)
outf.write( 'Permafrost Area (Mha): \n')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals2*'%15s ' % output_values2+'\n' )
outf.write( (42+novals2*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = tuple([ PreI_DICT[var]*1e-10 for var in output_values2 ])
    outf.write( '%25s %15s'%(config,scenario)+novals2*'%15.3f '%outdata+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        outdata = tuple( [np.mean(DATA_DICT[config][scenario][var])*1e-10 for var in output_values2] )
        outf.write( '%25s %15s'%(config,scenario)+novals2*'%15.3f '%outdata+'\n' )
    outf.write( (42+novals2*16)*'_'+'\n' )

outf.write( (42+novals*16)*'#'+'\n\n\n' )
outf.write( (42+novals*16)*'#'+'\n' )

centiles = [25,50,75,95,100]
output_values3=('25%','50%','75%','95%','100%')
novals3=len(centiles)
outf.write( 'Permafrost Area (Mha) at 3m: \n')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals3*'%15s ' % output_values3+'\n' )
outf.write( (42+novals3*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = tuple([PreI_DICT['PermaArea3m']*1e-10 for i in range(novals3)])
    outf.write( '%25s %15s'%(config,scenario)+novals3*'%15i '%outdata+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        mean_data = np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost3m'],axis=0)*100.
        outdata = tuple( [np.sum(AREA_1D[mean_data>=centile])*1e-10 for centile in centiles ] ) 
        outf.write( '%25s %15s'%(config,scenario)+novals3*'%15i '%outdata+'\n' )
    outf.write( (42+novals3*16)*'_'+'\n' )

outf.write( (novals*16)*'#'+'\n' )

centiles = [25,50,75,95,100]
output_values3=('25%','50%','75%','95%','100%')
novals3=len(centiles)
outf.write( 'Permafrost Area (Mha) at 1m: \n')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals3*'%15s ' % output_values3+'\n' )
outf.write( (42+novals3*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = tuple([PreI_DICT['PermaArea1m']*1e-10 for i in range(novals3)])
    outf.write( '%25s %15s'%(config,scenario)+novals3*'%15i '%outdata+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        mean_data = np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost1m'],axis=0)*100.
        outdata = tuple( [np.sum(AREA_1D[mean_data>=centile])*1e-10 for centile in centiles ] ) 
        outf.write( '%25s %15s'%(config,scenario)+novals3*'%15i '%outdata+'\n' )
    outf.write( (42+novals3*16)*'_'+'\n' )

outf.close()
#quit()
##################################################################################################
# Plot map of Permafrost Area Quartiles at 3m
#COL_MAP   = col.ListedColormap(['#dcc7b8','#c39f85','#a97773','#81536f'],'indexed')
COL_MAP   = col.ListedColormap(['#dcc7b8','#a97773','#81536f'],'indexed')
plt.cm.register_cmap(cmap=COL_MAP)
PALETTE   = COL_MAP
PREI_PAL  = col.ListedColormap(['silver','silver'])
plt.cm.register_cmap(cmap=PREI_PAL)
centiles = [25,75,95,100]
NORM      = col.BoundaryNorm(centiles, len(centiles)-1 ) #, clip=False)
parallels = np.arange(0.,90.,30.)
meridians = np.arange(-180.,180.,30)
LON_0=75
BOUNDINGLAT=50
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    
    preI_plot_data =  MapPreI_DICT['Permafrost3m'][land_index]
    preI_plot_data[land_index.mask==True]=0
    preI_plot_data =  np.ma.masked_equal(preI_plot_data,0)
    preI_plot_data.data[preI_plot_data.mask==True]=preI_plot_data.fill_value

    fig,axes=plt.subplots(nrows=1,ncols=nSCENARIOs,figsize=(6*nSCENARIOs,6) )
    fig.subplots_adjust(left=0.03,right=0.97,bottom=0.1)
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        ax=axes[iscenario]
        M = Basemap(projection='npstere',lon_0=LON_0,boundinglat=BOUNDINGLAT,ax=ax,resolution='i')
        xx,yy = M(lons_2d,lats_2d)
        mean_data = np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost3m'],axis=0)[land_index]*100.
        mean_data[land_index.mask==True]=0
        mean_data = np.ma.masked_less(mean_data,centiles[0])
        plot_data = mean_data
        M.pcolormesh(xx,yy,preI_plot_data,cmap=PREI_PAL)
        IMAGE=M.pcolormesh(xx,yy,plot_data,cmap=PALETTE,norm=NORM)
    
        M.drawcountries(linewidth =0.4)
        M.drawcoastlines(linewidth =0.4)
        M.drawparallels(parallels,labels =[0,0,0,0], linewidth = 0.5,fontsize=12)
        if iscenario==0:
            M.drawmeridians(meridians,labels =[1,0,1,1], linewidth = 0.5,fontsize=12)
        elif iscenario==nSCENARIOs-1:
            M.drawmeridians(meridians,labels =[0,1,1,1], linewidth = 0.5,fontsize=12)
        else:
            M.drawmeridians(meridians,labels =[0,0,1,1], linewidth = 0.5,fontsize=12)
    
        ax.set_title(scenario,fontsize=16)
    
    CBAR=plt.colorbar(IMAGE,ax=axes.flatten().tolist(),orientation='horizontal',fraction=0.05,pad=0.06) 
    CBAR.ax.axhline(linewidth=2,color='black')
    CBAR.ax.axhline(y=1,linewidth=2,color='black')
    for tickloc in CBAR.ax.get_xticks():  
        CBAR.ax.axvline(tickloc,linewidth=2,color='black')
    CBAR.ax.xaxis.set_tick_params(length=1,labelsize=12)
    CBAR.set_label('Percent of GCMs',fontsize=18)

    fig.savefig(PLOT_DIR+config+'_3m_Permafrost_Area_Quartiles.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+config+'_3m_Permafrost_Area_Quartiles.eps',bbox_inches='tight')
    plt.close()


##################################################################################################
# Plot map of Permafrost Area Quartiles at 1m
#COL_MAP   = col.ListedColormap(['#dcc7b8','#c39f85','#a97773','#81536f'],'indexed')
COL_MAP   = col.ListedColormap(['#dcc7b8','#a97773','#81536f'],'indexed')
plt.cm.register_cmap(cmap=COL_MAP)
PALETTE   = COL_MAP
PREI_PAL  = col.ListedColormap(['silver','silver'])
plt.cm.register_cmap(cmap=PREI_PAL)
centiles = [25,75,95,100]
NORM      = col.BoundaryNorm(centiles, len(centiles)-1 ) #, clip=False)
parallels = np.arange(0.,90.,30.)
meridians = np.arange(-180.,180.,30)
LON_0=75
BOUNDINGLAT=50
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    
    preI_plot_data =  MapPreI_DICT['Permafrost1m'][land_index]
    preI_plot_data[land_index.mask==True]=0
    preI_plot_data =  np.ma.masked_equal(preI_plot_data,0)
    preI_plot_data.data[preI_plot_data.mask==True]=preI_plot_data.fill_value

    fig,axes=plt.subplots(nrows=1,ncols=nSCENARIOs,figsize=(6*nSCENARIOs,6) )
    fig.subplots_adjust(left=0.03,right=0.97,bottom=0.1)
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        ax=axes[iscenario]
        M = Basemap(projection='npstere',lon_0=LON_0,boundinglat=BOUNDINGLAT,ax=ax,resolution='i')
        xx,yy = M(lons_2d,lats_2d)
        mean_data = np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost1m'],axis=0)[land_index]*100.
        mean_data[land_index.mask==True]=0
        mean_data = np.ma.masked_less(mean_data,centiles[0])
        plot_data = mean_data
        M.pcolormesh(xx,yy,preI_plot_data,cmap=PREI_PAL)
        IMAGE=M.pcolormesh(xx,yy,plot_data,cmap=PALETTE,norm=NORM)
    
        M.drawcountries(linewidth =0.4)
        M.drawcoastlines(linewidth =0.4)
        M.drawparallels(parallels,labels =[0,0,0,0], linewidth = 0.5,fontsize=12)
        if iscenario==0:
            M.drawmeridians(meridians,labels =[1,0,1,1], linewidth = 0.5,fontsize=12)
        elif iscenario==nSCENARIOs-1:
            M.drawmeridians(meridians,labels =[0,1,1,1], linewidth = 0.5,fontsize=12)
        else:
            M.drawmeridians(meridians,labels =[0,0,1,1], linewidth = 0.5,fontsize=12)
        ax.set_title(scenario,fontsize=16)
    
    CBAR=plt.colorbar(IMAGE,ax=axes.flatten().tolist(),orientation='horizontal',fraction=0.05,pad=0.06) 
    CBAR.ax.axhline(linewidth=2,color='black')
    CBAR.ax.axhline(y=1,linewidth=2,color='black')
    for tickloc in CBAR.ax.get_xticks():  
        CBAR.ax.axvline(tickloc,linewidth=2,color='black')
    CBAR.ax.xaxis.set_tick_params(length=1,labelsize=12)
    CBAR.set_label('Percent of GCMs',fontsize=18)

    fig.savefig(PLOT_DIR+config+'_1m_Permafrost_Area_Quartiles.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+config+'_1m_Permafrost_Area_Quartiles.eps',bbox_inches='tight')
    plt.close()
quit()

###################################################################################################
# Bar plot of Absolute Stocks:
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    # Open a nice big figure with one axes:
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[15,7])
    plotvar_list = ['CS','OldCS','PermaCS']
    plotcolour_list=['#8B4513','#8B8513','#8B4579']
    legend_names=['Soil Carbon','Defrosted Permafrost Carbon', 'Permafrost Carbon']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    #   add ocean and BECCS and Woody products to this at some point
    scenario_space = 1.5 # how much space all the bars should take up for a scenario 
    # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    iscenario=0
    scen_cen_pos = iscenario+0.5
    bar_list=[]  # Append the bar objects to list for legend
    for ibar in range(nbars):
        xpos = scen_cen_pos-(bar_width/2.)
        plotvar=plotvar_list[ibar]
        plotcolour = plotcolour_list[ibar]
        plotdata   = np.array(PreI_DICT[plotvar])*1e-12
        bar_list.append( ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)  )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+1.5
        for ibar in range(nbars):
            xpos = scen_cen_pos-(bar_width/2.)
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = np.array(DATA_DICT[config][scenario][plotvar])*1e-12
            ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    ax.set_xlim([0,nSCENARIOs+1.])
    ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    ax.set_xticks(np.arange(0.5,nSCENARIOs+1,1.))
    ax.set_xticklabels(['Pre-Industrial']+SCENARIO_names,fontsize=15)
    ax.set_ylabel('Carbon Store (GtC)',fontsize=20)
    ax.tick_params(axis='y',labelsize=10)
    ax.plot(ax.get_xlim(),[PreI_DICT['PermaCS']*1e-12,PreI_DICT['PermaCS']*1e-12],color='k')
    ax.grid(True)
    
    fig.legend(bar_list,legend_names,     
            loc='upper center',ncol=4,fontsize=20)
    fig.savefig(PLOT_DIR+config+'_Permafrost_CarbonStores.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+config+'_Permafrost_CarbonStores.eps',bbox_inches='tight')
    # Store as encasulated post script
    #plt.show()
    plt.close()


##################################################################################################
# Plot map of Permafrost Area
#COL_MAP   = col.ListedColormap(['white','#dcc7b8','#ad7c65','#7b4550','#7b4590'],'indexed')
# Don't bother with 2m
COL_MAP   = col.ListedColormap(['white','#dcc7b8','#7b4550','#7b4590'],'indexed')
plt.cm.register_cmap(cmap=COL_MAP)
PALETTE   = COL_MAP
NORM      = col.BoundaryNorm([0,0.8,1.6,2.4], 5, clip=False)
parallels = np.arange(0.,90.,15)
meridians = np.arange(-180.,180.,30)
LON_0=0
BOUNDINGLAT=30
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    
    #fig,axes=plt.subplots(nrows=nSCENARIOs+1,ncols=1,figsize=(24,(nSCENARIOs+1)*6) )
    fig,axes=plt.subplots(ncols=nSCENARIOs+1,nrows=1,figsize=((nSCENARIOs+1)*6,8))
    fig.subplots_adjust(bottom=0.1)
    scenario='Pre-Industrial'
    ax=axes[0]
    #plot_data =  MapPreI_DICT['Permafrost3m']+MapPreI_DICT['Permafrost2m'] \
    plot_data =  MapPreI_DICT['Permafrost3m'] \
               + MapPreI_DICT['Permafrost1m']+MapPreI_DICT['Permafrost0m']
    plot_data = np.ma.masked_array(plot_data[land_index],mask=land_index.mask,
                                    fill_value=-1)
    plot_data.data[plot_data.mask==True]=plot_data.fill_value
   
    M = Basemap(projection='npstere',lon_0=LON_0,boundinglat=BOUNDINGLAT,ax=ax,resolution='i')
    #M = Basemap(projection='cyl',
    #             llcrnrlat = 20 ,urcrnrlat=90,llcrnrlon =-180 ,urcrnrlon=180, 
    #             ax=ax)
    xx,yy = M(lons_2d,lats_2d)
    IMAGE=M.pcolormesh(xx,yy,plot_data,cmap=PALETTE) #,norm=NORM)
    M.drawcountries(linewidth =0.4)
    M.drawcoastlines(linewidth =0.4)
    #M.drawparallels(parallels, labels =[1,1,0,0], linewidth = 0.5,fontsize=20)
    #M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5,fontsize=20)
    M.drawparallels(parallels, labels =[0,0,0,0], linewidth = 0.5,fontsize=20)
    M.drawmeridians(meridians, labels =[0,0,1,1], linewidth = 0.5,fontsize=20)
    #ax.set_title(scenario,fontsize=35)


    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        ax=axes[iscenario+1]
        plot_data =  np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost3m'],axis=0) \
                    +np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost1m'],axis=0) \
                    +np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC'][scenario]['Permafrost0m'],axis=0) 
                    #+np.mean(MAPDATA_DICT['BASELINE_CONFIG_NoLULUC']['presentday']['Permafrost2m'],axis=0) \

        plot_data = np.ma.masked_array(plot_data[land_index],mask=land_index.mask,
                                        fill_value=-1)
        plot_data.data[plot_data.mask==True]=plot_data.fill_value
    
        M = Basemap(projection='npstere',lon_0=LON_0,boundinglat=BOUNDINGLAT,ax=ax,resolution='i')
        #M = Basemap(projection='cyl',
        #               llcrnrlat=20 ,urcrnrlat=90,llcrnrlon =-180 ,urcrnrlon=180, 
        #               ax=ax)
        M.pcolormesh(xx,yy,plot_data,cmap=PALETTE) #,norm=NORM)
        M.drawcountries(linewidth =0.4)
        M.drawcoastlines(linewidth =0.4)
        M.drawparallels(parallels, labels =[1,1,0,0], linewidth = 0.5,fontsize=20)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5,fontsize=20)
        #ax.set_title(SCENARIO_names[iscenario],fontsize=35)

    CBAR=plt.colorbar(IMAGE,ax=axes.flatten().tolist(),orientation='horizontal',fraction=0.05,pad=0.03,
                        ticks=np.array([0.4,1.1,1.9,2.6])) 
                        #ticks=np.array([0.3,0.9,1.5,2.1])) 
    CBAR.ax.axhline(linewidth=2,color='black')
    CBAR.ax.axhline(y=1,linewidth=2,color='black')
    #pdb.set_trace()
    nticks=len(CBAR.ax.get_xticks())
    print(nticks)
    for itick in range(nticks):
        print((1./nticks)+(float(itick)/nticks))
        CBAR.ax.axvline(x=(1./nticks)+(float(itick)/nticks) ,linewidth=2,color='black')
    #CBAR.ax.axvline(x=0.5,linewidth=2,color='black')
    #CBAR.ax.axvline(x=0.75,linewidth=2,color='black')
    #CBAR.ax.set_xticks(np.array([0.8,1.6,2.4,3.2,4.0]))
    #CBAR.ax.set_xticklabels(['None','3m','2m','1m','0m'],fontsize=20)
    CBAR.ax.set_xticklabels(['None','3m','1m','0m'],fontsize=20)
    CBAR.ax.xaxis.set_tick_params(length=0)
    
    fig.savefig(PLOT_DIR+config+'_Permafrost_Area.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+config+'_Permafrost_Area.eps',bbox_inches='tight')
    plt.close()
#quit()

