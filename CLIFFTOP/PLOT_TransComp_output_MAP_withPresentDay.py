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
CONFIG2 = optional_argparse('-config2', 'METHANE_FEEDBACK_NoLULUC')
RUNID2  = optional_argparse('-runid2',  'CH4_FB')
CONFIGS = [CONFIG1,CONFIG2]
RUNIDS  = [RUNID1,RUNID2]
nCONFIGS = len(CONFIGS)
print('Configs: ',CONFIGS)

# Directories containing JULES output and plot output directories:i
DATA_DIR = optional_argparse('-datadir','/prj/CLIFFTOP/ECP_output/')
print('DATA_DIR: '+DATA_DIR)
#DATA_DIR = '/group_workspaces/jasmin2/clifftop/CLIFFTOP/ECP_output/EQUILIBRIUM_OUTPUT/'+CONFIG+'/'

PLOT_TAG = optional_argparse('-plottag', 'BL_CH4FB_NoLULUC')
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
lats_2d = grinf.variables['latitude'][:]
lons_2d = grinf.variables['longitude'][:]
#Area_2d = grinf.variables['Area'][:]       # I don't actually use this but it's here
land_index = grinf.variables['land_index'][:]
grinf.close()

# 1Dimension grid cell area data for calculating totals etc.
AREA_file=ANCILS_DIR+'Area_in_iris_format.nc'
Ainf=nc.Dataset(AREA_file,'r')
AREA_1D = Ainf.variables['area'][:].squeeze()
#lats_1D = Ainf.variables['latitude'][:].squeeze()
#lons_1D = Ainf.variables['longitude'][:].squeeze()
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
    print(outdata)
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
outf.write( 'Permafrost Area: \n')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals2*'%15s ' % output_values2+'\n' )
outf.write( (42+novals2*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = tuple([ PreI_DICT[var]*1e-10 for var in output_values2 ])
    print(outdata)
    outf.write( '%25s %15s'%(config,scenario)+novals2*'%15.3f '%outdata+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        outdata = tuple( [np.mean(DATA_DICT[config][scenario][var])*1e-10 for var in output_values2] )
        print(novals2,len(outdata))
        outf.write( '%25s %15s'%(config,scenario)+novals2*'%15.3f '%outdata+'\n' )
    outf.write( (42+novals2*16)*'_'+'\n' )

outf.close()


#####################################################################################################
# C-Stocks
# From Pre-Industrial:
output_values=('Cveg','Csoil','Catm','delCveg','delCsoil','delCatm','delCocean','delC')
novals=len(output_values)
outf=open(PLOT_DIR+'Carbon_Stores_PI.txt','w')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals*'%15s ' % output_values+'\n' )
outf.write( (42+novals*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    scenario = 'Pre-Industrial'
    outdata = (PreI_DICT['CV']*1e-12,PreI_DICT['CS']*1e-12,PreI_DICT['AtmCO2_kg']*1e-12,
                 0,0,0,0,0)
    outf.write( '%25s %15s'%(config,scenario)+novals*'%15.3f '%(outdata)+'\n' )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        outdata = [ np.mean(DATA_DICT[config][scenario]['CV'])*1e-12,
                    np.mean(DATA_DICT[config][scenario]['CS'])*1e-12,
                    np.mean(DATA_DICT[config][scenario]['AtmCO2_kg'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['CV'])-PreI_DICT['CV'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['CS'])-PreI_DICT['CS'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['AtmCO2_kg'])-PreI_DICT['AtmCO2_kg'])*1e-12,
                    (np.mean(DATA_DICT[config][scenario]['OceanCO2'])-PreI_DICT['OceanCO2'])*1e-12,
                    ]
        outdata = tuple( outdata + [outdata[3]+outdata[4]+outdata[5]+outdata[6]] )
        outf.write( '%25s %15s'%(config,scenario)+novals*'%15.3f '%outdata+'\n' )
    outf.write( (42+novals*16)*'_'+'\n' )
outf.close()

# From Present Day:
outf=open(PLOT_DIR+'Carbon_Stores_PD.txt','w')
outf.write( '%25s %15s'%('Configuration','Scenario')+novals*'%15s ' % output_values+'\n' )
outf.write( (42+novals*16)*'_'+'\n' )
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        outdata = [ 
           np.mean(DATA_DICT[config][scenario]['CV'])*1e-12,
           np.mean(DATA_DICT[config][scenario]['CS'])*1e-12,
           np.mean(DATA_DICT[config][scenario]['AtmCO2_kg'])*1e-12,
           (np.mean(DATA_DICT[config][scenario]['CV'])-np.mean(DATA_DICT[config]['presentday']['CV']))*1e-12,
           (np.mean(DATA_DICT[config][scenario]['CS'])-np.mean(DATA_DICT[config]['presentday']['CS']))*1e-12,
           (np.mean(DATA_DICT[config][scenario]['AtmCO2_kg'])-np.mean(DATA_DICT[config]['presentday']['AtmCO2_kg']))*1e-12,
           (np.mean(DATA_DICT[config][scenario]['OceanCO2'])-np.mean(DATA_DICT[config]['presentday']['OceanCO2']))*1e-12,
                    ]
        outdata = tuple( outdata + [outdata[3]+outdata[4]+outdata[5]+outdata[6]] )
        print(outdata) 
        outf.write( '%25s %15s'%(config,scenario)+novals*'%15.3f '%outdata+'\n' )
    outf.write( (42+novals*16)*'_'+'\n' )
outf.close()



###################################################################################################
# Bar plot of Absolute Stocks:
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    # Open a nice big figure with one axes:
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
    plotvar_list = ['CV','CS','AtmCO2_kg']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki']
    legend_names=['Vegetation','Soil','Atmosphere']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    #   add ocean and BECCS and Woody products to this at some point
    scenario_space = 0.8 # how much space all the bars should take up for a scenario 
    # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    #Plot Pre-Industrial First
    iscenario=0
    scen_cen_pos = iscenario+0.5
    bar_list=[]  # Append the bar objects to list for legend
    for ibar in range(nbars):
        xpos = scen_cen_pos+bar_positions[ibar]
        plotvar=plotvar_list[ibar]
        plotcolour = plotcolour_list[ibar]
        plotdata   = np.array(PreI_DICT[plotvar])*1e-12
        bar_list.append( ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)  )
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+1.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = np.array(DATA_DICT[config][scenario][plotvar])*1e-12
            ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            print(plotvar,plotdata.shape,nGCMs)
            ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    ax.set_xlim([0,nSCENARIOs+1.])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(np.arange(0.5,nSCENARIOs+1,1.))
    ax.set_xticklabels(['Pre-Industrial']+SCENARIO_names,fontsize=25)
    ax.set_ylabel('Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
    #Plot present day CO2 in red:
    present_CO2=400/GtC_to_ppm
    ax2.plot(ax2.get_xlim(),[present_CO2,present_CO2],c='r',lw=3,ls='--')
        
    fig.legend(bar_list,legend_names,     
            loc='upper center',ncol=4,fontsize=30)
    fig.savefig(PLOT_DIR+config+'_Equilibrium_CarbonStores.png',bbox_inches='tight')  # Store as png
    fig.savefig(PLOT_DIR+config+'_Equilibrium_CarbonStores.png',bbox_inches='tight')  # Store as encasulated post script
    #plt.show()
    plt.close()

###################################################################################################
# Bar plot of Absolute Stocks:
fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
plotcolour_lists=[ ['darkolivegreen','saddlebrown','khaki'],['green','red','yellow'] ]
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    # Open a nice big figure with one axes:
    plotvar_list = ['CV','CS','AtmCO2_kg']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki']  #plotcolour_lists[iconfig]
    legend_names=['Vegetation','Soil','Atmosphere']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    #   add ocean and BECCS and Woody products to this at some point
    scenario_space = 0.8 # how much space all the bars should take up for a scenario 
    # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    #Plot Pre-Industrial First
    iscenario=0
    scen_cen_pos = iscenario+0.5
    bar_list=[]  # Append the bar objects to list for legend
    for ibar in range(nbars):
        xpos = scen_cen_pos+bar_positions[ibar]
        plotvar=plotvar_list[ibar]
        plotcolour = plotcolour_list[ibar]
        plotdata   = np.array(PreI_DICT[plotvar])*1e-12
        bar_list.append( ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)  )
        if iconfig==0:
            ax.bar(xpos,np.mean(plotdata),color='white',alpha=0.5,width=bar_width)
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+1.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = np.array(DATA_DICT[config][scenario][plotvar])*1e-12
            ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            if iconfig==1:
                ax.bar(xpos,np.mean(plotdata),color='white',alpha=0.5,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            print(plotvar,plotdata.shape,nGCMs)
            ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    ax.set_xlim([0,nSCENARIOs+1.])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(np.arange(0.5,nSCENARIOs+1,1.))
    ax.set_xticklabels(['Pre-Industrial']+SCENARIO_names,fontsize=25)
    ax.set_ylabel('Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
        
fig.legend(bar_list,legend_names,     
           loc='upper center',ncol=4,fontsize=30)
fig.savefig(PLOT_DIR+config+'_Equilibrium_CarbonStores.png',bbox_inches='tight')  # Store as png
fig.savefig(PLOT_DIR+config+'_Equilibrium_CarbonStores.png',bbox_inches='tight')  # Store as encasulated post script
#plt.show()
plt.close()

############################################################
# plot Delta Cstores
# Open a nice big figure with one axes:
fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    plotvar_list = ['CV','CS','AtmCO2_kg','OceanCO2']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki','cornflowerblue']
    legend_names=['Vegetation','Soil','Atmosphere','Ocean']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    scenario_space = 0.8 # how much space all the bars should take up for a scenario
                            # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    bar_list=[]
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+0.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = (np.array(DATA_DICT[config][scenario][plotvar])-PreI_DICT[plotvar])*1e-12
            if iscenario == 0:
                bar_list.append(ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width))
            else:
                ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            if iconfig==1:
                ax.bar(xpos,np.mean(plotdata),color='white',alpha=0.5,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            ax.plot([xpos-0.01+(iconfig/50.)+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    
    ax.set_xlim([0,nSCENARIOs])
    # plot zero line
    ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.tick_params(labelright=True)
    ax.set_xticks(np.arange(0.5,nSCENARIOs,1.))
    ax.set_xticklabels(SCENARIO_names,fontsize=25)
    ax.set_ylabel('Change in Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('$\Delta$CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
fig.legend(bar_list,legend_names,     
          loc='upper center',ncol=4,fontsize=30)
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores_2data.png',bbox_inches='tight')
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores_2data.png',bbox_inches='tight')
#plt.show()
plt.close()
############################################################
# plot Delta Cstores
# Open a nice big figure with one axes:
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
    plotvar_list = ['CV','CS','AtmCO2_kg','OceanCO2']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki','cornflowerblue']
    legend_names=['Vegetation','Soil','Atmosphere','Ocean']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    scenario_space = 0.8 # how much space all the bars should take up for a scenario
                            # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    bar_list=[]
    for iscenario in range(nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario+0.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = (np.array(DATA_DICT[config][scenario][plotvar])-PreI_DICT[plotvar])*1e-12
            if iscenario == 0:
                bar_list.append(ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width))
            else:
                ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            ax.plot([xpos-0.01+(iconfig/50.)+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    
    ax.set_xlim([0,nSCENARIOs])
    # plot zero line
    ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.tick_params(labelright=True)
    ax.set_xticks(np.arange(0.5,nSCENARIOs,1.))
    ax.set_xticklabels(SCENARIO_names,fontsize=25)
    ax.set_ylabel('Change in Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('$\Delta$CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
    fig.legend(bar_list,legend_names,     
              loc='upper center',ncol=4,fontsize=30)
    fig.savefig(PLOT_DIR+config+'_Equilibrium_DeltaCarbonStores.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+config+'_Equilibrium_DeltaCarbonStores.png',bbox_inches='tight')
    #plt.show()
    plt.close()
#quit()
############################################################
# plot Delta Cstores
# Open a nice big figure with one axes:
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
    plotvar_list = ['CV','CS','AtmCO2_kg','OceanCO2']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki','cornflowerblue']
    legend_names=['Vegetation','Soil','Atmosphere','Ocean']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    scenario_space = 0.8 # how much space all the bars should take up for a scenario
                            # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    bar_list=[]
    for iscenario in range(1,nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario-0.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = (np.array(DATA_DICT[config][scenario][plotvar])-
                        np.array(DATA_DICT[config]['presentday'][plotvar]))*1e-12
            if iscenario == 1:
                bar_list.append(ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width))
            else:
                ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            ax.plot([xpos-0.01+(iconfig/50.)+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    
    ax.set_xlim([0,nSCENARIOs-1])
    # plot zero line
    ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
    ax.spines['top'].set_visible(False)
    #ax.tick_params(labelright=True)
    ax.set_xticks(np.arange(0.5,nSCENARIOs-1,1.))
    ax.set_xticklabels(SCENARIO_names[1:],fontsize=25)
    ax.set_ylabel('Change in Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('$\Delta$CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
    fig.legend(bar_list,legend_names,     
              loc='upper center',ncol=4,fontsize=30)
    fig.savefig(PLOT_DIR+config+'Equilibrium_DeltaCarbonStores_frompres.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+config+'Equilibrium_DeltaCarbonStores_frompres.eps',bbox_inches='tight')
    #plt.show()
    plt.close()

############################################################
# plot Delta Cstores
# Open a nice big figure with one axes:
fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    plotvar_list = ['CV','CS','AtmCO2_kg','OceanCO2']
    plotcolour_list=['darkolivegreen','saddlebrown','khaki','cornflowerblue']
    legend_names=['Vegetation','Soil','Atmosphere','Ocean']
    nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
    scenario_space = 0.8 # how much space all the bars should take up for a scenario
                            # maximum is 1 where the bars will be touching the next scenario bars
    bar_width = scenario_space/nbars
    # position of the bars relative to the scenario central position
    #  this is the left start point of the bar:
    bar_positions = [ -(scenario_space/2.)+(i*bar_width) for i in range(nbars) ]
    bar_list=[]
    for iscenario in range(1,nSCENARIOs):
        scenario=SCENARIOs[iscenario]
        scen_cen_pos = iscenario-0.5
        for ibar in range(nbars):
            xpos = scen_cen_pos+bar_positions[ibar]
            plotvar=plotvar_list[ibar]
            plotcolour = plotcolour_list[ibar]
            plotdata = (np.array(DATA_DICT[config][scenario][plotvar])-
                        np.array(DATA_DICT[config]['presentday'][plotvar]))*1e-12
            if iscenario == 1:
                bar_list.append(ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width))
            else:
                ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
            if iconfig==1:
                ax.bar(xpos,np.mean(plotdata),color='white',alpha=0.5,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
            #  Plot GCM spread as points:
            ax.plot([xpos-0.01+(iconfig/50.)+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 
    
    ax.set_xlim([0,nSCENARIOs-1])
    # plot zero line
    ax.plot(ax.get_xlim(),[0,0],c='k',lw=2)
    ax.spines['top'].set_visible(False)
    #ax.tick_params(labelright=True)
    ax.set_xticks(np.arange(0.5,nSCENARIOs-1,1.))
    ax.set_xticklabels(SCENARIO_names[1:],fontsize=25)
    ax.set_ylabel('Change in Carbon Store (GtC)',fontsize=30)
    ax.tick_params(axis='y',labelsize=20)
    ax.grid(True)
    # right hand labels in ppm
    ax2=ax.twinx()
    ax2.set_ylabel('$\Delta$CO$_2$ (ppm)',fontsize=30)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks=ax.get_yticks()
    ax2.set_yticklabels([label*GtC_to_ppm for label in ax.get_yticks()])
    ax2.tick_params(axis='y',labelsize=20)
fig.legend(bar_list,legend_names,     
          loc='upper center',ncol=4,fontsize=30)
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores_frompres.png',bbox_inches='tight')
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores_frompres.png',bbox_inches='tight')
#plt.show()
plt.close()

quit()

###################################################################################################################
#  Plot Bar chart of the Cover Fractions

# Construct a dictionary of the tile cover from the frac
Tile_Dict = { config:{ scenario:{tile:[] for tile in Tile_names} for scenario in SCENARIOs } for config in CONFIGS}
for config in CONFIGS:
    for scenario in SCENARIOs:
        for iTile in range(nTiles):
            for i_gcm in range(nGCMs):
                Tile_Dict[config][scenario][Tile_names[iTile]].append(DATA_DICT[config][scenario]['frac'][i_gcm][iTile])
            Tile_Dict[config][scenario][Tile_names[iTile]] = np.array(Tile_Dict[config][scenario][Tile_names[iTile]])
#############################################
# First Absoulute Data
# Nice big figure with axis for each scenario and the preindustrial on different rows:
for iconfig in range(nCONFIGS):
    config=CONFIGS[iconfig]
    fig,axes=plt.subplots(ncols=1,nrows=1+nSCENARIOs,figsize=[25,4*(nSCENARIOs+1)+1])
    
    # Set Y limits so constant for all plots:
    YLIMS = [0,5000]
    
    # only plot the PFTs and Bare soil. Urban, Ice and Lake do not change so are boring
    plotTiles_locs = list(range(13))+[15]
    nPLOTtiles=len(plotTiles_locs)
    
    # set bar width to one so they touch.
    bar_width=1.
    
    # plot Pre industrial on the top row:
    ax=axes[0]
    for iTilea in range(nPLOTtiles):  # Tiles):
        xpos=iTilea+0.5
        iTile=plotTiles_locs[iTilea]
        ax.bar(xpos,PreI_DICT['frac'][iTile],color=Tile_colours[iTile],width=bar_width,label=Tile_names[iTile])
    ax.set_ylabel('Area Mha',fontsize=20)
    ax.set_xticklabels(['' for i in range(nPLOTtiles)])
    ax.set_xlim([0,nPLOTtiles+1])
    ax.set_ylim(YLIMS)
    ax.set_title('Pre-Industrial',fontsize=30)
    ax.grid(True)
    
    for iscenario in range(nSCENARIOs):
        ax=axes[iscenario+1]
        scenario=SCENARIOs[iscenario]
        for iTilea in range(nPLOTtiles):  # Tiles):
            xpos=iTilea+0.5
            iTile=plotTiles_locs[iTilea]
            tile=Tile_names[iTile]
            tilecolour=Tile_colours[iTile]
            plotdata=Tile_Dict[config][scenario][tile]
            ax.bar(xpos,np.mean(plotdata),color=tilecolour,width=bar_width)
            #  Plot GCM as line:
            #ax.plot([xpos+bar_width/2.,xpos+bar_width/2.],[np.min(plotdata),np.max(plotdata],lw=3,c='k')
            #  Plot GCM spread as points:
            ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.')
        
        ax.set_ylabel('Area Mha',fontsize=20)
        ax.set_xticklabels(['' for i in range(nPLOTtiles)])
        ax.set_xlim([0,nPLOTtiles+1])
        ax.set_ylim(YLIMS)
        ax.set_title(SCENARIO_names[iscenario],fontsize=30)
        ax.grid(True)
    
    handles,labels = axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=int(np.ceil(nPLOTtiles/2.)),loc=8,fontsize=18)
    
    fig.savefig(PLOT_DIR+'Equilibrium_CoverFractions.png',bbox_inches='tight')
    fig.savefig(PLOT_DIR+'Equilibrium_CoverFractions.png',bbox_inches='tight')
    plt.show()
    #plt.close()

quit()

#######################################################################
# Now plot the delta fractional cover
# Nice big figure with axis for each scenario and the preindustrial on different rows:
fig,axes=plt.subplots(ncols=1,nrows=1+nSCENARIOs,figsize=[25,4*(nSCENARIOs+1)+1])

# Set Y limits so constant for all plots:
YLIMS = [0,5000]
YLIMS_delta = [-1500,1500]

# only plot the PFTs and Bare soil. Urban, Ice and Lake do not change so are boring
plotTiles_locs = list(range(13))+[15]
nPLOTtiles=len(plotTiles_locs)

# set bar width to one so they touch.
bar_width=1.

# plot Pre industrial on the top row:
ax=axes[0]
for iTilea in range(nPLOTtiles):  # Tiles):
    xpos=iTilea+0.5
    iTile=plotTiles_locs[iTilea]
    ax.bar(xpos,PreI_DICT['frac'][iTile],color=Tile_colours[iTile],width=bar_width,label=Tile_names[iTile])
ax.set_ylabel('Area Mha',fontsize=20)
ax.set_xticklabels(['' for i in range(nPLOTtiles)])
ax.set_xlim([0,nPLOTtiles+1])
ax.set_ylim(YLIMS)
ax.set_title('Pre-Industrial',fontsize=30)
ax.grid(True)

for iscenario in range(nSCENARIOs):
    ax=axes[iscenario+1]
    scenario=SCENARIOs[iscenario]
    for iTilea in range(nPLOTtiles):  # Tiles):
        xpos=iTilea+0.5
        iTile=plotTiles_locs[iTilea]
        tile=Tile_names[iTile]
        tilecolour=Tile_colours[iTile]
        plotdata=Tile_Dict[scenario][tile]-PreI_DICT['frac'][iTile]
        ax.bar(xpos,np.mean(plotdata),color=tilecolour,width=bar_width)
        #  Plot GCM as line:
        #ax.plot([xpos+bar_width/2.,xpos+bar_width/2.],[np.min(plotdata),np.max(plotdata],lw=3,c='k')
        #  Plot GCM spread as points:
        ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.')
               
    ax.set_ylabel('$\Delta$Area km$^2$',fontsize=20)
    #ax.set_xticklabels(['' for i in range(nTiles)])
    ax.set_xticklabels(['' for i in range(nPLOTtiles)])
    ax.set_xlim([0,nPLOTtiles+1])
    ax.set_ylim(YLIMS_delta)
    ax.set_title(SCENARIO_names[iscenario],fontsize=30)
    ax.grid(True)

handles,labels = axes[0].get_legend_handles_labels()
fig.legend(handles,labels,ncol=int(np.ceil(nPLOTtiles/2.)),loc=8,fontsize=18)
    
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCoverFractions.png',bbox_inches='tight')
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCoverFractions.png',bbox_inches='tight')
#plt.show()
plt.close()
#quit()


# Plot differences between scenarios, THIS IS A START BUT NEEDS TIDYING
#fig,axes=plt.subplots(ncols=1,nrows=1,figsize=[20,10])
#
#plotTiles_locs = list(range(13))+[15]
#nPLOTtiles=len(plotTiles_locs)
#
#bar_width=1.
#OceanCO2 = np.array(DATA_DICT[scenario]['OceanCO2'])*1e-12
# Land position, iSCENARIO-0.2
#ax=axes
#for iTilea in range(nPLOTtiles):  # Tiles):
#    xpos=iTilea+0.5
#    iTile=plotTiles_locs[iTilea]
#    tile=Tile_names[iTile]
#    DIFF_DATA=Tile_Dict['1p5equi'][tile]-Tile_Dict['2equi'][tile]
#    ax.bar(xpos,np.mean(DIFF_DATA),color=Tile_colours[iTile],width=bar_width,label=Tile_names[iTile])
#    ax.plot([xpos+0.5,xpos+0.5],[np.min(DIFF_DATA),np.max(DIFF_DATA)],lw=3,c='k')
#ax.set_ylabel('Area km$^2$',fontsize=20)
#ax.set_xticklabels(['' for i in range(nPLOTtiles)])
#ax.set_xlim([0,nPLOTtiles+1])
#ax.set_yticklabels([str('%4i'%label) for label in ax.get_yticks()])
#ax.set_title('1.5 K - 2 K',fontsize=30)
#ax.grid(True)
#
#handles,labels = axes.get_legend_handles_labels()
#fig.legend(handles,labels,ncol=int(np.ceil(nPLOTtiles/2.)),loc=8,fontsize=18)
#plt.show()
#plt.close()


##############################################################################################
# PLOT MAP SECTION
##############################################################################################
# Map Variables:
VARIABLES  = ['CV','CS' ]
nVARIABLEs = len(VARIABLES)
LONG_NAMES = ['Vegetation Carbon','Soil Carbon' ]
MAP_UNITS  = ['kgC m$^{-2}$','kgC m$^{-2}$']

ABS_RANGES = [ (0,15), (0,150) ]
DIFF_RANGES = [ (-5,5),(-5,5) ]
DEV_RANGES  = [ (0,2),(0,2)   ]
                
# plot maps of CV, and the change in CV:
#variable='CV'
#long_name='Vegetation Carbon'
#map_unit='kgC m$^{-2}$'
#Abs_range=[0,15]
#Abs_Colours=['#f9f9ea','#fff68f','#fff043','#a8da61','#4fbd44']
#Diff_range=[-5,5]
#Diff_Colours=['#03487b','#eeeeee','#f23452']
#Dev_range=[0,100]
#Dev_range=[0,2]
#Dev_Colours=['#f9f9ea','#d98d0e']

#Abs_Colours=['#eeeeee','#fff68f','#ffa500','#5f1205']
#Abs_Colours=['#f9f9ea','#fff68f','#ffa500','#5f1205']
#Diff_range=[-5,5]
#Diff_Colours=['#03487b','#eeeeee','#f23452']
#Dev_range=[0,100]
#Dev_Colours=['#eeeeee','#f23452']
#Dev_range=[0,2]
#Dev_Colours=['#f9f9ea','#d98d0e']

#  My plot map routine can interploate a colour bar from a list of colours
#  Alternatively you can use a standard matplotlib one, see subroutine info ( help(PTs.plot_map) ),
#    or load a brewer one and get it into the right format, I can show you how if you really want to.
# Different colour bars for absolute data:
ABS_COLOURS = [ ('#f9f9ea','#fff68f','#fff043','#a8da61','#4fbd44') ,
                ('#f9f9ea','#fff68f','#ffa500','#5f1205'), 
                ]
# Same colour bars for difference and deviation bars
Diff_Colours=['#03487b','#eeeeee','#f23452']
Dev_Colours=['#f9f9ea','#d98d0e']

for ivar in range(nVARIABLEs):
    variable  = VARIABLES[ivar]
    long_name = LONG_NAMES[ivar]
    map_unit  = MAP_UNITS[ivar]
    Abs_range = ABS_RANGES[ivar]
    Diff_range= DIFF_RANGES[ivar]
    Dev_range = DEV_RANGES[ivar]
    Abs_Colours = ABS_COLOURS[ivar]

    PreI_map_Data=MapPreI_DICT[variable]

    # This is my way of converting the data to 2D for my plotting routine.
    #  I will one day update the routine to deal with 1D data but I haven't done yet.
    #  The mask ensures that spaces are left where there is no data
    PreI_plotdata=np.ma.masked_array(PreI_map_Data[land_index],mask=land_index.mask)
    for iscenario in range(nSCENARIOs):
        scenario  = SCENARIOs[iscenario]
        scen_name = SCENARIO_names[iscenario]
        print(scen_name)
    
        # Create figure with 4 axes:
        FIG,AXES=plt.subplots(ncols=2,nrows=2,figsize=[20,10])
        # Plot Preindustrial in top left (AXES[0,0])
        PTs.plot_map(PreI_plotdata,lons_2d,lats_2d,
                     DATA_RANGE=Abs_range,
                     COLOURS=Abs_Colours,INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,
                     RESOLUTION='c', MAP_TYPE='Contour',SET_OVER='r',
                     PLOT_TITLE='Pre-Industrial',AXIS=AXES[0,0],
                     CBAR_LABEL=map_unit,FONTSIZES=[12,12,14,18]
                    )
    
        # Plot mean of Scenario absolute data in top right (AXES[0,1])
        SCEN_mean_map_DATA=np.mean(np.array(MAPDATA_DICT[scenario][variable]),axis=0)
        # onto 2D grid:
        SCEN_mean_plotdata=np.ma.masked_array(SCEN_mean_map_DATA[land_index],mask=land_index.mask)
        PTs.plot_map(SCEN_mean_plotdata,lons_2d,lats_2d,
                DATA_RANGE=Abs_range,
                COLOURS=Abs_Colours,INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,
                RESOLUTION='c', MAP_TYPE='Contour',
                PLOT_TITLE=scen_name+' (Mean)',AXIS=AXES[0,1],
                CBAR_LABEL=map_unit,FONTSIZES=[12,12,14,18]
                )
        
        # Plot Mean Difference in bottom left (AXES[1,0])
        SCEN_DIFF_DATA = np.array(MAPDATA_DICT[scenario][variable])-PreI_map_Data
        SCEN_mean_diff_map = np.mean(SCEN_DIFF_DATA,axis=0)
        # onto 2D grid:
        SCEN_diff_plotdata=np.ma.masked_array(SCEN_mean_diff_map[land_index],mask=land_index.mask)
        PTs.plot_map(SCEN_diff_plotdata,lons_2d,lats_2d,
                DATA_RANGE=Diff_range,
                COLOURS=Diff_Colours,INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,
                RESOLUTION='c', MAP_TYPE='Contour',
                PLOT_TITLE='Pre-Industrial - '+scen_name+' (Mean)',AXIS=AXES[1,0],
                CBAR_LABEL=map_unit,FONTSIZES=[12,12,14,18]
                )
        
        # Plot Standard deviation of the GCMs in bottom right (AXES[1,0])
        SCEN_std_map = np.std(MAPDATA_DICT[scenario][variable],axis=0)
        # onto 2D grid:
        SCEN_std_plotdata=np.ma.masked_array(SCEN_std_map[land_index],mask=land_index.mask)
        PTs.plot_map(SCEN_std_plotdata,lons_2d,lats_2d,
                DATA_RANGE=Dev_range,
                COLOURS=Dev_Colours,INTERPOLATE_COLOURS=True,NLEVELS=250,NTICKS=11,
                RESOLUTION='c', MAP_TYPE='Contour',
                PLOT_TITLE=scen_name+'  (Standard Deviation)',AXIS=AXES[1,1],
                #CBAR_LABEL='\% of Mean',FONTSIZES=[12,12,14,18],
                CBAR_LABEL=map_unit,FONTSIZES=[12,12,14,18],
                )
        
        FIG.suptitle(scen_name+' Equilibrium '+long_name,fontsize=30)
        
        FIG.savefig(PLOT_DIR+'Map_'+variable+'_'+scenario+'.png',bbox_inches='tight')
        FIG.savefig(PLOT_DIR+'Map_'+variable+'_'+scenario+'.png',bbox_inches='tight')
        plt.close()




##########################################################################################################
# Plot Max Frac data:

# Nice big plot with an axis per row for each scenario + pre industrial
FIG,AXES=plt.subplots(ncols=1,nrows=nSCENARIOs+1,figsize=[15,(5*nSCENARIOs)+1])

#Plot Preindustiral in top row:
plotdata=np.ma.masked_array(MapPreI_DICT['Max_Frac'][land_index],mask=land_index.mask)
PTs.plot_map(plotdata,lons_2d,lats_2d,
             COLOURS=Tile_colours,DATA_RANGE=[0,nTiles+0.1],
             TickLEVELS=list(np.arange(0.5,nTiles+1,1)),TickLABELS=Tile_names,
             CBAR_TICK_LENGTH=0,CLEVELS=range(nTiles+1),CBAR_ORIENTATION='vertical',
             PLOT_TITLE='Pre-Industrial',FONTSIZES=[12,12,14,18],
             RESOLUTION='c',
             AXIS=AXES[0],
            )

for iscenario in range(nSCENARIOs):
    scenario=SCENARIOs[iscenario]
    scen_name=SCENARIO_names[iscenario]
    plotdata=np.ma.masked_array(MAPDATA_DICT[scenario]['Max_Frac'][land_index],mask=land_index.mask)
    PTs.plot_map(plotdata,lons_2d,lats_2d,
            COLOURS=Tile_colours,DATA_RANGE=[0,nTiles+0.1],
            TickLEVELS=list(np.arange(0.5,nTiles+1,1)),TickLABELS=Tile_names,
            CBAR_TICK_LENGTH=0,CLEVELS=range(nTiles+1),CBAR_ORIENTATION='vertical',
            PLOT_TITLE=scen_name,FONTSIZES=[12,12,14,18],
            RESOLUTION='c',
            AXIS=AXES[iscenario+1],
            )


FIG.suptitle('Maximum Gridcell Land-Cover',fontsize=30)
FIG.savefig(PLOT_DIR+'Map_MaxLandCover.png',bbox_inches='tight')   
FIG.savefig(PLOT_DIR+'Map_MaxLandCover.png',bbox_inches='tight')   
plt.close()


# In[25]:
#quit()

# The following line will give you the help notes I've made for the plotmap routine, 
#   they are far from complete but is something at least
#help(PTs.plot_map)

