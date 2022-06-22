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

# Directories containing JULES output and plot output directories:
CONFIG= optional_argparse('-config','S3') 
print('Config: ',CONFIG)
DATA_DIR = '/prj/CLIFFTOP/TRENDY/'+CONFIG+'/'

PLOT_DIR=DATA_DIR+'plots/'
os.system('mkdir '+PLOT_DIR)

# Directory of ancillary data:
ANCILS_DIR='/prj/CLIFFTOP/TRENDY/'
# 1Dimension grid cell area data for calculating totals etc.
AREA_file=ANCILS_DIR+'Area.nc'
Ainf=nc.Dataset(AREA_file,'r')
AREA = Ainf.variables['area'][:].squeeze()
lats_1D = Ainf.variables['latitude'][:].squeeze()
lons_1D = Ainf.variables['longitude'][:].squeeze()
Ainf.close()
lons_2D,lats_2D = np.meshgrid(lons_1D,lats_1D)

# Select Scenarios to plot:
SCENARIOs = ['presentday']
SCENARIOyears = [PRESENT_DAY_YEAR]
                                              #  (see filename construction later)
SCENARIO_names=['Present Day ('+str(PRESENT_DAY_YEAR)+')']
nSCENARIOs=len(SCENARIOs)
runid='CRU-NCEP_1p1.update.vn4.6.'   # runid from jules output (i.e. the start)

# File containing the pre-industrial conditions to use as a baseline:
PI_year=1850
PI_scen='1p5deg'

###################################################################################################
# Names of variables to store global total stock in dictionary:
stock_vars = ['CV','CS']  #,'AtmCO2_ppm','AtmCO2_kg','frac','Woody_Products','OceanCO2']

# Variables to store map data:
map_vars = ['CV','CS'] #,'Max_Frac']+[Tname for Tname in Tile_names ]

###################################################################################################
# Read in the Scenario data

# Dictionaries for storing data:
DATA_DICT={ var: [] for var in stock_vars}
PreI_DICT={ var: [] for var in stock_vars}
MAPDATA_DICT={ var: [] for var in map_vars} 
MapPreI_DICT={ var: [] for var in map_vars} 

# Read in Vegetation Carbon which  is just on land points:
CVfile=DATA_DIR+runid+CONFIG+'.Annual.cVeg.nc'
CVinf=nc.Dataset(CVfile,'r')
time_obj = nc.num2date(CVinf.variables['time'][:],
                        units=CVinf.variables['time'].units,
                        calendar=CVinf.variables['time'].calendar)
years = np.array([time.year for time in time_obj])
PD_index = np.where(years==PRESENT_DAY_YEAR)[0]
PD_CV = CVinf.variables['cVeg'][PD_index,:]
PI_CV = CVinf.variables['cVeg'][0,:]
CVinf.close()
MAPDATA_DICT['CV'].append(PD_CV)
# Sum CV over land points
PD_CV = np.sum(PD_CV*AREA)
DATA_DICT['CV'].append(CV)
        
CSfile=DATA_DIR+runid+CONFIG+'.Annual.cSoil.nc'
CSinf = nc.Dataset(CSfile)
PD_CS = CSinf.variables['cSoil'][PD_index,:]
PI_CS = CSinf.variables['cSoil'][0,:]
# Store map data:
MAPDATA_DICT['CS'].append(CS)
# sum over land points:
CS = np.sum( CS*AREA )
#print('CS = ', CS)
DATA_DICT['CS'].append(CS)
            
        
###################################################################################################
# Read in the Pre Industrial Data into it's own dictionaries:
PreI_DICT={}
MapPreI_DICT={}
print('GCM: ',gcm)
DUMP_FILE=DATA_DIR+PI_gcm+'/'+runid+'_'+PI_gcm+'_'+PI_scen+'.dump.'+str(PI_year)+'0101.0.nc'
print(DUMP_FILE)
DINF = nc.Dataset(DUMP_FILE,'r')
Ann_File=DATA_DIR+PI_gcm+'/'+runid+'_'+PI_gcm+'_'+PI_scen+'.Annual_carbon.'+str(PI_year)+'.nc'
print(Ann_File)
Ainf=nc.Dataset(Ann_File,'r')

# Read in CV, land points only:
#CV = DINF.variables['cv'][:]
CV = Ainf.variables['cv'][:].squeeze()
# store map data
MapPreI_DICT['CV']=CV
# total:
CV = np.sum(CV*AREA_1D)
print('CV = ',CV)
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
print('CS = ', CS)
PreI_DICT['CS']=CS

# Fill Woody Products with a dummy zero
PreI_DICT['Woody_Products']=0.0
        
#Atmospheric CO2 
AtmCO2_ppm = DINF.variables['co2_ppmv'][0]
AtmCO2_kg = AtmCO2_ppm*ppm_to_kgC
print('AtmCO2 = ',AtmCO2_kg)
PreI_DICT['AtmCO2_kg']=AtmCO2_kg
PreI_DICT['AtmCO2_ppm']=AtmCO2_ppm
        
# Ocean CO2 from dtemp_o
OCEAN_CO2 = 0
PreI_DICT['OceanCO2']=OCEAN_CO2
print('OCEAN_CO2 = ',OCEAN_CO2)

#FRAC = DINF.variables['frac'][:]
FRAC = Ainf.variables['frac'][:].squeeze()
for iTile in range(nTiles):
    MapPreI_DICT[Tile_names[iTile]]=FRAC[iTile,:]

MAX_FRAC = np.argmax(FRAC,axis=0)
MapPreI_DICT['Max_Frac']=MAX_FRAC
print(MAX_FRAC.shape)
print(MAX_FRAC)
FRAC = np.sum(FRAC*AREA_1D.squeeze()*1e-10,axis=1)  # m^2 to Mha
print(FRAC.shape)
PreI_DICT['frac']=FRAC
print(FRAC)
        
DINF.close()


###################################################################################################
# Bar plot of Absolute Stocks:
# I've tried to formulate the construction of the bar chart a little to prevent overly hard coding
#   Hopefully it's followable

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
        plotdata = np.array(DATA_DICT[scenario][plotvar])*1e-12
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


#fig.legend([cv_bar,cs_bar,atm_bar,ocean_bar],['Vegetation','Soil','Atmosphere','Ocean'],\
#          loc='upper center',ncol=4,fontsize=30)
fig.legend(bar_list,legend_names,     
            loc='upper center',ncol=4,fontsize=30)
fig.savefig(PLOT_DIR+'Equilibrium_CarbonStores.png')  # Store as png
fig.savefig(PLOT_DIR+'Equilibrium_CarbonStores.eps')  # Store as encasulated post script
#plt.show()
plt.close()


############################################################
# plot Delta Cstores

# Open a nice big figure with one axes:
fig,ax=plt.subplots(ncols=1,nrows=1,figsize=[25,10])

plotvar_list = ['CV','CS','AtmCO2_kg','OceanCO2']
plotcolour_list=['darkolivegreen','saddlebrown','khaki','aqua']
legend_names=['Vegetation','Soil','Atmosphere','Ocean']

nbars = len(plotvar_list)  # Number of bars per scenario (i.e. Veg, Soil, Amos)
            #   add ocean and BECCS and Woody products to this at some point
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
        plotdata = (np.array(DATA_DICT[scenario][plotvar])-PreI_DICT[plotvar])*1e-12
        if iscenario == 0:
            bar_list.append(ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width))
        else:
            ax.bar(xpos,np.mean(plotdata),color=plotcolour,width=bar_width)
        #  Plot GCM as line:
        #ax.plot([xpos+bar_width/2,xpos+bar_width/2],[np.min(CV),np.max(CV)],c='k',lw=2)
        #  Plot GCM spread as points:
        ax.plot([xpos+bar_width/2 for i in range(nGCMs)],plotdata,c='k',ls='',marker='.') 

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
#Plot present day CO2 in red:
present_CO2=(400/GtC_to_ppm)-(PreI_DICT['AtmCO2_kg']*1e-12)
print(present_CO2)
ax2.plot(ax2.get_xlim(),[present_CO2,present_CO2],c='r',lw=3,ls='--')


fig.legend(bar_list,legend_names,     
            loc='upper center',ncol=4,fontsize=30)

fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores.png')
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCarbonStores.eps')
#plt.show()
plt.close()
#quit()
###################################################################################################################
#  Plot Bar chart of the Cover Fractions

# Construct a dictionary of the tile cover from the frac
Tile_Dict = { scenario:{tile:[] for tile in Tile_names} for scenario in SCENARIOs }
for scenario in SCENARIOs:
    for iTile in range(nTiles):
        for i_gcm in range(nGCMs):
            Tile_Dict[scenario][Tile_names[iTile]].append(DATA_DICT[scenario]['frac'][i_gcm][iTile])
        Tile_Dict[scenario][Tile_names[iTile]] = np.array(Tile_Dict[scenario][Tile_names[iTile]])
#############################################
# First Absoulute Data
# Nice big figure with axis for each scenario and the preindustrial on different rows:
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
        plotdata=Tile_Dict[scenario][tile]
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

fig.savefig(PLOT_DIR+'Equilibrium_CoverFractions.png')
fig.savefig(PLOT_DIR+'Equilibrium_CoverFractions.eps')
#plt.show()
plt.close()
#quit()

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
    
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCoverFractions.png')
fig.savefig(PLOT_DIR+'Equilibrium_DeltaCoverFractions.eps')
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
        
        FIG.savefig(PLOT_DIR+'Map_'+variable+'_'+scenario+'.png')
        FIG.savefig(PLOT_DIR+'Map_'+variable+'_'+scenario+'.eps')
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
FIG.savefig(PLOT_DIR+'Map_MaxLandCover.png')   
FIG.savefig(PLOT_DIR+'Map_MaxLandCover.eps')   
plt.close()


# In[25]:
#quit()

# The following line will give you the help notes I've made for the plotmap routine, 
#   they are far from complete but is something at least
#help(PTs.plot_map)

