#!/bin/env python
import sys, os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import PlotTools.plot_tools as PTs

###################################################################
#Directories and Filenames
############################
SOIL_DIR='/users/eow/edwcom/GREENHOUSE/SOIL_PROPERTIES/'
# Data
Data_Dir=SOIL_DIR+'datasets/'
Comp_File=Data_Dir+'Merged_Soil_Composition_CHESSgrid.nc'
BC_file=Data_Dir+'SoilProperties_BC_CHESSgrid.nc'
VG_file=Data_Dir+'SoilProperties_VG_CHESSgrid.nc'

CHESS_landcover_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'
############################
# Output
plot_DIR=SOIL_DIR+'plots/'


Soil_Layer_Thick=np.array([0.1,0.25,0.65,2.0])   #(metres)
Soil_Layer_Depth=np.array([0.1,0.35,1.0,3.0]) #(metres)
nSD=len(Soil_Layer_Thick)


# Open and Read Composition File
COMP_Dict={}
COMPinf=nc.Dataset(Comp_File,'r')
for var in COMPinf.variables:
    COMP_Dict[var]=COMPinf.variables[var][:]
COMPinf.close()

# Open and Read Brooks and Corey File
BC_Dict={}
BCinf=nc.Dataset(BC_file,'r')
for var in BCinf.variables:
    BC_Dict[var]=BCinf.variables[var][:]
BCinf.close()

# Open and Read Van Genuchten File
VG_Dict={}
VGinf=nc.Dataset(VG_file,'r')
for var in VGinf.variables:
    VG_Dict[var]=VGinf.variables[var][:]
VGinf.close()

#read in latlon and landcover data from chess_landcover
LLinf=nc.Dataset(CHESS_landcover_file,'r')
lats=LLinf.variables['lat'][:]
lons=LLinf.variables['lon'][:]
landcover=LLinf.variables['frac'][:]
LLinf.close()


plot_titles=['Depth = '+str(depth)+' m' for depth in Soil_Layer_Depth]
plot_LATSLIST=[lats for iSD in range(nSD)]
plot_LONSLIST=[lons for iSD in range(nSD)]

# Plot Composition maps for 4 layers
# Create lists of lat and lon for plotting in multi_plot
plot_dict = {'sand': { 'DATA_RANGE':[0,100],                             \
                       'CBAR_TITLE':'Sand ($%$)' } ,                       \
             'silt': { 'DATA_RANGE':[0,50],                              \
                       'CBAR_TITLE':'Silt ($%$)' } ,                       \
             'clay': { 'DATA_RANGE':[0,50],                              \
                       'CBAR_TITLE':'Clay ($%$)' } ,                       \
             'carbon': { 'DATA_RANGE':[0,50],                              \
                         'CBAR_TITLE':'Carbon ($%$)' } ,                     \
               }

plot_comp=False
if plot_comp:
    for param in plot_dict:    #['sand','silt','clay','carbon']:
        print(param)
        plot_DATALIST=[COMP_Dict[param][iSD,:] for iSD in range(nSD)]
        Title='Soil Composition - '+param.title()
        FILEPLOT=plot_DIR+'SoilComposition_'+param+'.png'
        COLOURS=['linen','khaki','gold','saddlebrown']
        
        PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,\
                            Ncols=2,Nrows=2, FIGSIZE=(20,24), \
                            COMMON_CBAR=True, PLOT_TITLES=plot_titles, \
                            lCLOSE=True, \
                            RESOLUTION='i', SUPTITLE=Title, \
                            INTERPOLATE_COLOURS=True, COLOURS=COLOURS,NLEVELS=250, \
                            CBAR_LABEL=plot_dict[param]['CBAR_TITLE'], \
                            LATDEL=2,LONDEL=2, LON_RANGE=[-8,2], \
                            DATA_RANGE=plot_dict[param]['DATA_RANGE'],\
                            FONTSIZES=[15,15,18,18], \
                            FILEPLOT=FILEPLOT, \
                           )

plot_dict = {'hcap': { 'DATA_RANGE':[1.1,1.35],                    \
                       'PLOT_SCL_FAC':1e-6,                     \
                       'TITLE':'Heat Capacity',                 \
                       'CBAR_TITLE':'(MJ m$^{-2}$ K$^{-1}$)' }, \
             'hcon': { 'DATA_RANGE':[0.2,0.3],                   \
                       'PLOT_SCL_FAC':1.,                      \
                       'TITLE':'Thermal Conductivity',         \
                       'CBAR_TITLE':'(W m$^{-1}$ K$^{-1}$)' }, \
             'bexp': { 'DATA_RANGE':[1,12],  \
                       'PLOT_SCL_FAC':1.,    \
                       'TITLE':'b exponent', \
                       'CBAR_TITLE':'(-)' }, \
             'satcon': { 'DATA_RANGE':[0,20],              \
                         'PLOT_SCL_FAC':1e3,               \
                         'TITLE':'Hydraulic Conductivity', \
                         'CBAR_TITLE':'(mm s$^{-1}$)' },   \
             'sathh': { 'DATA_RANGE':[0,1],             \
                        'PLOT_SCL_FAC':1,              \
                        'TITLE':'Soil Matric Potential', \
                        'CBAR_TITLE':'(m)' },            \
             'sm_sat': { 'DATA_RANGE':[0.,0.5],                  \
                         'PLOT_SCL_FAC':1,                      \
                         'TITLE':'Soil Moisture at Saturation', \
                         'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },      \
             'sm_crit': { 'DATA_RANGE':[0.,0.5],                      \
                          'PLOT_SCL_FAC':1,                          \
                          'TITLE':'Soil Moisture at Critical Point', \
                          'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },          \
             'sm_wilt': { 'DATA_RANGE':[0.,0.5],                     \
                          'PLOT_SCL_FAC':1,                         \
                          'TITLE':'Soil Moisture at Wilting Point', \
                          'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },         \
             }

BC_plotparams=['bexp','sathh','sm_sat', 'sm_crit','sm_wilt']
BC_plotparams=[]
#for param in plot_dict: # ['hcap', 'sm_crit', 'sm_sat', 'hcon', 'sathh', 'satcon', 'sm_wilt', 'bexp']:
for param in BC_plotparams:
    plot_DATALIST=[BC_Dict[param][iSD,:]*plot_dict[param]['PLOT_SCL_FAC'] for iSD in range(nSD)]
    Title='Brooks and Corey '+plot_dict[param]['TITLE']
    print(param, Title)
    FILEPLOT=plot_DIR+'BC_Property_'+param+'.png'
    MPL_CBAR='RdYlBu_r'
    print(plot_dict[param]['CBAR_TITLE'])
    PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,\
                        Ncols=2,Nrows=2, FIGSIZE=(20,24), \
                        COMMON_CBAR=True, PLOT_TITLES=plot_titles, \
                        lCLOSE=True, \
                        RESOLUTION='i', SUPTITLE=Title, \
                        MPL_CBAR=MPL_CBAR,NLEVELS=250,\
                        CBAR_LABEL=plot_dict[param]['CBAR_TITLE'], \
                        LATDEL=2,LONDEL=2, LON_RANGE=[-8,2], \
                        DATA_RANGE=plot_dict[param]['DATA_RANGE'],\
                        FONTSIZES=[15,15,18,18], \
                        FILEPLOT=FILEPLOT, \
                       )
#quit()

plot_dict = {'hcap': { 'DATA_RANGE':[0.5,1.5],                    \
                       'PLOT_SCL_FAC':1e-6,                     \
                       'TITLE':'Heat Capacity',                 \
                       'CBAR_TITLE':'(MJ m$^{-2}$ K$^{-1}$)' }, \
             'hcon': { 'DATA_RANGE':[0.1,0.3],                \
                       'PLOT_SCL_FAC':1.,                      \
                       'TITLE':'Thermal Conductivity',         \
                       'CBAR_TITLE':'(W m$^{-1}$ K$^{-1}$)' }, \
             'oneovernminusone': { 'DATA_RANGE':[2,11],   \
                                    'PLOT_SCL_FAC':1.,    \
                                    'TITLE':'1/(n-1)',    \
                                    'CBAR_TITLE':'(-)' }, \
             'oneoveralpha': { 'DATA_RANGE':[0,1],          \
                               'PLOT_SCL_FAC':1,            \
                               'TITLE':'1/alpha',           \
                               'CBAR_TITLE':'(m$^{-1}$)' }, \
             'ksat': { 'DATA_RANGE':[0,70],              \
                       'PLOT_SCL_FAC':1e3,               \
                       'TITLE':'Hydraulic Conductivity', \
                       'CBAR_TITLE':'(mm s$^{-1}$)' },   \
             'sm_sat': { 'DATA_RANGE':[0.3,0.8],                  \
                         'PLOT_SCL_FAC':1.,                      \
                         'TITLE':'Soil Moisture at Saturation', \
                         'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },      \
             'sm_crit': { 'DATA_RANGE':[0.3,0.8],                      \
                          'PLOT_SCL_FAC':1.,                          \
                          'TITLE':'Soil Moisture at Critical Point', \
                          'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },          \
             'sm_wilt': { 'DATA_RANGE':[0.3,0.8],                     \
                          'PLOT_SCL_FAC':1.,                         \
                          'TITLE':'Soil Moisture at Wilting Point', \
                          'CBAR_TITLE':'(m$^{3}$ m$^{-3}$)' },         \
             }

VG_plotparams=['hcap','hcon','oneovernminusone','oneoveralpha','ksat','sm_crit','sm_sat','sm_wilt']
VG_plotparams=['sm_crit','sm_sat','sm_wilt']
VG_plotparams=[]

for param in VG_plotparams:
    print(param)
    plot_DATALIST=[VG_Dict[param][iSD,:]*plot_dict[param]['PLOT_SCL_FAC'] for iSD in range(nSD)]
    Title='Van Genuchten '+plot_dict[param]['TITLE']
    print(Title)
    FILEPLOT=plot_DIR+'VG_Property_'+param+'.png'
    MPL_CBAR='RdYlBu_r'
    print(plot_dict[param]['CBAR_TITLE'])
    PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,\
                        Ncols=2,Nrows=2, FIGSIZE=(20,24), \
                        COMMON_CBAR=True, PLOT_TITLES=plot_titles, \
                        lCLOSE=True, \
                        RESOLUTION='i', SUPTITLE=Title, \
                        MPL_CBAR=MPL_CBAR,NLEVELS=250,\
                        CBAR_LABEL=plot_dict[param]['CBAR_TITLE'], \
                        LATDEL=2,LONDEL=2, LON_RANGE=[-8,2], \
                        DATA_RANGE=plot_dict[param]['DATA_RANGE'],\
                        FONTSIZES=[15,15,18,18], \
                        FILEPLOT=FILEPLOT, \
                       )

#Plot VG Soil Textures
param='Soil_Texture'
print(param)
plot_DATALIST=[VG_Dict[param][iSD,:] for iSD in range(nSD)]
Title='Van Genuchten Soil Texture'
print(Title)
FILEPLOT=plot_DIR+'VG_Property_'+param+'.png'
MPL_CBAR='Set1'
PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,\
                   Ncols=2,Nrows=2, FIGSIZE=(20,24), \
                   COMMON_CBAR=True, PLOT_TITLES=plot_titles,NTICKS_COM=7, \
                   lCLOSE=True, \
                   RESOLUTION='i', SUPTITLE=Title, \
                   MPL_CBAR=MPL_CBAR,NLEVELS=7,\
                   CBAR_LABEL='1=CR, 2=MD, 3=MF, 4=FI, 5=VF, 6=OR', \
                   LATDEL=2,LONDEL=2, LON_RANGE=[-8,2], \
                   DATA_RANGE=[1,7],\
                   FONTSIZES=[15,15,18,18], \
                   FILEPLOT=FILEPLOT, \
                   )


