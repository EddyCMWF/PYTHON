#!/bin/env python3.5
import netCDF4 as nc
import numpy as np
#import matplotlib.pylab as plt
#import matplotlib.cm as cm
import sys
import brewer2mpl

import PlotTools.plot_tools as PTs

import iris
from iris_jules import jules
#import iris.plot as iplt
#import iris.quickplot as qplt

#from matplotlib.pyplot import rcParams
#get_ipython().magic('matplotlib inline')
#rcParams['figure.figsize'] = (15.,10)

#Directories, files, and options
WFDEI_DIR='/users/eow/edwcom/BVOCs/JULES_OUTPUT/WFDEI_global/yearly/'
runids = ['J4.6_WFDEI_VG', \
          'J4.6_WFDEI_VG_Guent',\
          'J4.6_WFDEI_VG_tefld_Guent',\
          'J4.6_WFDEI_VG_PHEN',\
          'J4.6_WFDEI_TRIFFID_VG',\
          'J4.6_WFDEI_TRIF_Guent',\
          'J4.6_WFDEI_TRIF_tefld_Guent',\
          ]

run_snames = ['std','Guent','tefld_Guent','PHEN','TRIFFID','TRIF_Guent','TRIF_tefld_Guent']
nruns=len(runids)
OUT_DIR=WFDEI_DIR+'plots/'
#year=1990
year=2000
START_YEAR=1990
END_YEAR=1999
profile = 'monthly'#+str(year)
invars=['isoprene','terpene','acetone','methanol','gpp','lai','tstar_gb']
FRAC_file='/users/eow/edwcom/WFD_EI/J4.6_WFDEI_FRAC.monthly.nc'
PFT_names=['Broadleaf Tree','Needleleaf Tree','C3 Grass','C4 Grass','Shrub']
nPFTs=len(PFT_names)

# Read in fraction of land surface cover and area of each grid cell
FRAC_constraint=iris.Constraint('Fractional cover of each surface type')
FRAC_cube=jules.load(FRAC_file,FRAC_constraint)[0]

AREA_constraint=iris.Constraint('Area of gridcell')
AREA_cube=jules.load(FRAC_file,AREA_constraint)[0][0,:]

PFT_frac_cube=FRAC_cube[0,:nPFTs,:,:]

# Isoprene constraints
iso_constraint=iris.Constraint('PFT isoprene emission flux')
terp_constraint=iris.Constraint('PFT (mono-)terpene emission flux')


# Read the data into a dictionary as a cube
ISO_cubes={}
TERP_cubes={}
for runid in runids:
    #J_fname=WFDEI_DIR+runid+'.'+profile+'.'+str(year)+'.nc'
    J_fnames=[WFDEI_DIR+runid+'.'+profile+'.'+str(YEAR)+'.nc' \
                for YEAR in range(START_YEAR,END_YEAR+1)      ]
    print(J_fnames[0])
    print(J_fnames[-1])
    print(len(J_fnames))
    ISO_cubes[runid]=jules.load(J_fnames,iso_constraint)[0]*PFT_frac_cube
    TERP_cubes[runid]=jules.load(J_fnames,terp_constraint)[0]*PFT_frac_cube


# Calculate annual emissions in Tg yr-1
ISO_AnnEmis_cubes={}
TERP_AnnEmis_cubes={}
for runid in runids:
    cube=ISO_cubes[runid].copy()
    AnnEmis=cube.collapsed('time',iris.analysis.MEAN)
    AnnEmis=AnnEmis*AREA_cube
    AnnEmis.convert_units('Gg yr-1')
    AnnEmis.long_name=runid
    ISO_AnnEmis_cubes[runid]=AnnEmis.copy()
    
    cube=TERP_cubes[runid].copy()
    AnnEmis=cube.collapsed('time',iris.analysis.MEAN)
    AnnEmis=AnnEmis*AREA_cube
    AnnEmis.convert_units('Gg yr-1')
    AnnEmis.long_name=runid
    TERP_AnnEmis_cubes[runid]=AnnEmis.copy()

lons_2d,lats_2d=np.meshgrid(AnnEmis.dim_coords[2].points,AnnEmis.dim_coords[1].points)


# Output total annual emissions to text file
outf=open(OUT_DIR+'statistics_'+str(START_YEAR)+'_'+str(END_YEAR)+'.txt','w')
outf.write('%20s   %40s\n'%('', 'Global Annual Isoprene Emission (Tg)'))
outf.write(('%20s - ' + nruns*'%15s '+'\n') % tuple(['PFT']+run_snames))
for iPFT in range(5):
    outf.write(('%20s -'+nruns*'%15.3f '+'\n')%tuple([PFT_names[iPFT]]+ \
                                      [np.sum(ISO_AnnEmis_cubes[runid].data[iPFT,:])*1e-3 \
                                        for runid in runids]))
            
outf.write(('%20s -'+nruns*'%15.3f '+'\n')%tuple(['Total']+ \
                                  [np.sum(ISO_AnnEmis_cubes[runid].data)*1e-3 \
                                    for runid in runids]))
outf.write('\n'+100*'#'+'\n\n')
outf.write('%20s   %40s\n'%('', 'Global Annual Mono-Terpene Emission (Tg)'))
outf.write(('%20s - ' + nruns*'%15s '+'\n') % tuple(['PFT']+run_snames))
for iPFT in range(5):
    outf.write(('%20s -'+nruns*'%15.3f '+'\n')%tuple([PFT_names[iPFT]]+ \
                                      [np.sum(TERP_AnnEmis_cubes[runid].data[iPFT,:])*1e-3 \
                                        for runid in runids]))
            
outf.write(('%20s -'+nruns*'%15.3f '+'\n')%tuple(['Total']+ \
                                  [np.sum(TERP_AnnEmis_cubes[runid].data)*1e-3 \
                                    for runid in runids]))

#
outf.close()

# Output total annual emissions to text file
outf=open(OUT_DIR+'statistics_update_'+str(START_YEAR)+'_'+str(END_YEAR)+'.txt','w')
outf.write('%20s   %40s\n'%('', 'Global Annual Isoprene Emission (Tg)'))
outf.write(('%20s - ' + nruns*'%15s '+'\n') % tuple(['PFT']+run_snames))
PFT_CFs=np.array([0.82,0.71,0.9,0.81,0.89])
for iPFT in range(5):
    run_CFs=[1]+[PFT_CFs[iPFT] for icnt in range(5)]
    outf.write(('%20s -'+nruns*'%15.3f '+'\n')%tuple([PFT_names[iPFT]]+ \
                                      [np.sum(ISO_AnnEmis_cubes[runid].data[iPFT,:])*1e-3*CF \
                                        for runid,CF in zip(runids,run_CFs)]))
            
#
outf.close()

# Plot Global maps of annual isoprene emissions
Ncols=int(np.ceil(nruns/2.))

CMAP=brewer2mpl.get_map('YlGnBu','Sequential','5',reverse=True).get_mpl_colormap(N=200,gamma=0.6)
plot_titles=PLOT_TITLES=[runid.replace('_',' ') for runid in runids]
plot_LONSLIST=[lons_2d for runid in runids]
plot_LATSLIST=[lats_2d for runid in runids]
ISO_plt_maxes=[150,50,20,10,20]
for iPFT in range(nPFTs):
    plot_DATALIST=[ISO_AnnEmis_cubes[runid][iPFT,:].data for runid in runids]
    PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,    \
                       Ncols=Ncols,Nrows=2, FIGSIZE=(28,10),             \
                       COMMON_CBAR=True, FONTSIZE=30, fraction=0.07, \
                       SUPTITLE='Annual Isoprene Emissions - '+PFT_names[iPFT] ,\
                       CMAP=CMAP,NLEVELS=250, \
                       LATDEL=30,LONDEL=45, LON_RANGE=[-180,180],LAT_RANGE=[-60,90], \
                       DATA_RANGE=[0,ISO_plt_maxes[iPFT]], RESOLUTION='c', FONTSIZES=[15,15,18,18], \
                       PLOT_TITLES=plot_titles, CBAR_LABEL='Isoprene GgC', \
                       FILEPLOT=OUT_DIR+'Global_Isoprene_Emissions_'+PFT_names[iPFT]+'_'+str(START_YEAR)+'_'+str(END_YEAR)+'.png',
                      )

ISOmax=150
plot_DATALIST=[ISO_AnnEmis_cubes[runid].collapsed('pft',iris.analysis.SUM).data for runid in runids]
PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,   \
                   Ncols=Ncols,Nrows=2, FIGSIZE=(28,10),            \
                   COMMON_CBAR=True, FONTSIZE=30, fraction=0.07,\
                   SUPTITLE='Annual Isoprene Emissions - Total' ,\
                   PLOT_TITLES=plot_titles, CBAR_LABEL='Isoprene GgC', \
                   FILEPLOT=OUT_DIR+'Global_Isoprene_Emissions_Total_'+str(START_YEAR)+'_'+str(END_YEAR)+'.png',
                   lCLOSE=True,\
                   CMAP=CMAP,NLEVELS=250, \
                   LATDEL=30,LONDEL=45, LON_RANGE=[-180,180],LAT_RANGE=[-60,90], \
                   DATA_RANGE=[0,ISOmax],RESOLUTION='c', FONTSIZES=[15,15,18,18], \
                  )

# Plot Global maps of annual terpene emissions
TERP_plt_maxes=[100,50,5,5,20]
for iPFT in range(nPFTs):
    plot_DATALIST=[TERP_AnnEmis_cubes[runid][iPFT,:].data for runid in runids]
    PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,    \
                       Ncols=Ncols,Nrows=2, FIGSIZE=(28,10),             \
                       COMMON_CBAR=True, FONTSIZE=30, fraction=0.07, \
                       SUPTITLE='Annual Monoterpene Emissions - '+PFT_names[iPFT] ,\
                       CMAP=CMAP,NLEVELS=250, \
                       LATDEL=30,LONDEL=45, LON_RANGE=[-180,180],LAT_RANGE=[-60,90], \
                       DATA_RANGE=[0,TERP_plt_maxes[iPFT]], RESOLUTION='c', FONTSIZES=[15,15,18,18], \
                       PLOT_TITLES=plot_titles, CBAR_LABEL='Monoterpene GgC', \
                       FILEPLOT=OUT_DIR+'Global_Terpene_Emissions_'+PFT_names[iPFT]+'_'+str(START_YEAR)+'_'+str(END_YEAR)+'.png',
                      )

TERPmax=100
plot_DATALIST=[TERP_AnnEmis_cubes[runid].collapsed('pft',iris.analysis.SUM).data for runid in runids]
PTs.plot_map_multi(plot_DATALIST,plot_LONSLIST,plot_LATSLIST,   \
                   Ncols=Ncols,Nrows=2, FIGSIZE=(28,10),            \
                   COMMON_CBAR=True, FONTSIZE=30, fraction=0.07, \
                   SUPTITLE='Annual Monoterpene Emissions - Total' ,\
                   CMAP=CMAP,NLEVELS=250, \
                   LATDEL=30,LONDEL=45, LON_RANGE=[-180,180],LAT_RANGE=[-60,90], \
                   DATA_RANGE=[0,TERPmax], RESOLUTION='c', FONTSIZES=[15,15,18,18], \
                   PLOT_TITLES=plot_titles, CBAR_LABEL='Monoterpene GgC', \
                   FILEPLOT=OUT_DIR+'Global_Terpene_Emissions_Total_'+str(START_YEAR)+'_'+str(END_YEAR)+'.png',
                  )


