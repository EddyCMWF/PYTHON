#!/usr/bin/python
#
# Edward Comyn-Platt (edwcom@ceh.ac.uk)
# Routine to plot BVOC comparison
#
#  ./EMEP4UK_plot_BVOCS_comparison.py [INTERACTIVE] [iDISPLAY] [SOURCES] [iTRES] [PLOTS]
#  ./EMEP4UK_plot_BVOCS_comparison.py Y
#  ./EMEP4UK_plot_BVOCS_comparison.py N N [0,1,2] 1 YYYYYY
#  
#   Optional parse inputs:  -start_year YYYY
#                           -end_year YYYY
#   
###################################################
#

import os,sys
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
import data_info_EMEP4UK as data_info
import plot_tools as PT
import datetime as dt
import maths_tools.TimeSeriesTools as TST
import pandas as pd
import matplotlib.pyplot as plt

OUT_DIR='/users/eow/edwcom/BVOCs/PLOTS/ECP_plots/'
iso_mean_range=[0,1]
iso_diff_range=[-0.5,0.5]
terp_mean_range=[0,2]
terp_diff_range=[-1,1]
lonrange=[-13,10.8]
latrange=[51.5,56.8]
TS_latrange=[50,60]
TS_lonrange=[-8,2]
colors=['r','b','g','y','k','c','indigo','darkorange']

# Deal with optional inputs first.
if '-start_year' in sys.argv:
    # Overwriting start_year. otherwise start year is taken from data sources
    argloc = sys.argv.index('-start_year')
    temp   = sys.argv.pop(argloc)
    START_YEAR = int(sys.argv.pop(argloc))
    del temp
else:
    START_YEAR=None

if '-end_year' in sys.argv:
    # Overwriting end_year. otherwise start year is taken from data sources
    argloc = sys.argv.index('-end_year')
    temp   = sys.argv.pop(argloc)
    END_YEAR = int(sys.argv.pop(argloc))
    del temp
else:
    END_YEAR=None

# If no input values provided assume INTERACTIVE='Y'
if len(sys.argv)<2:
    INTERACTIVE='Y'
else:
    INTERACTIVE=sys.argv[1] #'Y'


if INTERACTIVE=='Y':
    iDISPLAY=raw_input('Display images? (Y/N) ')
else:
    iDISPLAY=sys.argv[2]


BVOC_SOURCES = data_info.BVOC_sources()
nSOURCES_ALL = len(BVOC_SOURCES)
print 'Available sources: '
for iSOURCE in range(nSOURCES_ALL):
    print iSOURCE, BVOC_SOURCES[iSOURCE]['NAME']
if INTERACTIVE=='Y':
    SOURCE_NUMS = raw_input('Select sources seperated by commas, for all sources enter ALL: ')
else:
    SOURCE_NUMS = sys.argv[3]
    SOURCE_NUMS = SOURCE_NUMS.replace('[','').replace(']','')

if SOURCE_NUMS=='ALL':
    SOURCE_NUMS=range(nSOURCES)
else:
    SOURCE_NUMS=SOURCE_NUMS.split(',')
    SOURCE_NUMS = [ int(SOURCE_N) for SOURCE_N in SOURCE_NUMS ]

nSOURCES=len(SOURCE_NUMS)
SOURCE_NAMES = [ BVOC_SOURCES[iSOURCE]['NAME'] for iSOURCE in SOURCE_NUMS ]
print SOURCE_NUMS

print 'Available time resolutions:'
Tresolutions=['hourly','daily','monthly']
TRES_list   = [3600,86400,-1]
for iTRES in range(len(Tresolutions)):
    print iTRES, Tresolutions[iTRES]

if INTERACTIVE=='Y':
    iTRES = input('Select a time resolution for comparison: ')
else:
    iTRES = int(sys.argv[4])
TRES = TRES_list[iTRES]
print iTRES, Tresolutions[iTRES], ', '+str(TRES)+' seconds'

PLOTS=''
if INTERACTIVE=='Y':
    #PLOTS[0]
    PLOTS+=raw_input('Produce maps of mean emission rate for entire simulation period? (Y/N) ')
    #PLOTS[1]
    PLOTS+=raw_input('Produce maps of mean emission rate for each year of simulation? (Y/N) ')
    #PLOTS[2]
    PLOTS+=raw_input('Produce maps of mean difference in emission rate for entire simulation period? '+\
                     '(SOURCE[0]-SOURCE[1]) (Y/N) ')
    #PLOTS[3]
    PLOTS+=raw_input('Produce maps of mean difference in emission rate for each year of simulation? '+\
                     '(SOURCE[0]-SOURCE[1]) (Y/N) ')
    #PLOTS[4]
    PLOTS+=raw_input('Produce time series of the mean UK emission rate for simulation? (Y/N) ')
    #PLOTS[5]
    PLOTS+=raw_input('Produce time series of the Total UK annual emission for simulation? (Y/N) ')
    #PLOTS[6]
    PLOTS+=raw_input("Produce regional time series'of the Total annual emission for simulation? (Y/N) ")
    #PLOTS[7]
    PLOTS+=raw_input("Produce regional time series' of the Mean emission rate for simulation? (Y/N) ")
else:
    PLOTS=sys.argv[5]

# print out resubmit command in red text
if INTERACTIVE=='Y':
    print 'Resubmit Command: '
    print '\033[1;31m '+\
            __file__ + ' N ' + \
            iDISPLAY+' '+\
            str(SOURCE_NUMS).replace(' ','')+' '+\
            str(iTRES)+' '+\
            PLOTS+' '+\
            '\033[0m'


iso_data=[]
terp_data=[]
time_data=[]

for iSOURCE in SOURCE_NUMS:
    # Read isoprene data onto 2D grid and convert to same units using special funciton:
    iso_temp,terp_temp,time_vector = \
                    data_info.read_netCDF_BVOC_data_to_TRES( BVOC_SOURCES[iSOURCE],\
                                                             TRES=TRES, \
                                                             return_time_vector=True, \
                                                             start_year=START_YEAR, \
                                                             end_year=END_YEAR,     \
                                                             )
    print(iSOURCE,iso_temp.shape) 
    iso_data.append(iso_temp)
    terp_data.append(terp_temp)
    time_data.append(time_vector)

# read lat/lon data from index file
grid_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc'
grinf=nc.Dataset(grid_file,'r')
grimask=np.ones_like(grinf.variables['land_index'][:])
grinf.close()

# read Country Mask from file
country_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_CountryFile.nc'
cntinf=nc.Dataset(country_file,'r')
CountryMask=cntinf.variables['Country'][:]
CountryNote=cntinf.variables['Country'].note
cntinf.close()

CountryNote=CountryNote.split(',')
CountryList=[ tuple(str(note).split('=')) for note in CountryNote ]

AREA_file='/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_output/rv4.3/'+\
            'EMEP4UK_UK_webrun_emep_4.3_Area_Grid_km2.nc'
AREAinf=nc.Dataset(AREA_file,'r')
lats_2D=AREAinf.variables['lat'][:].squeeze()
lons_2D=AREAinf.variables['lon'][:].squeeze()
AREA_grid=AREAinf.variables['Area_Grid_km2'][:].squeeze()*1000*1000
AREAinf.close()

###################################################################################################
# Plot map of mean emission rate for full simulation
if PLOTS[0]=='Y':
    for jSOURCE in range(nSOURCES): 
        iSOURCE=SOURCE_NUMS[jSOURCE]
        print 'Producing maps of mean emission rates for: '+BVOC_SOURCES[iSOURCE]['NAME']
        start_year,end_year=time_data[jSOURCE][0].year,time_data[jSOURCE][-1].year
        cbar_title= '$g.m^{-2}.yr^{-1}$'
        
        plot_map_data = np.mean(iso_data[jSOURCE],axis=0)
        
        data_range=iso_mean_range
        plot_title='Mean Isoprene Emission Rate - '   \
                    +  BVOC_SOURCES[iSOURCE]['NAME']  \
                    + ' - '+str(start_year)+','+str(end_year)
        plot_title=plot_title.replace('_','-')
        file_name=OUT_DIR+'Isoprene_Mean_Emission_Map_'      \
                   +BVOC_SOURCES[iSOURCE]['NAME']+'_'        \
                   +str(start_year)+'_'+str(end_year)+'.png'        
        PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                    DATA_RANGE=data_range,                          \
                    COLOURS=['beige','greenyellow','darkgreen'],    \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,            \
                    CBAR_ORIENTATION='vertical',                    \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,   \
                    iDISPLAY=iDISPLAY, FILE_PLOT=file_name,         \
                    LATDEL=2., LONDEL=2., RESOLUTION='i',           \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
        plt.close()

        del plot_map_data
        
        plot_map_data = np.mean(terp_data[jSOURCE],axis=0)
        
        data_range=terp_mean_range
        plot_title='Mean Terpene Emission Rate - '   \
                    +  BVOC_SOURCES[iSOURCE]['NAME']  \
                    + ' - '+str(start_year)+','+str(end_year)
        plot_title=plot_title.replace('_','-')
        file_name=OUT_DIR+'Terpene_Mean_Emission_Map_'      \
                   +BVOC_SOURCES[iSOURCE]['NAME']+'_'        \
                   +str(start_year)+'_'+str(end_year)+'.png'        

        PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                    DATA_RANGE=data_range,                          \
                    COLOURS=['beige','greenyellow','darkgreen'],    \
                    INTERPOLATE_COLOURS=True,NLEVELS=11,            \
                    CBAR_ORIENTATION='vertical',                    \
                    WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                    PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,   \
                    iDISPLAY=iDISPLAY, FILE_PLOT=file_name,         \
                    LATDEL=2., LONDEL=2., RESOLUTION='i',           \
                    PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
        plt.close()
        
        del plot_map_data
###################################################################################################

        
###################################################################################################
# Plot map of mean emission rate for each year of simulation
if PLOTS[1]=='Y':
    os.system('mkdir -p '+OUT_DIR+'Annual_Emission_Maps/')
    for jSOURCE in range(nSOURCES):
        iSOURCE=SOURCE_NUMS[jSOURCE]
        print 'Producing maps of mean emission rates for: '+BVOC_SOURCES[iSOURCE]['NAME']
        start_year,end_year=time_data[jSOURCE][0].year,time_data[jSOURCE][-1].year
        cbar_title= '$g.m^{-2}.yr^{-1}$'
        
        for year in range(start_year,end_year+1):
            sp = np.where(time_data[jSOURCE]==dt.datetime(year,1,1,0,0))[0]
            ep = np.where(time_data[jSOURCE]==dt.datetime(year+1,1,1,0,0))[0]
            
            if year==end_year:
                plot_map_data = np.mean(iso_data[jSOURCE][sp:,:],axis=0)
            else:
                plot_map_data = np.mean(iso_data[jSOURCE][sp:ep,:],axis=0)
            data_range=iso_mean_range
            plot_title= 'Mean Isoprene Emission Rate - ' \
                        + BVOC_SOURCES[iSOURCE]['NAME']  \
                        + ' - '+str(year)
            plot_title=plot_title.replace('_','-')
            file_name=OUT_DIR+'Annual_Emission_Maps/Isoprene_Mean_Emission_Map_'      \
                       +BVOC_SOURCES[iSOURCE]['NAME']+'_'        \
                       +str(year)+'.png' 
 
            PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                        DATA_RANGE=data_range,                          \
                        COLOURS=['beige','greenyellow','darkgreen'],    \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,            \
                        CBAR_ORIENTATION='vertical',                    \
                        WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                        PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,   \
                        iDISPLAY=iDISPLAY, FILE_PLOT=file_name,              \
                        LATDEL=2., LONDEL=2., RESOLUTION='i',           \
                        PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
            plt.close()
            
            del plot_map_data

            if year==end_year:
                plot_map_data = np.mean(terp_data[jSOURCE][sp:,:],axis=0)
            else:
                plot_map_data = np.mean(terp_data[jSOURCE][sp:ep,:],axis=0)
            data_range=terp_mean_range
            plot_title= 'Mean Terpene Emission Rate - '   \
                        +  BVOC_SOURCES[iSOURCE]['NAME']   \
                        + ' - '+str(year)
            plot_title=plot_title.replace('_','-')
            file_name=OUT_DIR+'Annual_Terpene_Mean_Emission_Map_'      \
                       +BVOC_SOURCES[iSOURCE]['NAME']+'_'       \
                       +str(year)+'.png'
            PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                        DATA_RANGE=data_range,                          \
                        COLOURS=['beige','greenyellow','darkgreen'],    \
                        INTERPOLATE_COLOURS=True,NLEVELS=11,            \
                        CBAR_ORIENTATION='vertical',                    \
                        WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],                 \
                        PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,  \
                        iDISPLAY=iDISPLAY, FILE_PLOT=file_name,                           \
                        LATDEL=2., LONDEL=2., RESOLUTION='i',                       \
                        PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
            plt.close()

            del plot_map_data
###################################################################################################


###################################################################################################
# Plot map of mean difference in emission rates (SOURCE[0]-SOURCE[1]) for full simulation
if PLOTS[2]=='Y':
    print 'Producing maps of mean emission rate difference for: '+ \
             BVOC_SOURCES[SOURCE_NUMS[0]]['NAME']+ ' - ' +         \
             BVOC_SOURCES[SOURCE_NUMS[1]]['NAME']

    # Get time overlap indexes of the 2 datasets

    overlapA,overlapB=TST.overlap_indexes(time_data[0],time_data[1])
    print(time_data[0])
    print(time_data[1])

    start_year=max(time_data[0][0],time_data[1][0]).year
    end_year  =min(time_data[0][-1],time_data[1][-1]).year
    cbar_title= '$g.m^{-2}.yr^{-1}$'
        
    plot_map_data = np.mean( iso_data[0][overlapA,:] - \
                             iso_data[1][overlapB,:],  \
                             axis=0                    )
    data_range=iso_diff_range
    plot_title='Mean Isoprene Emission Rate Difference - \n '        \
                    +  BVOC_SOURCES[SOURCE_NUMS[0]]['NAME'] +'-'  \
                    +  BVOC_SOURCES[SOURCE_NUMS[1]]['NAME']       \
                    + ' - '+str(start_year)+','+str(end_year)
    plot_title=plot_title.replace('_','-')
    file_name=OUT_DIR+'Isoprene_Mean_Emission_Map_'      \
               +BVOC_SOURCES[SOURCE_NUMS[0]]['NAME']+'_'  \
               +BVOC_SOURCES[SOURCE_NUMS[1]]['NAME']+'_'   \
               +str(start_year)+'_'+str(end_year)+'.png'       
    print(plot_map_data.shape)
    PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                DATA_RANGE=data_range,                          \
                COLOURS=['red','pink','lightgrey','lightblue','blue'],    \
                INTERPOLATE_COLOURS=True,NLEVELS=210,           \
                CBAR_ORIENTATION='vertical',                    \
                WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,   \
                TICK_FORMAT='%0.2f', NTICKS=11,                 \
                iDISPLAY=iDISPLAY, FILE_PLOT=file_name,         \
                LATDEL=2., LONDEL=2., RESOLUTION='i',           \
                PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
    plt.close()


    plot_map_data = np.mean( terp_data[0][overlapA,:] - \
                             terp_data[1][overlapB,:],  \
                             axis=0                    )
    data_range=terp_diff_range
    plot_title='Mean Terpene Emission Rate Difference - \n '        \
                    +  BVOC_SOURCES[SOURCE_NUMS[0]]['NAME'] +'-'  \
                    +  BVOC_SOURCES[SOURCE_NUMS[1]]['NAME']       \
                    + ' - '+str(start_year)+','+str(end_year)
    plot_title=plot_title.replace('_','-')
    file_name=OUT_DIR+'Terpene_Mean_Emission_Map_'      \
               +BVOC_SOURCES[SOURCE_NUMS[0]]['NAME']+'_'  \
               +BVOC_SOURCES[SOURCE_NUMS[1]]['NAME']+'_'   \
               +str(start_year)+'_'+str(end_year)+'.png' 
    PT.plot_map(plot_map_data,lons_2D,lats_2D,                  \
                DATA_RANGE=data_range,                          \
                COLOURS=['red','pink','lightgrey','lightblue','blue'],    \
                INTERPOLATE_COLOURS=True,NLEVELS=210,           \
                CBAR_ORIENTATION='vertical',                    \
                WIDTH=8, HEIGHT=8, FONTSIZES=[12,12,12,18],     \
                PLOT_TITLE=plot_title, CBAR_LABEL=cbar_title,   \
                TICK_FORMAT='%0.2f', NTICKS=11,                 \
                iDISPLAY=iDISPLAY, FILE_PLOT=file_name,         \
                LATDEL=2., LONDEL=2., RESOLUTION='i',           \
                PROJECTION='stere', LON_RANGE=lonrange,LAT_RANGE=latrange)
    plt.close()
###################################################################################################

   
###################################################################################################
# Plot time-series of mean emission rate
if PLOTS[4]=='Y':
    # Geographic Index
    TS_index = (CountryMask>=1) & (CountryMask<=4)
    
    TS_iso_ps_list  = []
    TS_terp_ps_list = []
    iso_tot_ps_list = []
    terp_tot_ps_list = []
    for iSOURCE in range(nSOURCES):
        iso_dat=iso_data[iSOURCE]*grimask
        terp_dat=terp_data[iSOURCE]*grimask
        TS_iso_ps_list.append( pd.Series(np.mean(iso_dat[:,TS_index],axis=1),\
                index=time_data[iSOURCE],  \
                name=SOURCE_NAMES[iSOURCE] ) )

        TS_terp_ps_list.append( pd.Series(np.mean(terp_dat[:,TS_index],axis=1),\
                index=time_data[iSOURCE],  \
                name=SOURCE_NAMES[iSOURCE] ) )


        iso_TS=[]
        terp_TS=[]
        time_TS=[]
        start_year,end_year=time_data[iSOURCE][0].year,time_data[iSOURCE][-1].year
        for year in range(start_year,end_year+1):
            sp = np.where(time_data[iSOURCE]==dt.datetime(year,1,1,0,0))[0]
            ep = np.where(time_data[iSOURCE]==dt.datetime(year+1,1,1,0,0))[0]
            time_TS.append(time_data[iSOURCE][sp][0].year)
            if year==end_year:
                iso_dat=iso_data[iSOURCE][sp:,:]*AREA_grid
                terp_dat=terp_data[iSOURCE][sp:,:]*AREA_grid
            else:
                iso_dat=iso_data[iSOURCE][sp:ep,:]*AREA_grid
                terp_dat=terp_data[iSOURCE][sp:ep,:]*AREA_grid
            iso_TS.append(np.sum(np.mean(iso_dat,axis=0)[TS_index])/1e9)
            terp_TS.append(np.sum(np.mean(terp_dat,axis=0)[TS_index])/1e9)
        iso_tot_ps_list.append( pd.Series(np.array(iso_TS),index=np.array(time_TS),\
                name=SOURCE_NAMES[iSOURCE] ) )
        terp_tot_ps_list.append( pd.Series(np.array(terp_TS),index=np.array(time_TS),\
                name=SOURCE_NAMES[iSOURCE] ) )

    iso_df=pd.concat(TS_iso_ps_list,axis=1)
    terp_df=pd.concat(TS_terp_ps_list,axis=1)

    iso_tot_df=pd.concat(iso_tot_ps_list,axis=1)
    terp_tot_df=pd.concat(terp_tot_ps_list,axis=1)
    
    # Create figure for 2x2 plots
    FIG=plt.figure(figsize=[18,12])
    
    # left column is for time-series of mean emission rate
    AX=FIG.add_subplot(2,2,1)
    iso_df.plot(ax=AX,lw=1.5,legend=False,color=colors)
    AX.set_title('Emission Rates',fontsize=18) 
    AX.set_ylabel('Isoprene ($g.m^{-2}.yr^{-1}$)')
    AX.set_ylim(ymin=0)
    AX=FIG.add_subplot(2,2,3)
    terp_df.plot(ax=AX,lw=1.5,legend=False,color=colors)
    AX.set_ylabel('Terpene ($g.m^{-2}.yr^{-1}$)')
    AX.set_ylim(ymin=0)
    handles,labels=AX.get_legend_handles_labels()
    
    # right column for bar charts of total annual emission
    AX=FIG.add_subplot(2,2,2)
    iso_tot_df.plot(ax=AX,kind='bar',color=colors,legend=False)
    AX.set_title('Total Annual Emission',fontsize=18)
    AX.set_ylabel('Isoprene G$g$')
    AX=FIG.add_subplot(2,2,4)
    terp_tot_df.plot(ax=AX,kind='bar',color=colors,legend=False)
    AX.set_ylabel('Terpene G$g$')

    FIG.legend( handles,labels, loc=8, ncol=nSOURCES)

    FIG.suptitle('BVOC UK Emissions',fontsize=14)
    FIG.savefig(OUT_DIR+'BVOC_UK_Emissions.png', bbox_inches='tight')

    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()
###################################################################################################


###################################################################################################
# Plot time-series of UK annual emission
if PLOTS[5]=='Y':
    print "Annual emission plot been combined with time series PLOTS[4]"

###################################################################################################
# Plot time-series of annual emission for UK countries
if PLOTS[6]=='Y':
    # Geographic Index - Use country file mask
    TS_indexes =  [ (CountryMask>=1) & (CountryMask<=4), \
                    CountryMask==1, \
                    CountryMask==2, \
                    CountryMask==3, \
                    CountryMask==4, \
                    #CountryMask==5, \
                    ]

    TS_index_names = [ 'UK',               \
                       'England',          \
                       'Scotland',         \
                       'Wales',            \
                       'Northern Ireland', \
                       #'Ireland',          \
                       ]

    nTS = len(TS_indexes)
    
    # Create empty list for each TS_index
    iso_df_list=[]
    terp_df_list=[]
    iso_tot_df_list  = []
    terp_tot_df_list = []
    iso_max = 0.
    terp_max = 0.
    for iSOURCE in range(nSOURCES):
        iso_dat=iso_data[iSOURCE]*grimask
        terp_dat=terp_data[iSOURCE]*grimask
        TS_iso_ps_list  = [ [] for iTS in range(nTS) ]
        TS_terp_ps_list = [ [] for iTS in range(nTS) ]
        for iTS in range(nTS):
            TS_iso_ps_list[iTS] = pd.Series( np.mean(iso_dat[:,TS_indexes[iTS]],axis=1),  \
                                              index=time_data[iSOURCE], \
                                              name=TS_index_names[iTS]) 
            TS_terp_ps_list[iTS] = pd.Series( np.mean(terp_dat[:,TS_indexes[iTS]],axis=1),\
                                               index=time_data[iSOURCE], \
                                               name=TS_index_names[iTS]) 
    
        iso_df_list.append(pd.concat(TS_iso_ps_list,axis=1))
        terp_df_list.append(pd.concat(TS_terp_ps_list,axis=1))

        iso_TS, terp_TS = {}, {} 
        for name in TS_index_names:
            iso_TS[name],terp_TS[name]=[],[]

        time_TS=[]
        start_year,end_year=time_data[iSOURCE][0].year,time_data[iSOURCE][-1].year
        for year in range(start_year,end_year+1):
            sp = np.where(time_data[iSOURCE]==dt.datetime(year,1,1,0,0))[0]
            ep = np.where(time_data[iSOURCE]==dt.datetime(year+1,1,1,0,0))[0]
            time_TS.append(time_data[iSOURCE][sp][0].year)
            if year==end_year:
                iso_dat=iso_data[iSOURCE][sp:,:]*AREA_grid
                terp_dat=terp_data[iSOURCE][sp:,:]*AREA_grid
            else:
                iso_dat=iso_data[iSOURCE][sp:ep,:]*AREA_grid
                terp_dat=terp_data[iSOURCE][sp:ep,:]*AREA_grid

            for iTS in range(nTS):
                iso_append_data=np.sum(np.mean(iso_dat,axis=0)[TS_indexes[iTS]])/1e9
                iso_max=max(iso_max,np.max(iso_append_data))
                iso_TS[TS_index_names[iTS]].append(iso_append_data)

                terp_append_data=np.sum(np.mean(terp_dat,axis=0)[TS_indexes[iTS]])/1e9
                terp_max=max(terp_max,np.max(terp_append_data))
                terp_TS[TS_index_names[iTS]].append(terp_append_data)
        
        for iTS in iso_TS:
            iso_TS[iTS]=np.array(iso_TS[iTS])
            terp_TS[iTS]=np.array(terp_TS[iTS])
            
        iso_tot_df_list.append( pd.DataFrame(iso_TS,index=time_TS ).reindex_axis(TS_index_names,axis=1) )
        terp_tot_df_list.append( pd.DataFrame(terp_TS,index=time_TS ).reindex_axis(TS_index_names,axis=1) )
    
    iso_tot_max  = np.ceil(iso_max/10.)*10.
    terp_tot_max = np.ceil(terp_max/10.)*10.
    iso_mean_max= np.ceil(max([ np.max(temp.max()) for temp in iso_df_list ])*2.)/2.
    terp_mean_max= np.ceil(max([ np.max(temp.max()) for temp in terp_df_list ])*2.)/2.

    # Plot time series of emission rates (2xnSOURCES, top row Isoprene, bottom row Terpene)
    FIG=plt.figure(figsize=[10*nSOURCES,10])
    for iSOURCE in range(nSOURCES):
        # plot Isoprene emission time series
        AX=FIG.add_subplot(2,nSOURCES,iSOURCE+1)
        iso_df_list[iSOURCE].plot(ax=AX,lw=1.5,legend=False,color=colors)
        AX.set_ylabel('Mean Isoprene Emission Rate ($g.m^{-2}.yr^{-1}$)',fontsize=14)
        AX.set_ylim([0,iso_mean_max])
        AX.set_title(SOURCE_NAMES[iSOURCE],fontsize=18)
        # Plot Terpene emission time series
        AX=FIG.add_subplot(2,nSOURCES,iSOURCE+nSOURCES+1)
        terp_df_list[iSOURCE].plot(ax=AX,lw=1.5,legend=False,color=colors)
        AX.set_ylabel('Mean Terpene Emission Rate ($g.m^{-2}.yr^{-1}$)',fontsize=14)
        AX.set_ylim([0,terp_mean_max])
    handles,labels=AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nTS)
    FIG.suptitle('BVOC mean emission rate by region',fontsize=30)
    FIG.savefig(OUT_DIR+'BVOC_MeanEmissionRateTimeSeries_Regional.png', bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()

    # Plot bar chart of annual emission (2xnSOURCES, top row Isoprene, bottom row Terpene)
    FIG=plt.figure(figsize=[10*nSOURCES,10])
    for iSOURCE in range(nSOURCES):
        # plot Isoprene emission time series
        AX=FIG.add_subplot(2,nSOURCES,iSOURCE+1)
        iso_tot_df_list[iSOURCE].plot(kind='bar',ax=AX,lw=1.5,legend=False,color=colors)
        AX.set_ylabel('Total Isoprene Emission G$g$',fontsize=14)
        AX.set_ylim([0,iso_tot_max])
        AX.set_title(SOURCE_NAMES[iSOURCE],fontsize=18)
        # Plot Terpene emission time series
        AX=FIG.add_subplot(2,nSOURCES,iSOURCE+nSOURCES+1)
        terp_tot_df_list[iSOURCE].plot(kind='bar',ax=AX,lw=1.5,legend=False,color=colors)
        AX.set_ylabel('Total Terpene Emission G$g$',fontsize=14)
        AX.set_ylim([0,terp_tot_max])
    handles,labels=AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nTS)
    FIG.suptitle('BVOC total annual emission by region',fontsize=30)
    FIG.savefig(OUT_DIR+'BVOC_TotalEmissionTimeSeries_Regional.png', bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()

    # Plot bar chart and time series for each source (2xnSOURCES, top row Time Series, bottom row Bar)
    #    Seperate plot for each species
    TS_DF_LIST=[iso_df_list,terp_df_list]
    TOT_DF_LIST=[iso_tot_df_list,terp_tot_df_list]
    SPECIES=['Isoprene','Terpene']
    TOT_MAXES=[iso_tot_max,terp_tot_max]
    MEAN_MAXES=[iso_mean_max,terp_mean_max]
    for ts_df_list,tot_df_list,species,tot_max,mean_max in \
            zip(TS_DF_LIST,TOT_DF_LIST,SPECIES,TOT_MAXES,MEAN_MAXES):
        FIG=plt.figure(figsize=[10*nSOURCES,10])
        for iSOURCE in range(nSOURCES):
            # plot emission time series
            AX=FIG.add_subplot(2,nSOURCES,iSOURCE+1)
            ts_df_list[iSOURCE].plot(ax=AX,lw=1.5,legend=False,color=colors)
            AX.set_ylabel('Mean Emission Rate ($g.m^{-2}.yr^{-1}$)',fontsize=14)
            AX.set_ylim([0,mean_max])
            AX.set_title(SOURCE_NAMES[iSOURCE],fontsize=18)
            # Plot Total annual emission bar chart
            AX=FIG.add_subplot(2,nSOURCES,iSOURCE+nSOURCES+1)
            tot_df_list[iSOURCE].plot(kind='bar',ax=AX,lw=1.5,legend=False,color=colors)
            AX.set_ylabel('Total Annual Emission G$g$',fontsize=14)
            AX.set_ylim([0,tot_max])
        handles,labels=AX.get_legend_handles_labels()
        FIG.legend( handles,labels, loc=8, ncol=nTS)
        FIG.suptitle(species+' emissions by region',fontsize=30)
        FIG.savefig(OUT_DIR+species+'_EmissionSeriesBar_Regional.png', bbox_inches='tight')
        if iDISPLAY=='Y':
            plt.show()
        else:
            plt.close()

   
###################################################################################################
if PLOTS[6]=='Y':
    print "Mean emission plots have been combined with Total emmisions plots (PLOTS[6])"


