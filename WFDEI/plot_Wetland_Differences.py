#!/bin/env python

#
import netCDF4 as nc
import netcdftime as nctime
import numpy as np
import matplotlib.pyplot as plt
from PlotTools import plot_tools as PT
from GEO_tools import TRANSCOM
import datetime as dt
import maths_tools.DateTimeTools as DTT

start_time=dt.datetime(2000,1,1,0,0,0)
end_time=dt.datetime(2015,1,1,0,0,0)

BASE_DATA_DIR='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JASMIN/WFD_EI_global/'
plot_dir='/users/eow/edwcom/Wetland_Differences/wfdei/Wetlands/'

TRANS_data=TRANSCOM.load_data()
TRANS_names=TRANSCOM.land_region_Snames()

file1=BASE_DATA_DIR+'MAY_2016_results/Jv4.5_WFDEI_nti_TRIFFID_gridded_monthly_ch4.1980_2014.nc'
varname1='fwetl'
timename1='time'

inf1=nc.Dataset(file1,'r')
data1=inf1.variables[varname1][:]
time1=nctime.num2date(inf1.variables[timename1][:],units=inf1.variables[timename1].units)
lats=inf1.variables['latitude'][:]
lons=inf1.variables['longitude'][:]
inf1.close()

nlats=len(lats)
nlons=len(lons)
lons_2d,lats_2d=np.meshgrid(lons,lats)

index1= (time1>=start_time)&(time1<end_time)
#data1.mask[data1.data==0]=True
#data1.data[data1.data==0]=data1.fill_value
data1=data1[index1,:]

file2=BASE_DATA_DIR+'JOEY/PARAMS.nc'
varname2='WETL'
timename2='time'

inf2=nc.Dataset(file2,'r')
data2=inf2.variables[varname2][:]
temp2=inf2.variables['TEMP'][:,:,:]
#time2=nctime.num2date(inf2.variables[timename2][:],units='months since 1993-01-01 00:00:00')
#time2=time1[158:]
time2=np.array(DTT.DTarange_months(dt.datetime(1993,1,1),dt.datetime(2014,12,01) ))
inf2.close()

#data2[:]=np.append(data2[:,:,360:],data2[:,:,:360],axis=2)
index2= (time2>=start_time)&(time2<end_time)
data2=data2[index2,:,:]

#data2=np.ma.masked_array(data2,mask=(data2==0)|(temp2<=0),fill_value=data1.fill_value)
#mask2=data2==0
data2=np.ma.masked_array(data2,mask=data1.mask,fill_value=data1.fill_value)
data2.data[data2.mask==True]=data2.fill_value

time_actual=time1[index1]
ntsteps=len(time_actual)
nyears=ntsteps/12.

mask3=(data1.mask==True)|(data2.mask==True)
ratio_data=np.ma.masked_array(data1/data2,mask=mask3,fill_value=-999.)
diff_data=np.ma.masked_array((data1-data2)/data2,mask=mask3,fill_value=-999.)

climat1=np.mean(data1.reshape(nyears,12,nlats,nlons),axis=0)
climat2=np.mean(data2.reshape(nyears,12,nlats,nlons),axis=0)


#############################################################################################
# plot maps of maximum wetland area
map_datarange=[0,0.3]
nticks=7
maxdata1=np.max(data1,axis=0)
PT.plot_map(maxdata1,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Maximum Wetland Fraction - ECP run', \
            FILE_PLOT=plot_dir+'Global_Max_fwetl_ECP.png', \
            )

maxdata2=np.max(data2,axis=0)
PT.plot_map(maxdata2,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Maximum Wetland Fraction - JRM run', \
            FILE_PLOT=plot_dir+'Global_Max_fwetl_JRM.png', \
            )



#map of mean(data1/data2) ratio
map_datarange=[0.5,1.5]
nticks=11
plotdata=np.median(ratio_data,axis=0)
PT.plot_map(plotdata,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['r','w','b'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Median Wetland Fraction Ratio - ECP/JRM', \
            FILE_PLOT=plot_dir+'Global_Median_fwetl_Ratio.png', \
            )

plotdata=np.mean(ratio_data,axis=0)
PT.plot_map(plotdata,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['r','w','b'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Mean Wetland Fraction Ratio - ECP/JRM', \
            FILE_PLOT=plot_dir+'Global_Mean_fwetl_Ratio.png', \
            )


# map of percentage difference
map_datarange=[-0.5,0.5]
plotdata=np.nanmean(diff_data,axis=0)
PT.plot_map(plotdata,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['r','w','b'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Mean Wetland Fraction Ratio - ECP/JRM', \
            FILE_PLOT=plot_dir+'Global_Mean_fwetl_PercDiff.png', \
            )


zonalmax1=np.mean(maxdata1,axis=1)
zonalmax2=np.mean(maxdata2,axis=1)
plt.plot(zonalmax1,lats,label='ECP (tiwetl=1.5)',lw=2.)
plt.plot(zonalmax2,lats,label='JRM (tiwetl=2.0)',lw=2.)
plt.title('Zonal distribution of maximum wetland',fontsize=24)
plt.legend()
plt.savefig(plot_dir+'Zonal_Max_fwetl_comparison.png')
plt.close()


plt.plot(zonalmax1/zonalmax1.max(),lats,label='ECP (tiwetl=1.5)',lw=2.)
plt.plot(zonalmax2/zonalmax2.max(),lats,label='JRM (tiwetl=2.0)',lw=2.)
plt.title('Normalised Zonal distribution of maximum wetland',fontsize=24)
plt.legend()
plt.savefig(plot_dir+'NormZonal_Max_fwetl_comparison.png')
plt.close()


#MAMM plots
lat_range=[55,75]
MAMM_latdex=(lats>=lat_range[0])&(lats<=lat_range[1])

lon_range=[0,40]
MAMM_londex=(lons>=lon_range[0])&(lons<=lon_range[1])

MAMM_ST=dt.datetime(2012,7,1,0,0,0)
MAMM_ET=dt.datetime(2012,9,1,0,0,0)
MAMM_Tindex = (time_actual>=MAMM_ST)&(time_actual<MAMM_ET)

MAMM_data1=data1[MAMM_Tindex,:,:]
#MAMM_data1=MAMM_data1[:,MAMM_latdex,:]
#MAMM_data1=MAMM_data1[:,:,MAMM_londex]

MAMM_data2=data2[MAMM_Tindex,:,:]
#MAMM_data2=MAMM_data2[:,MAMM_latdex,:]
#MAMM_data2=MAMM_data2[:,:,MAMM_londex]

map_datarange=[0,0.3]
nticks=7
maxdata1=np.max(MAMM_data1,axis=0)
PT.plot_map(maxdata1,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            LAT_RANGE=lat_range,LON_RANGE=lon_range, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Maximum Wetland Fraction - ECP run', \
            FILE_PLOT=plot_dir+'MAMM_Max_fwetl_ECP.png', \
            )

maxdata2=np.max(MAMM_data2,axis=0)
PT.plot_map(maxdata2,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            LAT_RANGE=lat_range,LON_RANGE=lon_range, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Maximum Wetland Fraction - JRM run', \
            FILE_PLOT=plot_dir+'MAMM_Max_fwetl_JRM.png', \
            )





# Read grid file
grid_file='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
grinf=nc.Dataset(grid_file,'r')
grindex=grinf.variables['land_index'][:]-1
grinf.close()

NG_file='/prj/ALANIS/UM_Modelling/EMISSIONS/a_JULES_Gedney/WFDEI79103h.monthly.201007.nc'
NGinf=nc.Dataset(NG_file,'r')
NG_fwetl=NGinf.variables['fwetl'][:].squeeze()
NG_fwetl_2d=np.ma.masked_array(NG_fwetl[grindex],mask=grindex.mask,fill_value=-999.)
NG_fwetl_2d.data[NG_fwetl_2d.mask==True]=NG_fwetl_2d.fill_value

NG_time=dt.datetime(2010,7,1,0,0,0)
ECP_fwetl=data1[time_actual==NG_time,:].squeeze()
JRM_fwetl=data2[time_actual==NG_time,:].squeeze()


# plot maps of  wetland area for 2010-07
map_datarange=[0,0.3]
nticks=7
PT.plot_map(ECP_fwetl,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction  07-2010 - ECP run', \
            FILE_PLOT=plot_dir+'Global_201007_fwetl_ECP.png', \
            )

PT.plot_map(ECP_fwetl,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            LAT_RANGE=lat_range,LON_RANGE=lon_range, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction 07-2010 - ECP run', \
            FILE_PLOT=plot_dir+'MAMM_201007_fwetl_ECP.png', \
            )

PT.plot_map(JRM_fwetl,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction 07-2010 - JRM run', \
            FILE_PLOT=plot_dir+'GLOBAL_201007_fwetl_JRM.png', \
            )

PT.plot_map(JRM_fwetl,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            LAT_RANGE=lat_range,LON_RANGE=lon_range, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction 07-2010 - JRM run', \
            FILE_PLOT=plot_dir+'MAMM_201007_fwetl_JRM.png', \
            )

PT.plot_map(NG_fwetl_2d,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction 07-2010 - NG run', \
            FILE_PLOT=plot_dir+'GLOBAL_201007_fwetl_NG.png', \
            )

PT.plot_map(NG_fwetl_2d,lons_2d,lats_2d,\
            DATA_RANGE=map_datarange, \
            LAT_RANGE=lat_range,LON_RANGE=lon_range, \
            INTERPOLATE_COLOURS=True, \
            COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
            NLEVELS=100, NTICKS=nticks, \
            PLOT_TITLE='Wetland Fraction 07-2010 - NG run', \
            FILE_PLOT=plot_dir+'MAMM_201007_fwetl_NG.png', \
            )


runs=['ECP','JRM']
mov_dats=[data1,data2]
for run,mov_data in zip(runs,mov_dats):
    # Create movie of full fwetl time-series
    for itstep in range(len(time_actual)):
        time_str=str(time_actual[itstep].date())
        fwetl=mov_data[itstep,:]
        PT.plot_map(fwetl,lons_2d,lats_2d,\
                    DATA_RANGE=map_datarange, \
                    INTERPOLATE_COLOURS=True, \
                    COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
                    NLEVELS=20, NTICKS=nticks, \
                    PLOT_TITLE='Wetland Fraction - '+run+' run - '+time_str, \
                    FILE_PLOT=plot_dir+'mpeg_plots/'+run+'_run/Global_'+time_str+'_fwetl_'+run+'.png', \
                    )




# Climatologies
runs=['ECP','JRM']
mov_dats=[climat1,climat2]
map_datarange=[0,0.3]
nticks=7
for run,mov_data in zip(runs,mov_dats):
    # Create movie of full fwetl time-series
    for itstep in range(12):
        time_str='%02d'% (itstep+1)
        fwetl=mov_data[itstep,:]
        PT.plot_map(fwetl,lons_2d,lats_2d,\
                    DATA_RANGE=map_datarange, \
                    INTERPOLATE_COLOURS=True, \
                    COLOURS=['#e1f6fa','#9ee2f0','#09d96c','#def616'], \
                    NLEVELS=20, NTICKS=nticks, TICK_FORMAT='%.2f', \
                    PLOT_TITLE='Wetland Fraction - '+run+' run - '+time_str, \
                    FILE_PLOT=plot_dir+'mpeg_plots/'+run+'_run_climat/Global_'+time_str+'_fwetl_'+run+'.png', \
                    )


runs=['ECP','JRM']
dats=[data1,data2]
cols=['b','g']
#plot TRANSOM time-series
FIG,AXES = plt.subplots(4,3,figsize=[15,15])
AX=(AXES.flat)[0]
AX.set_title('Global')
for run,dat,col in zip(runs,dats,cols):
    TS=np.mean(dat.reshape(ntsteps,-1),axis=1)
    AX.plot(time_actual,TS,c=col,lw=2,label=run)
for iTRANS in range(len(TRANS_names)):
    AX=(AXES.flat)[iTRANS+1]
    AX.set_title(TRANS_names[iTRANS])
    index=TRANS_data==iTRANS+1
    for run,dat,col in zip(runs,dats,cols):
        TS=np.mean(dat[:,index],axis=1)
        AX.plot(time_actual,TS,c=col,lw=2,label=run)

handles,labels=AX.get_legend_handles_labels()
FIG.legend( handles,labels, loc=8, ncol=len(handles),fontsize=16)
FIG.suptitle('Wetland fraction time-series by TRANSCOM region',fontsize=24)
FIG.savefig(plot_dir+'Timeseries_ByTranscomRegion.png',bbox_inxhes='tight')
plt.close()


runs=['ECP','JRM']
dats=[data1.reshape(nyears,12,nlats,nlons),data2.reshape(nyears,12,nlats,nlons)]
cols=['b','g']
#plot TRANSOM time-series
FIG,AXES = plt.subplots(4,3,figsize=[15,15])
AX=(AXES.flat)[0]
AX.set_title('Global')
for run,dat,col in zip(runs,dats,cols):
    TS=np.mean(dat.reshape(nyears,-1),axis=1)
    AX.plot(time_actual[::12],TS,c=col,lw=2,label=run)
for iTRANS in range(len(TRANS_names)):
    AX=(AXES.flat)[iTRANS+1]
    AX.set_title(TRANS_names[iTRANS])
    index=TRANS_data==iTRANS+1
    for run,dat,col in zip(runs,dats,cols):
        TS=np.mean(dat[:,:,index].reshape(nyears,-1),axis=1)
        AX.plot(time_actual[::12],TS,c=col,lw=2,label=run)

handles,labels=AX.get_legend_handles_labels()
FIG.legend( handles,labels, loc=8, ncol=len(handles),fontsize=16)
FIG.suptitle('Wetland fraction time-series by TRANSCOM region',fontsize=24)
FIG.savefig(plot_dir+'AnnualTimeseries_ByTranscomRegion.png',bbox_inxhes='tight')
plt.close()


runs=['ECP','JRM']
dats=[climat1,climat2]
cols=['b','g']
#plot TRANSOM time-series
FIG,AXES = plt.subplots(4,3,figsize=[15,15])
AX=(AXES.flat)[0]
AX.set_title('Global')
AX.set_xlim(1,12)
AX.set_xticks([1,4,7,10])
AX.set_xticklabels(['Jan','Apr','Jul','Oct'])
for run,dat,col in zip(runs,dats,cols):
    TS=np.mean(dat.reshape(12,-1),axis=1)
    AX.plot(np.arange(1,13,dtype='int'),TS,c=col,lw=2,label=run)
for iTRANS in range(len(TRANS_names)):
    AX=(AXES.flat)[iTRANS+1]
    AX.set_title(TRANS_names[iTRANS])
    AX.set_xlim(1,12)
    AX.set_xticks([1,4,7,10])
    AX.set_xticklabels(['Jan','Apr','Jul','Oct'])
    index=TRANS_data==iTRANS+1
    for run,dat,col in zip(runs,dats,cols):
        TS=np.mean(dat[:,index],axis=1)
        AX.plot(np.arange(1,13,dtype='int'),TS,c=col,lw=2,label=run)

handles,labels=AX.get_legend_handles_labels()
FIG.legend( handles,labels, loc=8, ncol=len(handles),fontsize=16)
FIG.suptitle('Wetland fraction climatology by TRANSCOM region',fontsize=24)
FIG.savefig(plot_dir+'Climatology_ByTranscomRegion.png',bbox_inxhes='tight')
plt.close()




