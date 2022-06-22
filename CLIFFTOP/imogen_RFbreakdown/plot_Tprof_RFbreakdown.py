# This plots up the output from inversion with the EBM part of IMOGEN. 
# C. Huntingford (17th March 2017)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys,os
import math
#from matplotlib.patches import Polygon
#from matplotlib.collections import PatchCollection

import brewer2mpl
CMAP=brewer2mpl.get_map('Paired','Qualitative','12')

DATA_DIR='/prj/CLIFFTOP/Tprofiles/'
emiss_invert_running = np.load(DATA_DIR+'emiss_invert_running.npy')
delta_temp_global = np.load(DATA_DIR+'delta_temp_global.npy')
dq_all = np.load(DATA_DIR+'dq_all.npy')
co2_ppm_all = np.load(DATA_DIR+'co2_ppm_all.npy')


n_yr = dq_all.shape[0]
n_gcm = dq_all.shape[1]
yr_plot=np.zeros(n_yr)
for i_yr in range(0, n_yr):
    yr_plot[i_yr]=1850.0+np.float(i_yr)


# read in the RF from csv
rcp3pd_file=DATA_DIR+'RCP3PD_MIDYEAR_RADFORCING.csv'
inf=open(rcp3pd_file)
data_lines=inf.readlines()
inf.close()
headers=data_lines.pop(0)[:-1].split(',')
RF_DICT={ hdr:[] for hdr in headers }
for line in data_lines:
    vals=line[:-1].split(',')
    for hdr,val in zip(headers,vals):
        RF_DICT[hdr].append(float(val))

for hdr in headers:
    RF_DICT[hdr]=np.array(RF_DICT[hdr])

RF_DICT['Q_GHG_NON_CO2']=RF_DICT['GHG_RF']-RF_DICT['CO2_RF']
RF_DICT['Q_NON_GHG_incVOLC']=RF_DICT['TOTAL_INCLVOLCANIC_RF']-RF_DICT['GHG_RF']

mpl.rcParams['axes.labelsize']=12
mpl.rcParams['xtick.labelsize']=12
mpl.rcParams['ytick.labelsize']=12

x_lower=1850
x_upper=2015

fig=plt.figure(figsize=(15.,8))
# Radiative forcing 
ax=plt.subplot(111)
ax.set_xlim(x_lower, x_upper)
ax.set_ylim(-2.0, 3.5)
ax.set_xlabel('Year', fontsize=15)
ax.set_ylabel('Radiative Forcing (W/m2)', fontsize=15)

qco2_polygon_y=np.append(RF_DICT['CO2_RF'],np.zeros_like(RF_DICT['CO2_RF']))
qco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qco2_polygon_x,qco2_polygon_y,c=CMAP.hex_colors[2],label='CO2 RF')
qnonco2_polygon_y=np.append(RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF'],RF_DICT['CO2_RF'][::-1])
qnonco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qnonco2_polygon_x,qnonco2_polygon_y,c=CMAP.hex_colors[-2],label='non CO2 GHG RF')
qsolar_polygon_y=np.append((RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF']), 
                           (RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF']+RF_DICT['SOLAR_RF'])[::-1])
qsolar_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qsolar_polygon_x,qsolar_polygon_y,c=CMAP.hex_colors[1],label='Solar RF')
qtotaer_polygon_y=np.append(RF_DICT['TOTAER_DIR_RF'],  #+RF_DICT['VOLCANIC_ANNUAL_RF'],
                            np.zeros_like(RF_DICT['TOTAER_DIR_RF']))
qtotaer_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qtotaer_polygon_x,qtotaer_polygon_y,c=CMAP.hex_colors[-4],label='Aerosol RF')
qother_polygon_y=np.append(RF_DICT['TOTAER_DIR_RF'],  
        (RF_DICT['TOTAER_DIR_RF']+RF_DICT['Other'])[::-1])
qother_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qother_polygon_x,qother_polygon_y,c=CMAP.hex_colors[4],label='Other RF')
qvolc_polygon_y=np.append(RF_DICT['TOTAER_DIR_RF']+RF_DICT['Other'],  
        (RF_DICT['TOTAER_DIR_RF']+RF_DICT['Other']+RF_DICT['VOLCANIC_ANNUAL_RF'])[::-1])
qvolc_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qvolc_polygon_x,qvolc_polygon_y,c=CMAP.hex_colors[0],label='Other RF')
for i_gcm in range(0, n_gcm):
    RFline=ax.plot(yr_plot, dq_all[:,i_gcm], linewidth=1.0, color='blue')

    qtotaer_infered = RF_DICT['TOTAER_DIR_RF'][85:] + RF_DICT['Other'][85:] +\
            (RF_DICT['TOTAL_INCLVOLCANIC_RF'][85:]-dq_all[:-15,i_gcm])
            #(RF_DICT['TOTAL_ANTHRO_RF'][85:]+RF_DICT['SOLAR_RF'][85:]-dq_all[:-15,i_gcm])
    AERline=ax.plot(yr_plot[:-15], qtotaer_infered , linewidth=1.0, color=CMAP.hex_colors[-3])

ax.plot(RF_DICT['YEAR'],RF_DICT['TOTAL_INCLVOLCANIC_RF'],linewidth=2.0,color='red',label='rcp3PD estiamte')
#ax.plot(RF_DICT['YEAR'],RF_DICT['TOTAL_ANTHRO_RF']+RF_DICT['SOLAR_RF'],\
#        linewidth=2.0,color='red',label='rcp3PD estiamte')

handles,labels=ax.get_legend_handles_labels()
handles+=RFline
handles+=AERline
labels+=['GCM-RF']
labels+=['infered GCM-AER-RF']
fig.legend(handles,labels,loc=8,ncol=8)

fig.suptitle('Breakdown of Radiative Forcing',fontsize=25)

fig.savefig(DATA_DIR+'RF_Breakdown_withVolcanoes.png')

fig=plt.figure(figsize=(15.,8))
# Radiative forcing 
ax=plt.subplot(111)
ax.set_xlim(x_lower, x_upper)
ax.set_ylim(-2.0, 3.5)
ax.set_xlabel('Year', fontsize=15)
ax.set_ylabel('Radiative Forcing (W/m2)', fontsize=15)

qco2_polygon_y=np.append(RF_DICT['CO2_RF'],np.zeros_like(RF_DICT['CO2_RF']))
qco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qco2_polygon_x,qco2_polygon_y,c=CMAP.hex_colors[2],label='CO2 RF')
qnonco2_polygon_y=np.append(RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF'],RF_DICT['CO2_RF'][::-1])
qnonco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qnonco2_polygon_x,qnonco2_polygon_y,c=CMAP.hex_colors[-2],label='non CO2 GHG RF')
qsolar_polygon_y=np.append((RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF']), 
                           (RF_DICT['Q_GHG_NON_CO2']+RF_DICT['CO2_RF']+RF_DICT['SOLAR_RF'])[::-1])
qsolar_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qsolar_polygon_x,qsolar_polygon_y,c=CMAP.hex_colors[1],label='Solar RF')
qtotaer_polygon_y=np.append(RF_DICT['TOTAER_DIR_RF'], 
                            np.zeros_like(RF_DICT['TOTAER_DIR_RF']))
qtotaer_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qtotaer_polygon_x,qtotaer_polygon_y,c=CMAP.hex_colors[-4],label='Aerosol RF')
qother_polygon_y=np.append(RF_DICT['TOTAER_DIR_RF'],  
        (RF_DICT['TOTAER_DIR_RF']+RF_DICT['Other'])[::-1])
qother_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qother_polygon_x,qother_polygon_y,c=CMAP.hex_colors[4],label='Other RF')
for i_gcm in range(0, n_gcm):
    RFline=ax.plot(yr_plot, dq_all[:,i_gcm], linewidth=1.0, color='blue')

    qtotaer_infered = RF_DICT['TOTAER_DIR_RF'][85:] + RF_DICT['Other'][85:] +\
            (RF_DICT['TOTAL_ANTHRO_RF'][85:]+RF_DICT['SOLAR_RF'][85:]-dq_all[:-15,i_gcm])
    AERline=ax.plot(yr_plot[:-15], qtotaer_infered , linewidth=1.0, color=CMAP.hex_colors[-3])

ax.plot(RF_DICT['YEAR'],RF_DICT['TOTAL_ANTHRO_RF']+RF_DICT['SOLAR_RF'],\
        linewidth=2.0,color='red',label='rcp3PD estiamte')

handles,labels=ax.get_legend_handles_labels()
handles+=RFline
handles+=AERline
labels+=['GCM-RF']
labels+=['infered GCM-AER-RF']
fig.legend(handles,labels,loc=8,ncol=8)

fig.suptitle('Breakdown of Radiative Forcing',fontsize=25)

fig.savefig(DATA_DIR+'RF_Breakdown_noVolcanoes.png')


fig=plt.figure(figsize=(15.,8))
# Radiative forcing 
ax=plt.subplot(111)
ax.set_xlim(x_lower, x_upper)
ax.set_ylim(-2.0, 3.5)
ax.set_xlabel('Year', fontsize=15)
ax.set_ylabel('Radiative Forcing (W/m2)', fontsize=15)

baseline=RF_DICT['Q_NON_CO2'].copy()
baseline[baseline<0]=0
qco2_polygon_y=np.append(RF_DICT['CO2_RF']+baseline,np.zeros_like(baseline))
qco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qco2_polygon_x,qco2_polygon_y,c=CMAP.hex_colors[2],label='CO2 RF')

qnonco2_polygon_y=np.append(RF_DICT['Q_NON_CO2'],np.zeros_like(RF_DICT['Q_NON_CO2']))
qnonco2_polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
ax.fill(qnonco2_polygon_x,qnonco2_polygon_y,c=CMAP.hex_colors[7],label='non CO2 RF')

for i_gcm in range(0, n_gcm):
    RFline=ax.plot(yr_plot, dq_all[:,i_gcm], linewidth=1.0, color='blue')
    qnonco2_infered = RF_DICT['Q_NON_CO2'][85:] - \
                      (RF_DICT['TOTAL_INCLVOLCANIC_RF'][85:]-dq_all[:-15,i_gcm]) 
    AERline=ax.plot(yr_plot[:-15], qnonco2_infered , linewidth=1.0, color=CMAP.hex_colors[-3])

ax.plot(RF_DICT['YEAR'],RF_DICT['TOTAL_INCLVOLCANIC_RF'],\
        linewidth=2.0,color='red',label='rcp3PD estiamte')

handles,labels=ax.get_legend_handles_labels()
handles+=RFline
handles+=AERline
labels+=['GCM-RF']
labels+=['infered NON-CO2-RF']
fig.legend(handles,labels,loc=8,ncol=8)

fig.suptitle('Breakdown of Radiative Forcing',fontsize=25)

fig.savefig(DATA_DIR+'RF_Breakdown_CO2_NonCO2.png')


fig=plt.figure(figsize=(15.,8))
# Radiative forcing 
ax=plt.subplot(111)
ax.set_xlim(x_lower, x_upper)
ax.set_ylim(-2.0, 3.5)
ax.set_xlabel('Year', fontsize=15)
ax.set_ylabel('Radiative Forcing (W/m2)', fontsize=15)
ax.tick_params(labelright=True,labelsize=12)

polygon_x=np.append(RF_DICT['YEAR'],RF_DICT['YEAR'][::-1])
baseline=RF_DICT['Q_NON_GHG_incVOLC'].copy()
baseline[baseline<0.]=0.0

qnonco2_polygon_y=np.append((RF_DICT['CO2_RF']+RF_DICT['Q_GHG_NON_CO2']+baseline),\
                            (RF_DICT['CO2_RF']+baseline)[::-1])
ax.fill(polygon_x,qnonco2_polygon_y,c=CMAP.hex_colors[-2],label='rcp3PD-non CO2 GHG')

qco2_polygon_y=np.append(RF_DICT['CO2_RF']+baseline,baseline[::-1])
ax.fill(polygon_x,qco2_polygon_y,c=CMAP.hex_colors[2],label='rcp3PD-CO2')

qnonghg_polygon_y=np.append(RF_DICT['Q_NON_GHG_incVOLC'],\
                            np.zeros_like(RF_DICT['Q_NON_GHG_incVOLC']))
ax.fill(polygon_x,qnonghg_polygon_y,c=CMAP.hex_colors[6],label='rcp3PD-Non GHG')

ax.plot(yr_plot,np.zeros_like(yr_plot),c='k')
for i in range(-5,5):
    ax.plot(yr_plot,np.zeros_like(yr_plot)+i,c='k',ls=':')

for i_gcm in range(0, n_gcm):
    RFline=ax.plot(yr_plot, dq_all[:,i_gcm], linewidth=1.0, color='blue')
    qnonco2_infered = RF_DICT['Q_NON_GHG_incVOLC'][85:] - \
                      (RF_DICT['TOTAL_INCLVOLCANIC_RF'][85:]-dq_all[:-15,i_gcm]) 
    AERline=ax.plot(yr_plot[:-15], qnonco2_infered , linewidth=1.0, color=CMAP.hex_colors[-3])

ax.plot(RF_DICT['YEAR'],RF_DICT['TOTAL_INCLVOLCANIC_RF'],\
        linewidth=2.0,color='red',label='rcp3PD Net estimate')

handles,labels=ax.get_legend_handles_labels()
handles+=RFline
handles+=AERline
labels+=['imogen GCM']
labels+=['infered NON-GHG']
fig.legend(handles,labels,loc=8,ncol=8)

fig.suptitle('Breakdown of Radiative Forcing',fontsize=25)

fig.savefig(DATA_DIR+'RF_Breakdown_GHG_nonGHG.png')




