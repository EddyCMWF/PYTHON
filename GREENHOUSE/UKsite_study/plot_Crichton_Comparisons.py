#!/bin/env python
#
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import PlotTools.plot_tools as PT
import netCDF4 as nc
import netcdftime as nctime
import datetime as dt
import glob
import pandas as pd

import PlotTools.plot_tools as PT
#import data_info_UKsites as data_info

INTERACTIVE= sys.argv[1]
C_flux_units='$\mu$mol $m^{-2} s^{-1}$'
kgC_to_umolsCO2_factor = (1./12.)* 1e9

SITE_col=['r']
J_cols=['darkgreen','b','orange','cyan','pink','y','g']

if INTERACTIVE=='Y':
    iDISPLAY=raw_input('Display Images? (Y/N) ')
else:
    iDISPLAY=sys.argv[2]

Managed_N_file='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Mark_Data/'+\
               'Crichton_Managed_plus_AtmosBG_Ndeposition.nc'
N_con_fact = 1e4 * 86400.
inf=nc.Dataset(Managed_N_file,'r')
N_dep_log=inf.variables['N_dep'][:].squeeze()
N_dep_time=nctime.num2date(inf.variables['time'][:], \
                           units=inf.variables['time'].units )
N_dep_pd = pd.Series(N_dep_log,index=N_dep_time)*N_con_fact
N_dep_pd.columns='N_dep'
inf.close()

Management_file='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Mark_Data/'+\
                 'Crichton_Management_Log_Data.nc'
N_con_fact = 1e4 * 86400.
inf=nc.Dataset(Management_file,'r')

Management_dict = { 'N input':inf.variables['N_addition'][:].squeeze()*N_con_fact, \
                    'Grass Harvest':inf.variables['GH_flag'][:].squeeze(), \
                    }
Management_dict['Grass Harvest']=Management_dict['Grass Harvest']*np.max(Management_dict['N input'])
N_dep_time=nctime.num2date(inf.variables['time'][:], \
                           units=inf.variables['time'].units )
Management_pdf = pd.DataFrame(Management_dict,index=N_dep_time)
N_dep_units='$kg$N $ha^{-2}$ $day^{-1}$'
inf.close()

JULES_output_dir='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/output/Crichton/'

SITE_data_dir = '/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Rob_Data/'
SITE_flux_file = SITE_data_dir+'Crichton_Fluxes_ECP.nc'

JULES_drive_dir = '/prj/GREENHOUSE/GREENHOUSE_sites/data/Crichton/Albmar_Data/'
JULES_drive_file = JULES_drive_dir+'drive_crich.nc'

plot_dir='/users/eow/edwcom/GREENHOUSE/GREENHOUSE_sites/plots/Crichton/'

ALL_JULES_runs_names=os.listdir(JULES_output_dir)

print 'Available JULES runs: '
for iJ in range(len(ALL_JULES_runs_names)):
    print str(iJ)+' - '+ALL_JULES_runs_names[iJ]

if INTERACTIVE=='Y':
    JULES_runs_input=raw_input('Select JULES runs for comparison seperated by commas, '+\
                               'for all sims enter ALL: \n')
else:
    JULES_runs_input=sys.argv[3]
#remove square bracket from input string
JULES_runs_input=JULES_runs_input.replace('[','').replace(']','')
if JULES_runs_input=='ALL':
    JULES_runs=range(len(ALL_JULES_runs_names))
else:
    JULES_runs = [ int(Jrun) for Jrun in JULES_runs_input.split(',') ]

nJs=len(JULES_runs)
JULES_runs_names=[ ALL_JULES_runs_names[iJ] for iJ in JULES_runs]

TRES_name_list=['tstep','day']
for iT in range(len(TRES_name_list)):
    print str(iT)+' - '+TRES_name_list[iT]
if INTERACTIVE=='Y':
    iTRES=int(raw_input('Select a time resolution: '))
else: 
    iTRES=int(sys.argv[4])

TRES=TRES_name_list[iTRES]

PLOTS=''
if INTERACTIVE=='Y':
    PLOTS+=raw_input('Produce GPP and NEE time-series plot? (Y/N) ')        # [0]
    PLOTS+=raw_input('Produce GPP and NEE scatter plots? (Y/N) ')           # [1]
    if iTRES==0: 
        PLOTS+=raw_input('Produce dirunal cycle plots? (Y/N) ')             # [2]
    else:
        PLOTS+='N'
    PLOTS+=raw_input('Produce GPP, NEE and TER time-series plot? (Y/N) ')   # [3]
    PLOTS+=raw_input('Produce GPP, NEE and TER scatter plots? (Y/N) ')      # [4]
    if iTRES==0: 
        PLOTS+=raw_input('Produce dirunal cycle plots of GPP, NEE and TER? (Y/N) ')   # [5]
        PLOTS+=raw_input('Produce dirunal cycle plots of resp_s and resp_p? (Y/N) ')  # [6]
    else:
        PLOTS+='N'
    PLOTS+=raw_input('Produce time-series plots of driving met data? (Y/N) ') # [7]
    PLOTS+=raw_input('Produce scatter plots of GPP vs met data? (Y/N) '     ) # [8]
    PLOTS+=raw_input('Produce scatter plots of NEE vs met data? (Y/N) '     ) # [8]
    PLOTS+=raw_input('Produce scatter plots of TER vs met data? (Y/N) '     ) # [8]
else:
    PLOTS=sys.argv[5]

if (INTERACTIVE=='Y'):
    print 'Resubmit Command: '
    #print './plot_PALS_sites_day_pandas.py N '+\
    print '\033[1;31m '+\
            __file__ + ' N ' + \
            iDISPLAY+' '\
            '['+JULES_runs_input+'] '+\
            str(iTRES)+' '+\
            PLOTS+\
            '\033[0m'

# Read in SITE data here
Sinf=nc.Dataset(SITE_flux_file,'r')
S_data={}
S_data['NEE']=Sinf.variables['NEE'][:].squeeze()*-1.
S_data['GPP']=Sinf.variables['GPP'][:].squeeze()
S_data['TER']=Sinf.variables['TER'][:].squeeze()
S_time=nctime.num2date(Sinf.variables['time'][:].squeeze(), \
                       units=Sinf.variables['time'].units  )
Sinf.close()
S_pdf=pd.DataFrame(S_data,index=S_time)
if iTRES==1:
    S_pdf=S_pdf.resample('D')
    
# Read in Drive data 
Dinf=nc.Dataset(JULES_drive_file,'r')
D_data={}
D_data['sw_down']=Dinf.variables['sw_down'][:].squeeze()
D_data['lw_down']=Dinf.variables['lw_down'][:].squeeze()
D_data['pstar']=Dinf.variables['pstar'][:].squeeze()
D_data['precip']=Dinf.variables['precip'][:].squeeze()
D_data['t']=Dinf.variables['t'][:].squeeze()
D_data['q']=Dinf.variables['q'][:].squeeze()
D_time=nctime.num2date(Dinf.variables['time'][:].squeeze(), \
                       units=Dinf.variables['time'].units  )
Dinf.close()
D_pdf=pd.DataFrame(D_data,index=D_time)
if iTRES==1:
    D_pdf=D_pdf.resample('D')
    

Jrun_panda_list=[]
# Read in JULES data
for iJ in range(nJs):
    # constrct filename:
    run_name=JULES_runs_names[iJ]
    Jfname=JULES_output_dir+run_name+'/'+run_name+'_crich.'+TRES+'.nc'
    print Jfname
    Jinf = nc.Dataset(Jfname,'r')
    
    # Run dependent paramter names
    if ('alloff' in run_name) | ('Jvn4.3.1-E-F' not in run_name):
        NPP_diur_name='npp_gb'
        NPP_in_name='npp_gb'
        NPP_out_name='npp_gb'
        resp_s_name='resp_s_gb'
        resp_p_name='resp_p_gb'
    else:
        NPP_diur_name='npp_gb'
        NPP_in_name='npp_nuptake_in_gb'
        NPP_out_name='npp_nuptake_out_gb'
        resp_s_name='co2_soil_gb'
        resp_p_name='resp_p_gb'
    
    J_data={}
    J_data['GPP']=Jinf.variables['gpp_gb'][:].squeeze()*kgC_to_umolsCO2_factor
    J_data['resp_s']=Jinf.variables[resp_s_name][:].squeeze()*kgC_to_umolsCO2_factor
    J_data['resp_p']=Jinf.variables[resp_p_name][:].squeeze()*kgC_to_umolsCO2_factor
    J_data['lai']=Jinf.variables['lai'][:,2,:].squeeze()
    J_data['canht']=Jinf.variables['canht'][:,2,:].squeeze()
    J_data['smcl']=np.sum(Jinf.variables['smcl'][:,:2,:].squeeze(),axis=1)/250.
    J_data['t_soil']=Jinf.variables['t_soil'][:,0,:].squeeze()
    J_data['esoil_gb']=Jinf.variables['esoil_gb'][:].squeeze()

    if (iTRES==0):
        # apply the diurnal cycle from the standard NPP
        J_data['NPP']= (Jinf.variables[NPP_diur_name][:] * \
                        (Jinf.variables[NPP_out_name][:]/Jinf.variables[NPP_in_name][:]) \
                        ).squeeze() * kgC_to_umolsCO2_factor
        J_time=nctime.num2date(np.mean(Jinf.variables['time_bounds'][:].squeeze(),axis=1), \
                               units=Jinf.variables['time'].units  )
    else:
        J_data['NPP']=Jinf.variables[NPP_out_name][:].squeeze()*kgC_to_umolsCO2_factor
        J_time=nctime.num2date(Jinf.variables['time_bounds'][:,0].squeeze(), \
                               units=Jinf.variables['time'].units  )
            
    J_data['NEE'] = J_data['NPP']-J_data['resp_s']
    J_data['TER'] = J_data['GPP']-J_data['NEE']
    
    J_df = pd.DataFrame(J_data,index=J_time)
    Jrun_panda_list.append(J_df)
    Jinf.close()


GPP_pdf=pd.concat([S_pdf['GPP']]+[Jrun_panda['GPP'] for Jrun_panda in Jrun_panda_list],axis=1)
GPP_pdf.columns=['Site']+JULES_runs_names

NPP_pdf=pd.concat([Jrun_panda['NPP'] for Jrun_panda in Jrun_panda_list],axis=1)
NPP_pdf.columns=JULES_runs_names

TER_pdf=pd.concat([S_pdf['TER']]+[Jrun_panda['TER'] for Jrun_panda in Jrun_panda_list],axis=1)
TER_pdf.columns=['Site']+JULES_runs_names

NEE_pdf=pd.concat([S_pdf['NEE']]+[Jrun_panda['NEE'] for Jrun_panda in Jrun_panda_list],axis=1)
NEE_pdf.columns=['Site']+JULES_runs_names

resp_s_pdf=pd.concat([Jrun_panda['resp_s'] for Jrun_panda in Jrun_panda_list],axis=1)
resp_s_pdf.columns=JULES_runs_names

resp_p_pdf=pd.concat([Jrun_panda['resp_p'] for Jrun_panda in Jrun_panda_list],axis=1)
resp_p_pdf.columns=JULES_runs_names

LAI_pdf=pd.concat([Jrun_panda['lai'] for Jrun_panda in Jrun_panda_list],axis=1)
LAI_pdf.columns=JULES_runs_names
LAI_series=LAI_pdf.mean(axis=1)
LAI_series.rename='LAI'

CANHT_pdf=pd.concat([Jrun_panda['canht'] for Jrun_panda in Jrun_panda_list],axis=1)
CANHT_pdf.columns=JULES_runs_names
CANHT_series=CANHT_pdf.mean(axis=1)
CANHT_series.rename='CanHt'

SMCL_pdf=pd.concat([Jrun_panda['smcl'] for Jrun_panda in Jrun_panda_list],axis=1)
SMCL_pdf.columns=JULES_runs_names
SMCL_series=SMCL_pdf.mean(axis=1)
SMCL_series.rename='smcl'

TSOIL_pdf=pd.concat([Jrun_panda['t_soil'] for Jrun_panda in Jrun_panda_list],axis=1)
TSOIL_pdf.columns=JULES_runs_names
TSOIL_series=TSOIL_pdf.mean(axis=1)
TSOIL_series.rename='t_soil'

ESOIL_pdf=pd.concat([Jrun_panda['esoil_gb'] for Jrun_panda in Jrun_panda_list],axis=1)
ESOIL_pdf.columns=JULES_runs_names
ESOIL_series=ESOIL_pdf.mean(axis=1)
ESOIL_series.rename='esoil'


# PLOT SECTION
PLOT_COLS=SITE_col+J_cols[:nJs]#
PLOT_lws = [3]+[2.3 for i in range(nJs)]
PLOT_ls  = ['--']+['-' for i in range(nJs)]
PLOT_names = ['Site']+JULES_runs_names

date_limits=(dt.datetime(2015,4,01),dt.datetime(2015,8,01))

if iTRES==0:
    NEE_Scatter_Range=[-20,30]
else:
    NEE_Scatter_Range=[-5,10]

if PLOTS[0]=='Y':
    FIG,AXES=plt.subplots(nrows=2,ncols=1,figsize=(18,12))
    #FIG.subplots_adjust(right=0.85)
    FIG.subplots_adjust(top=0.93)
    FIG.subplots_adjust(bottom=0.07)
    
    # PLOT GPP Time-Series
    GPP_AX=AXES[0]
    N_AX=AXES[0].twinx()
    for dat,c,lw,ls in zip(GPP_pdf.columns,PLOT_COLS,PLOT_lws,PLOT_ls):
        GPP_pdf[dat].plot(ax=GPP_AX,color=c,legend=False,lw=lw,style=ls)
    
    GPP_AX.text(0.01,0.94,'Mean Bias: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.87,'Std Dev: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.8,"Pearson's $r^2$: ",transform=GPP_AX.transAxes, fontsize=15 )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        GPP_AX.text(0.11+(iJ*0.05),0.94,'%4.2f '%(np.mean(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.11+(iJ*0.05),0.87,'%4.2f '%(np.std(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.11+(iJ*0.05),0.8,'%4.2f '%(GPP_pdf[J_name].corr(GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    # Plot N_deposition in the background
    Management_pdf['N input'].plot(ax=N_AX,color='saddlebrown',lw=3,zorder=2)
    Management_pdf['Grass Harvest'].plot(ax=N_AX,color='chartreuse',lw=3,zorder=1)
    N_dep_pd.plot(ax=N_AX,color='saddlebrown',lw=3,ls=':',zorder=0)
    
    GPP_AX.set_ylabel('GPP ('+C_flux_units+')',fontsize=30)
    GPP_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    GPP_AX.set_ylim( (0,20) ) 
    GPP_AX.patch.set_visible(False)
    GPP_AX.set_zorder(2)

    N_AX.set_ylabel('N deposition ('+N_dep_units+')',fontsize=26)
    N_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    N_AX.set_yticks( (0,15,30,45,60) )
    N_AX.set_frame_on(True)
    N_AX.set_zorder(0)

    ## Plot JULES LAI in the background
    #LAI_AX=AXES[0].twinx()                  
    #LAI_series.plot(ax=LAI_AX,color='lime',lw=2,ls=':',zorder=1)
    #LAI_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    #LAI_AX.set_ylim( (3.5,4) ) 
    #LAI_AX.spines['right'].set_position(('axes',1.05))
    #LAI_AX.set_frame_on(True)
    #LAI_AX.patch.set_visible(False)
    #for sp in LAI_AX.spines.values():
    #    sp.set_visible(False)
    #LAI_AX.spines['right'].set_visible(True)
    #LAI_AX.set_ylabel('LAI',fontsize=20)
    #LAI_AX.set_zorder(1)
    

    # PLOT NEE Time-Series
    NEE_AX=AXES[1]
    N_AX=AXES[1].twinx()
    for dat,c,lw,ls in zip(NEE_pdf.columns,PLOT_COLS,PLOT_lws,PLOT_ls):
        NEE_pdf[dat].plot(ax=NEE_AX,color=c,legend=False,lw=lw,style=ls)
    
    NEE_AX.text(0.01,0.94,'Mean Bias: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.87,'Std Dev: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.8,"Pearson's $r^2$: ",transform=NEE_AX.transAxes, fontsize=15 )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        NEE_AX.text(0.11+(iJ*0.05),0.94,'%4.2f '%(np.mean(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.11+(iJ*0.05),0.87,'%4.2f '%(np.std(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.11+(iJ*0.05),0.8,'%4.2f '%(NEE_pdf[J_name].corr(NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    
    # Plot N_deposition in the background
    Management_pdf['N input'].plot(ax=N_AX,color='saddlebrown',lw=3,zorder=2)
    Management_pdf['Grass Harvest'].plot(ax=N_AX,color='palegreen',lw=3,zorder=1)
    N_dep_pd.plot(ax=N_AX,color='saddlebrown',lw=3,ls='-.',zorder=0)
    
    NEE_AX.set_ylabel('NEE ('+C_flux_units+')',fontsize=30)
    NEE_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2015,07,01)) )
    NEE_AX.set_ylim( (-6,10) )
    NEE_AX.set_yticks( (-6,-2,2,6,10) )
    NEE_AX.patch.set_visible(False)
    NEE_AX.set_zorder(2)

    N_AX.set_ylabel('N deposition ('+N_dep_units+')',fontsize=26)
    N_AX.set_xlim( date_limits ) 
    N_AX.set_yticks( (0,15,30,45,60) )
    N_AX.set_frame_on(True)
    N_AX.set_zorder(0)
    
    ## Plot JULES LAI in the background
    #LAI_AX=AXES[1].twinx()
    #LAI_series.plot(ax=LAI_AX,color='lime',lw=2,ls=':',zorder=1)
    #LAI_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    #LAI_AX.set_ylim( (3.5,4) ) 
    #LAI_AX.spines['right'].set_position(('axes',1.05))
    #LAI_AX.set_frame_on(True)
    #LAI_AX.patch.set_visible(False)
    #for sp in LAI_AX.spines.values():
    #    sp.set_visible(False)
    #LAI_AX.spines['right'].set_visible(True)
    #LAI_AX.set_ylabel('LAI',fontsize=20)
    #LAI_AX.set_zorder(1)
   

    data_handles,data_labels=AXES[0].get_legend_handles_labels()
    N_handles,N_labels=N_AX.get_legend_handles_labels()
    #LAI_handles,LAI_labels=LAI_AX.get_legend_handles_labels()
    handles = list(data_handles)+list(N_handles[:-1])#+list(LAI_handles)
    labels  = list(data_labels)+list(N_labels[:-1])#+['JULES LAI']
    FIG.legend( handles,labels, loc=8, ncol=nJs+4)
    
    FIG.suptitle('Carbon Fluxes at Crichton, Tower vs Model Simulations',fontsize=30)
    FIG.savefig(plot_dir+'/GPP_NEE_timeseries_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


if PLOTS[1]=='Y':
    print 'Plotting scatter of site versus jules simulations'
    FIG= plt.figure(figsize=(25,10))
    GPP_AX = FIG.add_subplot(plt.subplot2grid( (1,2), (0,0) ) )
    NEE_AX = FIG.add_subplot(plt.subplot2grid( (1,2), (0,1) ) )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        GPP_AX.scatter(GPP_pdf['Site'],GPP_pdf[J_name], \
                       c=J_cols[iJ],lw=0,s=25,label=J_name)
        
        GPP_AX.text(0.18+(iJ*0.09),0.96,'%4.2f '%(np.mean(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.18+(iJ*0.09),0.92,'%4.2f '%(np.std(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.18+(iJ*0.09),0.88,'%4.2f '%(GPP_pdf[J_name].corr(GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        NEE_AX.scatter(NEE_pdf['Site'],NEE_pdf[J_name], \
                   c=J_cols[iJ],lw=0,s=25)
        NEE_AX.text(0.18+(iJ*0.09),0.96,'%4.2f '%(np.mean(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.18+(iJ*0.09),0.92,'%4.2f '%(np.std(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.18+(iJ*0.09),0.88,'%4.2f '%(NEE_pdf[J_name].corr(NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    
    # GPP text, axis and 1-to-1 line
    GPP_AX.text(0.01,0.96,'Mean Bias: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.92,'Std Dev: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.88,"Pearson's $r^2$: ",transform=GPP_AX.transAxes, fontsize=15 )
    GPP_limits=[max(min(GPP_AX.get_xlim()[0],GPP_AX.get_ylim()[0],0),0),\
                max(GPP_AX.get_xlim()[1],GPP_AX.get_ylim()[1]) ]
    GPP_AX.plot(np.array(GPP_limits),np.array(GPP_limits),c='k')
    GPP_AX.set_title('Gross Primary Productivity',fontsize=20)
    GPP_AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    GPP_AX.set_xlim(GPP_limits)
    GPP_AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    GPP_AX.set_ylim(GPP_limits)
    GPP_AX.grid(True)
    
    # NEE text, axis and 1-to-1 line
    NEE_AX.text(0.01,0.96,'Mean Bias: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.92,'Std Dev: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.88,"Pearson's $r^2$: ",transform=NEE_AX.transAxes, fontsize=15 )
    NEE_limits=NEE_Scatter_Range
    NEE_AX.plot(np.array(NEE_limits),np.array(NEE_limits),c='k')
    NEE_AX.set_title('Net Ecosystem Exchange',fontsize=20)
    NEE_AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    NEE_AX.set_xlim(NEE_limits)
    NEE_AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    NEE_AX.set_ylim(NEE_limits)
    NEE_AX.grid(True)
    
    # FIGURE - Legend and title and save
    handles,labels=GPP_AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs)
    FIG.suptitle('Carbon Fluxes at Crichton, Tower vs Model Simulations',fontsize=30)
    FIG.savefig(plot_dir+'/GPP_NEE_Scatter_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()

if PLOTS[2]=='Y' and iTRES==0:
    print 'Plotting diurnal cycle data'
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    NEE_diur_pdf=NEE_pdf.dropna(0).copy()
    NEE_diur_pdf['Time']=NEE_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    NEE_diur_pdf = NEE_diur_pdf.groupby('Time').describe().unstack()
    NEE_diur_pdf.index = pd.to_datetime(NEE_diur_pdf.index.astype(str))
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    GPP_diur_pdf=GPP_pdf.dropna(0).copy()
    GPP_diur_pdf['Time']=GPP_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    GPP_diur_pdf = GPP_diur_pdf.groupby('Time').describe().unstack()
    GPP_diur_pdf.index = pd.to_datetime(GPP_diur_pdf.index.astype(str))
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    NPP_diur_pdf=NPP_pdf.dropna(0).copy()
    NPP_diur_pdf['Time']=NPP_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    NPP_diur_pdf = NPP_diur_pdf.groupby('Time').describe().unstack()
    NPP_diur_pdf.index = pd.to_datetime(NPP_diur_pdf.index.astype(str))

    # Plot Diurnal Cycle
    FIG,AXES=plt.subplots( ncols=2, nrows=1, figsize=(20,8) )
    
    for dat,c,lw,ls in zip(PLOT_names,PLOT_COLS,PLOT_lws,PLOT_ls):
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

    AXES[0].set_title('Gross Primary Productivity',fontsize=20)
    AXES[0].set_ylabel('GPP ('+C_flux_units+')',fontsize=20)
    AXES[1].set_title('Net Ecosystem Exchange',fontsize=20)
    AXES[1].set_ylabel('NEE ('+C_flux_units+')',fontsize=20)
    
    FIG.suptitle('Diurnal Cycle of Carbon Fluxes, Crichton',fontsize=30)
    handles,labels=AXES[0].get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs+1)
    FIG.savefig(plot_dir+'/GPP_NEE_Diurnal_Cycle_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()



if PLOTS[3]=='Y':
    FIG,AXES=plt.subplots(nrows=3,ncols=1,figsize=(22,22))
    #FIG.subplots_adjust(right=0.85)
    FIG.subplots_adjust(top=0.93)  
    FIG.subplots_adjust(bottom=0.05)
    # PLOT GPP Time-Series
    GPP_AX=AXES[0]
    N_AX=AXES[0].twinx()
    #LAI_AX=AXES[0].twinx()                  
    for dat,c,lw,ls in zip(GPP_pdf.columns,PLOT_COLS,PLOT_lws,PLOT_ls):
        GPP_pdf[dat].plot(ax=GPP_AX,color=c,legend=False,lw=lw,style=ls)
    
    GPP_AX.text(0.01,0.94,'Mean Bias: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.87,'Std Dev: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.8,"Pearson's $r^2$: ",transform=GPP_AX.transAxes, fontsize=15 )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        GPP_AX.text(0.11+(iJ*0.05),0.94,'%4.2f '%(np.mean(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.11+(iJ*0.05),0.87,'%4.2f '%(np.std(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.11+(iJ*0.05),0.8,'%4.2f '%(GPP_pdf[J_name].corr(GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    # Plot N_deposition in the background
    Management_pdf['N input'].plot(ax=N_AX,color='saddlebrown',lw=2,zorder=2)
    Management_pdf['Grass Harvest'].plot(ax=N_AX,color='chartreuse',lw=2,zorder=1)
    N_dep_pd.plot(ax=N_AX,color='saddlebrown',lw=2,ls=':',zorder=0)
    
    GPP_AX.set_ylabel('GPP ('+C_flux_units+')',fontsize=20)
    GPP_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    GPP_AX.set_ylim( (0,20) ) 
    GPP_AX.patch.set_visible(False)
    GPP_AX.set_zorder(2)

    N_AX.set_ylabel('N deposition ('+N_dep_units+')',fontsize=16)
    N_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    N_AX.set_yticks( (0,15,30,45,60) )
    N_AX.set_frame_on(True)
    N_AX.set_zorder(0)

    ## Plot JULES LAI in the background
    #LAI_series.plot(ax=LAI_AX,color='lime',lw=2,ls=':',zorder=1)
    #LAI_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    #LAI_AX.set_ylim( (3,5) ) 
    #LAI_AX.spines['right'].set_position(('axes',1.05))
    #LAI_AX.set_frame_on(True)
    #LAI_AX.patch.set_visible(False)
    #for sp in LAI_AX.spines.values():
    #    sp.set_visible(False)
    #LAI_AX.spines['right'].set_visible(True)
    #LAI_AX.set_ylabel('LAI',fontsize=20)
    #LAI_AX.set_zorder(1)
    

    # PLOT NEE Time-Series
    NEE_AX=AXES[1]
    N_AX=AXES[1].twinx()
    #LAI_AX=AXES[1].twinx()
    for dat,c,lw,ls in zip(NEE_pdf.columns,PLOT_COLS,PLOT_lws,PLOT_ls):
        NEE_pdf[dat].plot(ax=NEE_AX,color=c,legend=False,lw=lw,style=ls)
    
    NEE_AX.text(0.01,0.94,'Mean Bias: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.87,'Std Dev: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.8,"Pearson's $r^2$: ",transform=NEE_AX.transAxes, fontsize=15 )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        NEE_AX.text(0.11+(iJ*0.05),0.94,'%4.2f '%(np.mean(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.11+(iJ*0.05),0.87,'%4.2f '%(np.std(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.11+(iJ*0.05),0.8,'%4.2f '%(NEE_pdf[J_name].corr(NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    
    # Plot N_deposition in the background
    Management_pdf['N input'].plot(ax=N_AX,color='saddlebrown',lw=2,zorder=2)
    Management_pdf['Grass Harvest'].plot(ax=N_AX,color='palegreen',lw=2,zorder=1)
    N_dep_pd.plot(ax=N_AX,color='saddlebrown',lw=2,ls='-.',zorder=0)
    
    NEE_AX.set_ylabel('NEE ('+C_flux_units+')',fontsize=20)
    NEE_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    NEE_AX.set_ylim( (-6,10) )
    NEE_AX.set_yticks( (-6,-2,2,6,10) )
    NEE_AX.patch.set_visible(False)
    NEE_AX.set_zorder(2)

    N_AX.set_ylabel('N deposition ('+N_dep_units+')',fontsize=16)
    N_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    N_AX.set_yticks( (0,15,30,45,60) )
    N_AX.set_frame_on(True)
    N_AX.set_zorder(0)

    
    ## Plot JULES LAI in the background
    #LAI_series.plot(ax=LAI_AX,color='lime',lw=2,ls=':',zorder=1)
    #LAI_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    #LAI_AX.set_ylim( (3,5) ) 
    #LAI_AX.spines['right'].set_position(('axes',1.05))
    #LAI_AX.set_frame_on(True)
    #LAI_AX.patch.set_visible(False)
    #for sp in LAI_AX.spines.values():
    #    sp.set_visible(False)
    #LAI_AX.spines['right'].set_visible(True)
    #LAI_AX.set_ylabel('LAI',fontsize=20)
    #LAI_AX.set_zorder(1)
    

    # PLOT TER Time-Series
    TER_AX=AXES[2]
    N_AX=AXES[2].twinx()
    #LAI_AX=AXES[2].twinx()
    for dat,c,lw,ls in zip(TER_pdf.columns,PLOT_COLS,PLOT_lws,PLOT_ls):
        TER_pdf[dat].plot(ax=TER_AX,color=c,legend=False,lw=lw,style=ls)
    
    TER_AX.text(0.01,0.94,'Mean Bias: ',transform=TER_AX.transAxes, fontsize=15)
    TER_AX.text(0.01,0.87,'Std Dev: ',transform=TER_AX.transAxes, fontsize=15)
    TER_AX.text(0.01,0.8,"Pearson's $r^2$: ",transform=TER_AX.transAxes, fontsize=15 )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        TER_AX.text(0.11+(iJ*0.05),0.94,'%4.2f '%(np.mean(TER_pdf[J_name]-TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        TER_AX.text(0.11+(iJ*0.05),0.87,'%4.2f '%(np.std(TER_pdf[J_name]-TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        TER_AX.text(0.11+(iJ*0.05),0.8,'%4.2f '%(TER_pdf[J_name].corr(TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    
    # Plot N_deposition in the background
    Management_pdf['N input'].plot(ax=N_AX,color='saddlebrown',lw=2,zorder=2)
    Management_pdf['Grass Harvest'].plot(ax=N_AX,color='palegreen',lw=2,zorder=1)
    N_dep_pd.plot(ax=N_AX,color='saddlebrown',lw=2,ls='-.',zorder=0)
    
    TER_AX.set_ylabel('TER ('+C_flux_units+')',fontsize=20)
    TER_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    TER_AX.set_ylim( (0,16) )
    TER_AX.set_yticks( (0,4,8,12,16) )
    TER_AX.patch.set_visible(False)
    TER_AX.set_zorder(2)

    N_AX.set_ylabel('N deposition ('+N_dep_units+')',fontsize=16)
    N_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    N_AX.set_yticks( (0,15,30,45,60) )
    N_AX.set_frame_on(True)
    N_AX.set_zorder(0)
    
    ## Plot JULES LAI in the background
    #LAI_series.plot(ax=LAI_AX,color='lime',lw=2,ls=':',zorder=1)
    #LAI_AX.set_xlim( date_limits ) #(dt.datetime(2015,01,01),dt.datetime(2016,01,01)) )
    #LAI_AX.set_ylim( (3,5) ) 
    #LAI_AX.spines['right'].set_position(('axes',1.05))
    #LAI_AX.set_frame_on(True)
    #LAI_AX.patch.set_visible(False)
    #for sp in LAI_AX.spines.values():
    #    sp.set_visible(False)
    #LAI_AX.spines['right'].set_visible(True)
    #LAI_AX.set_ylabel('LAI',fontsize=20)
    #LAI_AX.set_zorder(1)
    

    data_handles,data_labels=AXES[0].get_legend_handles_labels()
    N_handles,N_labels=N_AX.get_legend_handles_labels()
    handles = list(data_handles)+list(N_handles[:-1])#+list(LAI_handles)
    labels  = list(data_labels)+list(N_labels[:-1])#+['JULES LAI']
    FIG.legend( handles,labels, loc=8, ncol=nJs+3)
    
    FIG.suptitle('Carbon Fluxes at Crichton, Tower vs Model Simulations',fontsize=40)
    FIG.savefig(plot_dir+'/GPP_NEE_TER_timeseries_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()



if PLOTS[4]=='Y':
    print 'Plotting scatter of site versus jules simulations'
    FIG= plt.figure(figsize=(37,10))
    GPP_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,0) ) )
    NEE_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,1) ) )
    TER_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,2) ) )
    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        GPP_AX.scatter(GPP_pdf['Site'],GPP_pdf[J_name], \
                       c=J_cols[iJ],lw=0,s=25,label=J_name)
        
        GPP_AX.text(0.18+(iJ*0.09),0.96,'%4.2f '%(np.mean(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.18+(iJ*0.09),0.92,'%4.2f '%(np.std(GPP_pdf[J_name]-GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        GPP_AX.text(0.18+(iJ*0.09),0.88,'%4.2f '%(GPP_pdf[J_name].corr(GPP_pdf['Site'])), \
                     transform=GPP_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        NEE_AX.scatter(NEE_pdf['Site'],NEE_pdf[J_name], \
                   c=J_cols[iJ],lw=0,s=25)
        NEE_AX.text(0.18+(iJ*0.09),0.96,'%4.2f '%(np.mean(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.18+(iJ*0.09),0.92,'%4.2f '%(np.std(NEE_pdf[J_name]-NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        NEE_AX.text(0.18+(iJ*0.09),0.88,'%4.2f '%(NEE_pdf[J_name].corr(NEE_pdf['Site'])), \
                     transform=NEE_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
        TER_AX.scatter(TER_pdf['Site'],TER_pdf[J_name], \
                   c=J_cols[iJ],lw=0,s=25)
        TER_AX.text(0.18+(iJ*0.09),0.96,'%4.2f '%(np.mean(TER_pdf[J_name]-TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        TER_AX.text(0.18+(iJ*0.09),0.92,'%4.2f '%(np.std(TER_pdf[J_name]-TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
        TER_AX.text(0.18+(iJ*0.09),0.88,'%4.2f '%(TER_pdf[J_name].corr(TER_pdf['Site'])), \
                     transform=TER_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )
    
    # GPP text, axis and 1-to-1 line
    GPP_AX.text(0.01,0.96,'Mean Bias: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.92,'Std Dev: ',transform=GPP_AX.transAxes, fontsize=15)
    GPP_AX.text(0.01,0.88,"Pearson's $r^2$: ",transform=GPP_AX.transAxes, fontsize=15 )
    GPP_limits=[max(min(GPP_AX.get_xlim()[0],GPP_AX.get_ylim()[0],0),0),\
                max(GPP_AX.get_xlim()[1],GPP_AX.get_ylim()[1]) ]
    GPP_AX.plot(np.array(GPP_limits),np.array(GPP_limits),c='k')
    GPP_AX.set_title('Gross Primary Productivity',fontsize=20)
    GPP_AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    GPP_AX.set_xlim(GPP_limits)
    GPP_AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    GPP_AX.set_ylim(GPP_limits)
    GPP_AX.grid(True)
    
    # NEE text, axis and 1-to-1 line
    NEE_AX.text(0.01,0.96,'Mean Bias: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.92,'Std Dev: ',transform=NEE_AX.transAxes, fontsize=15)
    NEE_AX.text(0.01,0.88,"Pearson's $r^2$: ",transform=NEE_AX.transAxes, fontsize=15 )
    NEE_limits=NEE_Scatter_Range
    NEE_AX.plot(np.array(NEE_limits),np.array(NEE_limits),c='k')
    NEE_AX.set_title('Net Ecosystem Exchange',fontsize=20)
    NEE_AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    NEE_AX.set_xlim(NEE_limits)
    NEE_AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    NEE_AX.set_ylim(NEE_limits)
    NEE_AX.grid(True)
    
    # TER text, axis and 1-to-1 line
    TER_AX.text(0.01,0.96,'Mean Bias: ',transform=TER_AX.transAxes, fontsize=15)
    TER_AX.text(0.01,0.92,'Std Dev: ',transform=TER_AX.transAxes, fontsize=15)
    TER_AX.text(0.01,0.88,"Pearson's $r^2$: ",transform=TER_AX.transAxes, fontsize=15 )
    TER_limits=[max(min(TER_AX.get_xlim()[0],TER_AX.get_ylim()[0],0),0),\
                max(TER_AX.get_xlim()[1],TER_AX.get_ylim()[1]) ]
    TER_AX.plot(np.array(TER_limits),np.array(TER_limits),c='k')
    TER_AX.set_title('Total Ecosystem Respiration',fontsize=20)
    TER_AX.set_xlabel('Site ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    TER_AX.set_xlim(TER_limits)
    TER_AX.set_ylabel('JULES ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    TER_AX.set_ylim(TER_limits)
    TER_AX.grid(True)
    
    # FIGURE - Legend and title and save
    handles,labels=GPP_AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs)
    FIG.suptitle('Carbon Fluxes at Crichton, Tower vs Model Simulations',fontsize=30)
    FIG.savefig(plot_dir+'/GPP_NEE_TER_Scatter_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()




if PLOTS[5]=='Y' and iTRES==0:
    print 'Plotting diurnal cycle data'
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    NEE_diur_pdf=NEE_pdf.dropna(0).copy()
    NEE_diur_pdf['Time']=NEE_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    NEE_diur_pdf = NEE_diur_pdf.groupby('Time').describe().unstack()
    NEE_diur_pdf.index = pd.to_datetime(NEE_diur_pdf.index.astype(str))
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    GPP_diur_pdf=GPP_pdf.dropna(0).copy()
    GPP_diur_pdf['Time']=GPP_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    GPP_diur_pdf = GPP_diur_pdf.groupby('Time').describe().unstack()
    GPP_diur_pdf.index = pd.to_datetime(GPP_diur_pdf.index.astype(str))
    
    # Create Diurnal Panda Data Frames for NEE and TER
    TER_diur_pdf=TER_pdf.dropna(0).copy()
    TER_diur_pdf['Time']=TER_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    TER_diur_pdf = TER_diur_pdf.groupby('Time').describe().unstack()
    TER_diur_pdf.index = pd.to_datetime(TER_diur_pdf.index.astype(str))

    # Create Diurnal Panda Data Frames for NEE and GPP
    NPP_diur_pdf=NPP_pdf.dropna(0).copy()
    NPP_diur_pdf['Time']=NPP_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    NPP_diur_pdf = NPP_diur_pdf.groupby('Time').describe().unstack()
    NPP_diur_pdf.index = pd.to_datetime(NPP_diur_pdf.index.astype(str))

    # Plot Diurnal Cycle
    FIG,AXES=plt.subplots( ncols=3, nrows=1, figsize=(30,8) )
    
    for dat,c,lw,ls in zip(PLOT_names,PLOT_COLS,PLOT_lws,PLOT_ls):
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[0].plot(GPP_diur_pdf.index,GPP_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[1].plot(NEE_diur_pdf.index,NEE_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

        AXES[2].plot(TER_diur_pdf.index,TER_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[2].plot(TER_diur_pdf.index,TER_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[2].plot(TER_diur_pdf.index,TER_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

    AXES[0].set_title('Gross Primary Productivity',fontsize=20)
    AXES[0].set_ylabel('GPP ('+C_flux_units+')',fontsize=20)
    AXES[1].set_title('Net Ecosystem Exchange',fontsize=20)
    AXES[1].set_ylabel('NEE ('+C_flux_units+')',fontsize=20)
    AXES[2].set_title('Total Ecosystem Respiration',fontsize=20)
    AXES[2].set_ylabel('TER ('+C_flux_units+')',fontsize=20)
    
    FIG.suptitle('Diurnal Cycle of Carbon Fluxes, Crichton',fontsize=30)
    handles,labels=AXES[0].get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs+1)
    FIG.savefig(plot_dir+'/GPP_NEE_TER_Diurnal_Cycle_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


if PLOTS[6]=='Y' and iTRES==0:
    print 'Plotting diurnal cycle of plant and soil respiration'
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    resp_s_diur_pdf=resp_s_pdf.dropna(0).copy()
    resp_s_diur_pdf['Time']=resp_s_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    resp_s_diur_pdf = resp_s_diur_pdf.groupby('Time').describe().unstack()
    resp_s_diur_pdf.index = pd.to_datetime(resp_s_diur_pdf.index.astype(str))
    
    # Create Diurnal Panda Data Frames for NEE and GPP
    resp_p_diur_pdf=resp_p_pdf.dropna(0).copy()
    resp_p_diur_pdf['Time']=resp_p_diur_pdf.index.map(lambda x: x.strftime("%H:%M"))
    resp_p_diur_pdf = resp_p_diur_pdf.groupby('Time').describe().unstack()
    resp_p_diur_pdf.index = pd.to_datetime(resp_p_diur_pdf.index.astype(str))
    
    # Plot Diurnal Cycle
    FIG,AXES=plt.subplots( ncols=2, nrows=1, figsize=(20,8) )
    
    for dat,c,lw,ls in zip(PLOT_names[1:],PLOT_COLS[1:],PLOT_lws[1:],PLOT_ls[1:]):
        AXES[0].plot(resp_s_diur_pdf.index,resp_s_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[0].plot(resp_s_diur_pdf.index,resp_s_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[0].plot(resp_s_diur_pdf.index,resp_s_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

        AXES[1].plot(resp_p_diur_pdf.index,resp_p_diur_pdf[dat]['50%'],color=c,label=dat,lw=lw+1)
        AXES[1].plot(resp_p_diur_pdf.index,resp_p_diur_pdf[dat]['25%'],color=c,lw=lw,ls=':')
        AXES[1].plot(resp_p_diur_pdf.index,resp_p_diur_pdf[dat]['75%'],color=c,lw=lw,ls=':')

    AXES[0].set_title('Soil Respiration',fontsize=20)
    AXES[0].set_ylabel('resp_s ('+C_flux_units+')',fontsize=20)
    AXES[0].set_ylim( (0,5) )
    AXES[1].set_title('Plant Respiration',fontsize=20)
    AXES[1].set_ylabel('resp_p ('+C_flux_units+')',fontsize=20)
    AXES[1].set_ylim( (0,5) )
    
    FIG.suptitle('Diurnal Cycle of Respiration, Crichton',fontsize=30)
    handles,labels=AXES[0].get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs)
    FIG.savefig(plot_dir+'/resp_s_resp_p_Diurnal_Cycle_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()

if PLOTS[7]=='Y':
    FIG,AXES=plt.subplots(nrows=1,ncols=1,figsize=(20,7))
    FIG.subplots_adjust(right=0.88,left=0.12,top=0.9,bottom=0.15)
    Precip_AX=AXES
    T_AX=AXES.twinx()
    RAD_AX=AXES.twinx()
    Q_AX=AXES.twinx()
    LAI_AX=AXES.twinx()

    tick_locs=[ np.where( D_pdf[date_limits[0]:date_limits[1]].index == \
                          dt.datetime(2015,mnth,01) )[0] \
                for mnth in range(date_limits[0].month,date_limits[1].month+1) ]
    tick_labels=[ dt.datetime(2015,mnth,01).strftime('%b') \
                  for mnth in range(date_limits[0].month,date_limits[1].month+1) ]

    (D_pdf['precip']*86400.)[date_limits[0]:date_limits[1]].plot(ax=Precip_AX,color='dodgerblue',\
                                                                 legend=False,kind='bar')
    Precip_AX.set_ylabel('Precip. (mm hr$^{-1}$)',fontsize=20)
    Precip_AX.yaxis.tick_right()
    Precip_AX.yaxis.set_label_position('right')
    Precip_AX.set_yticks( (0,8,16,24,32,40) )
    Precip_AX.xaxis.set_ticks(tick_locs)
    Precip_AX.xaxis.set_ticklabels(tick_labels)
    Precip_AX.xaxis.set_visible(False)
    Precip_AX.patch.set_visible(False)
    Precip_AX.set_zorder(1)

    Tplotdata=D_pdf['t'][date_limits[0]:date_limits[1]].values-273.15
    T_AX.plot(Tplotdata,color='brown',label='Air T',lw=2)
    T_AX.set_ylabel('T ($^{o}C$)',fontsize=20)
    T_AX.set_yticks( (0,5,10,15,20,25) )
    T_AX.set_zorder(3)
    T_AX.yaxis.tick_left()
    T_AX.yaxis.set_label_position('left')
    T_AX.xaxis.set_visible(True)

    Qplotdata=D_pdf['q'][date_limits[0]:date_limits[1]].values*1000.0
    Q_AX.plot(Qplotdata,color='seagreen',label='Spec. Hum.',lw=2)
    Q_AX.yaxis.tick_left()
    #Q_AX.set_yticks( (0,3,6,9,12,15) )
    Q_AX.spines['left'].set_position(('axes',-.07))
    Q_AX.set_frame_on(True)
    Q_AX.patch.set_visible(False)
    for sp in Q_AX.spines.values():
        sp.set_visible(True)
    Q_AX.spines['right'].set_visible(True)
    Q_AX.yaxis.set_label_position('left')
    Q_AX.set_ylabel('Spec. Hum. (g kg$^{-1}$)',fontsize=20)
    Q_AX.set_zorder(2)

    SWplotdata=D_pdf['sw_down'][date_limits[0]:date_limits[1]].values
    LWplotdata=D_pdf['lw_down'][date_limits[0]:date_limits[1]].values
    TOTRAD = SWplotdata#+LWplotdata
    RAD_AX.plot(TOTRAD,color='gold',label='SW down',lw=1.7)
    RAD_AX.fill_between(np.arange(len(TOTRAD)),TOTRAD,\
                            where=TOTRAD>np.zeros_like(TOTRAD),\
                            color='khaki')
    RAD_AX.set_ylim(0,RAD_AX.get_ylim()[1])
    #RAD_AX.set_yticks( (400,460,520,580,640,700) )
    RAD_AX.spines['right'].set_position(('axes',1.07))
    RAD_AX.set_frame_on(True)
    RAD_AX.patch.set_visible(True)
    for sp in RAD_AX.spines.values():
        sp.set_visible(False)
    RAD_AX.spines['right'].set_visible(True)
    RAD_AX.set_ylabel('Radiation (W m$^{2}$)',fontsize=20)
    RAD_AX.set_zorder(0)
    
    T_handles,T_labels=T_AX.get_legend_handles_labels()
    Q_handles,Q_labels=Q_AX.get_legend_handles_labels()
    Precip_handles,Precip_labels=Precip_AX.get_legend_handles_labels()
    RAD_handles,RAD_labels=RAD_AX.get_legend_handles_labels()
    handles = list(T_handles)+list(Q_handles)+list(Precip_handles)+list(RAD_handles)
    labels  = list(T_labels) +list(Q_labels) +list(Precip_labels) +list(RAD_labels) 
    FIG.legend( handles,labels, loc=8, ncol=nJs+4)
    
    FIG.suptitle('JULES Driving Met Data for Crichton Flux Tower',fontsize=30)
    FIG.savefig(plot_dir+'/JULES_DrivingData_timeseries_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()

if PLOTS[8]=='Y':
    print 'Plotting scatter of GPP versus diagnostics data'
    FIG= plt.figure(figsize=(24,7))
    FIG.subplots_adjust(bottom=0.15,top=0.85,left=0.05,right=0.98)
    TSOIL_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,0) ) )
    SMCL_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,1) ) )
    LAI_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,2) ) )

    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        TSOIL_AX.scatter(TSOIL_pdf[date_limits[0]:date_limits[1]][J_name]-273.15,\
                        GPP_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25,label=J_name)
        TSOIL_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%\
                        (TSOIL_pdf[date_limits[0]:date_limits[1]][J_name].corr(   \
                                GPP_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=TSOIL_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        SMCL_AX.scatter(SMCL_pdf[date_limits[0]:date_limits[1]][J_name],\
                        GPP_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25)
        SMCL_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%\
                        (SMCL_pdf[date_limits[0]:date_limits[1]][J_name].corr(    \
                                GPP_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=SMCL_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        LAI_AX.scatter(LAI_pdf[date_limits[0]:date_limits[1]][J_name],\
                        GPP_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25)
        LAI_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%(\
                        LAI_pdf[date_limits[0]:date_limits[1]][J_name].corr(\
                            GPP_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=LAI_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

    TSOIL_AX.scatter(TSOIL_series-273.15,GPP_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,label='Site',marker='x')
    TSOIL_AX.text(0.05,0.88,'%4.2f '%(TSOIL_series.corr(GPP_pdf['Site'])), \
                  transform=TSOIL_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    SMCL_AX.scatter(SMCL_series,GPP_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,marker='x')
    SMCL_AX.text(0.05,0.88,'%4.2f '%(SMCL_series.corr(GPP_pdf['Site'])), \
                  transform=SMCL_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    LAI_AX.scatter(LAI_series,GPP_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,marker='x')
    LAI_AX.text(0.05,0.88,'%4.2f '%(LAI_series.corr(GPP_pdf['Site'])), \
                  transform=LAI_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    
    # GPP text, axis and 1-to-1 line
    GPP_limits=[0, TSOIL_AX.get_ylim()[1]]
    TSOIL_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=TSOIL_AX.transAxes, fontsize=15 )
    TSOIL_AX.grid(True)
    TSOIL_AX.set_title('GPP vs Soil Temperature',fontsize=20)
    TSOIL_AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    TSOIL_AX.set_ylim(GPP_limits)
    TSOIL_AX.set_xlim((0,20))
    TSOIL_AX.set_xlabel('Soil Temperature (K)',fontsize=15)
    SMCL_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=SMCL_AX.transAxes, fontsize=15 )
    SMCL_AX.grid(True)
    SMCL_AX.set_title('GPP vs Soil Moisture',fontsize=20)
    SMCL_AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    SMCL_AX.set_ylim(GPP_limits)
    SMCL_AX.set_xlabel('Soil Moisture (kg $m^{2}$)',fontsize=15)
    LAI_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=LAI_AX.transAxes, fontsize=15 )
    LAI_AX.grid(True)
    LAI_AX.set_title('GPP vs LAI',fontsize=20)
    LAI_AX.set_ylabel('GPP ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    LAI_AX.set_ylim(GPP_limits)
    LAI_AX.set_xlim((3.5,4))
    LAI_AX.set_xlabel('Leaf Area Index ($m^{2}$ $m^{-2}$)',fontsize=15)
    
    # FIGURE - Legend and title and save
    handles,labels=TSOIL_AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs+1)
    FIG.suptitle('GPP vs Enviromental Conditions',fontsize=30)
    FIG.savefig(plot_dir+'/GPP_vs_Environment_Scatter_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


if PLOTS[9]=='Y':
    FIG,AXES=plt.subplots(nrows=1,ncols=1,figsize=(20,7))
    FIG.subplots_adjust(right=0.9,top=0.9,bottom=0.15,left=0.05)
    SMCL_AX=AXES.twinx()
    TSOIL_AX=AXES
    LAI_AX=AXES.twinx()
    
    tick_locs=[ np.where( D_pdf[date_limits[0]:date_limits[1]].index == \
                          dt.datetime(2015,mnth,01) )[0] \
                for mnth in range(date_limits[0].month,date_limits[1].month+1) ]
    tick_labels=[ dt.datetime(2015,mnth,01).strftime('%b') \
                  for mnth in range(date_limits[0].month,date_limits[1].month+1) ]
    
    (TSOIL_series[date_limits[0]:date_limits[1]]-273.15).plot(\
            ax=TSOIL_AX,color='brown',label='Soil T',lw=3)
    TSOIL_AX.set_ylabel('T ($^{o}C$)',fontsize=20)
    TSOIL_AX.set_yticks( (0,5,10,15,20) )
    TSOIL_AX.set_zorder(1)
    #TSOIL_AX.yaxis.tick_left()
    #TSOIL_AX.yaxis.set_label_position('left')
    TSOIL_AX.patch.set_visible(False)

    (SMCL_series[date_limits[0]:date_limits[1]]).plot(\
            ax=SMCL_AX,color='dodgerblue',label='Soil Moisture',lw=3)
    SMCL_AX.set_ylabel('SM ($m^{3} m^{-3}$)',fontsize=20)
    SMCL_AX.set_yticks( (0,0.25,0.5,0.75,1.0) )
    SMCL_AX.set_zorder(2)

    LAI_plotseries=LAI_series[date_limits[0]:date_limits[1]]
    LAI_plotseries.plot(ax=LAI_AX,color='#cbf68c',label='LAI',lw=2)
    LAI_AX.set_zorder(0)
    LAI_AX.set_ylim( (3.5,3.9) )
    LAI_AX.spines['right'].set_position(('axes',1.07))
    LAI_AX.set_frame_on(True)
    LAI_AX.patch.set_visible(True)
    for sp in LAI_AX.spines.values():
        sp.set_visible(False)
    LAI_AX.spines['right'].set_visible(True)
    LAI_AX.set_ylabel('LAI ($m^{2} m^{-2}$)',fontsize=20)
    LAI_AX.set_yticks( (3.5,3.6,3.7,3.8,3.9) )
    LAI_AX.fill_between(LAI_plotseries.index, \
                        LAI_plotseries.values, \
                        where=LAI_plotseries>np.zeros_like(LAI_plotseries),\
                        color='#e5fac5')
    
    LAI_handles,LAI_labels=LAI_AX.get_legend_handles_labels()
    TSOIL_handles,TSOIL_labels=TSOIL_AX.get_legend_handles_labels()
    SMCL_handles,SMCL_labels=SMCL_AX.get_legend_handles_labels()
    handles = list(LAI_handles)+list(TSOIL_handles)+list(SMCL_handles)
    labels  = list(LAI_labels) +list(TSOIL_labels) +list(SMCL_labels)
    FIG.legend( handles,labels, loc=8, ncol=nJs+4)
    
    FIG.suptitle('JULES Diagnostics for Crichton Flux Tower',fontsize=30)
    FIG.savefig(plot_dir+'/JULES_DiagnosticData_timeseries_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


if PLOTS[10]=='Y':
    print 'Plotting scatter of NEE versus diagnostics data'
    FIG= plt.figure(figsize=(24,7))
    FIG.subplots_adjust(bottom=0.15,top=0.85,left=0.05,right=0.98)
    TSOIL_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,0) ) )
    SMCL_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,1) ) )
    LAI_AX = FIG.add_subplot(plt.subplot2grid( (1,3), (0,2) ) )

    for iJ in range(nJs):
        J_name=JULES_runs_names[iJ]
        TSOIL_AX.scatter(TSOIL_pdf[date_limits[0]:date_limits[1]][J_name]-273.15,\
                        NEE_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25,label=J_name)
        TSOIL_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%\
                        (TSOIL_pdf[date_limits[0]:date_limits[1]][J_name].corr(   \
                                NEE_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=TSOIL_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        SMCL_AX.scatter(SMCL_pdf[date_limits[0]:date_limits[1]][J_name],\
                        NEE_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25)
        SMCL_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%\
                        (SMCL_pdf[date_limits[0]:date_limits[1]][J_name].corr(    \
                                NEE_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=SMCL_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

        LAI_AX.scatter(LAI_pdf[date_limits[0]:date_limits[1]][J_name],\
                        NEE_pdf[date_limits[0]:date_limits[1]][J_name], \
                         c=J_cols[iJ],lw=0,s=25)
        LAI_AX.text(0.20+(iJ*0.15),0.88,'%4.2f '%(\
                        LAI_pdf[date_limits[0]:date_limits[1]][J_name].corr(\
                            NEE_pdf[date_limits[0]:date_limits[1]][J_name])), \
                     transform=LAI_AX.transAxes, fontsize=15,color=PLOT_COLS[iJ+1] )

    TSOIL_AX.scatter(TSOIL_series-273.15,NEE_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,label='Site',marker='x')
    TSOIL_AX.text(0.05,0.88,'%4.2f '%(TSOIL_series.corr(NEE_pdf['Site'])), \
                  transform=TSOIL_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    SMCL_AX.scatter(SMCL_series,NEE_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,marker='x')
    SMCL_AX.text(0.05,0.88,'%4.2f '%(SMCL_series.corr(NEE_pdf['Site'])), \
                  transform=SMCL_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    LAI_AX.scatter(LAI_series,NEE_pdf['Site'], \
                     c=PLOT_COLS[0],s=25,marker='x')
    LAI_AX.text(0.05,0.88,'%4.2f '%(LAI_series.corr(NEE_pdf['Site'])), \
                  transform=LAI_AX.transAxes, fontsize=15,color=PLOT_COLS[0] )
    
    # NEE text, axis and 1-to-1 line
    NEE_limits=TSOIL_AX.get_ylim()
    TSOIL_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=TSOIL_AX.transAxes, fontsize=15 )
    TSOIL_AX.grid(True)
    TSOIL_AX.set_title('NEE vs Soil Temperature',fontsize=20)
    TSOIL_AX.set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    TSOIL_AX.set_ylim(NEE_limits)
    TSOIL_AX.set_xlim((0,20))
    TSOIL_AX.set_xlabel('Soil Temperature (K)',fontsize=15)
    SMCL_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=SMCL_AX.transAxes, fontsize=15 )
    SMCL_AX.grid(True)
    SMCL_AX.set_title('NEE vs Soil Moisture',fontsize=20)
    SMCL_AX.set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    SMCL_AX.set_ylim(NEE_limits)
    SMCL_AX.set_xlabel('Soil Moisture (kg $m^{2}$)',fontsize=15)
    LAI_AX.text(0.01,0.95,"Pearson's $r^2$: ",transform=LAI_AX.transAxes, fontsize=15 )
    LAI_AX.grid(True)
    LAI_AX.set_title('NEE vs LAI',fontsize=20)
    LAI_AX.set_ylabel('NEE ($\mu$mol $m^{-2} s^{-2}$)',fontsize=15)
    LAI_AX.set_ylim(NEE_limits)
    LAI_AX.set_xlim((3.5,4))
    LAI_AX.set_xlabel('Leaf Area Index ($m^{2}$ $m^{-2}$)',fontsize=15)
    
    # FIGURE - Legend and title and save
    handles,labels=TSOIL_AX.get_legend_handles_labels()
    FIG.legend( handles,labels, loc=8, ncol=nJs+1)
    FIG.suptitle('NEE vs Enviromental Conditions',fontsize=30)
    FIG.savefig(plot_dir+'/NEE_vs_Environment_Scatter_'+TRES+'.png',bbox_inches='tight')
    if iDISPLAY=='Y':
        plt.show()
    else:
        plt.close()


