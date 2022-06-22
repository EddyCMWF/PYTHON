#!/bin/env python2.7
##############################################################################################
#
# Optional inputs:
# -data_dir = Data directory, data for GCMs should be stored in a subdirectory with the 
#                              name of the GCM as in data_info.GCMs()
#               Default='/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/'
# -out_dir  = output directory to store the plots.
#               Default=DATA_DIR+'plots/'

import pdb
import numpy as np
import netCDF4 as nc
import sys,os 
#import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from scipy import stats

import data_info

def optional_argparse(arg,default):
    if arg in sys.argv:
        temp_loc=sys.argv.index(arg)
        temp_arg=sys.argv.pop(temp_loc)
        value=sys.argv.pop(temp_loc)
    else:
        value=default
    return value

def JULES_Q10(q10,T_SOIL,T0=273.15):
    return q10**(0.1*(T0/T_SOIL)*(T_SOIL-T0))

# Constants:
secs_to_year=3600.*24.*360. 
secs_to_month=3600.*24.*30. 
min_corre=0.95

TEST_MODE=True
DATA_DIR=optional_argparse('-data__dir', \
                           #'/work/scratch/ecomynplatt/CLIFFTOP/BASELINE_CONFIG/')
                           '/prj/CLIFFTOP/ECP_output/METHANE_FEEDBACK/')
while not os.path.isdir(DATA_DIR):
    DATA_DIR=raw_input('data_dir "'+DATA_DIR+'" not found, please input existing data_dir: ')
if not DATA_DIR[-1]=='/':
    DATA_DIR=DATA_DIR+'/'
print('data_dir: ',DATA_DIR)

OUT_DIR=optional_argparse('-out_dir','./')
if not OUT_DIR[-1]=='/':
    OUT_DIR=OUT_DIR+'/'
if not os.path.isdir(OUT_DIR):
    os.system('mkdir '+OUT_DIR)

# Soil Options from JULES
# kappa resp factors
Kappa_rothC = [ 3.22e-7,9.65e-9,2.12e-8,6.43e-10 ]
# and weights:
Kappa_weights = [ kappa/np.sum(Kappa_rothC) for kappa in Kappa_rothC ]
npools=len(Kappa_rothC)
# decay factor of resp with depth
tau_resp = 2.0
# soil layer thicknesses
dz_soil= np.array([0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293,
          0.21517585,0.23784142,0.25980762,0.28117066,0.30200527,
          0.32237098,0.34231625,0.36188121])
nz = len(dz_soil)
# Soil layer depths:
ztot = np.zeros_like(dz_soil)
ztot[0]=dz_soil[0]*0.5
for iz in range(1,nz):
    ztot[iz] = ztot[iz-1] + 0.5*(dz_soil[iz-1]+dz_soil[iz])

# CH4 depth factor:
CH4_depth_factor = np.exp(-tau_resp*ztot)


scenario=optional_argparse('-scenario','1p5deg')
GCMs=optional_argparse('-GCMs','ALL')
if GCMs.upper()=='ALL':
    GCMs=data_info.GCMs()
else:
    GCMs=GCMs.split(',')

YEAR=int(optional_argparse('-year','2000'))
GLOBAL_TOTAL_ref=float(optional_argparse('-global_total_ref','180'))
#time_res = 'Annual'
time_res = 'Monthly'

RUNID=optional_argparse('-runid','vn4.8_imogen')

AREA_file=optional_argparse('-area_file',\
                         '/prj/CLIFFTOP/COMMON_DATA/ANCILS/Area_in_iris_format.nc')
AREAinf=nc.Dataset(AREA_file,'r')
AREA=AREAinf.variables['area'][:].squeeze()

if TEST_MODE:
    GCMs=['CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0','CEN_MOHC_MOD_HadGEM2-ES']# ,'CEN_NOAA-GFDL_MOD_GFDL-ESM2G']

# Open File
outf=open(OUT_DIR+'fch4_TuretskyQ10_CS.txt','w')

Q10_ensemble=np.arange(1.0,15.1,0.1)
nQ10s=len(Q10_ensemble)
GlobTot_ensemble=np.arange(50.,250.1,10.)
nGlobTots=len(GlobTot_ensemble)
glob_tot_index = np.argmin(np.abs(GlobTot_ensemble-GLOBAL_TOTAL_ref))

#T_sensitvity_ensemble = [ [] for i in range(nQ10s) ]
T0=273.15
# First read in the data, GCM by GCM, then Year by Year (then month by month is required)

#Turetsky_Q10s = [ 2.6,   2.0,        1.7]
#Turetsky_As   = [ 17.4,  53.5,       50.2]
#Turetsky_type = ['Bog', 'Rich Fen', 'Poor Fen']
Turetsky_Q10 = 2.0
Turetsky_A   = 53.5
mgperday_to_kgperyear=(1e-6)*(360.)

scheme_names=['CS','NPP','RESP_S']
outf.write( 35*' '+'%17s%8s'%('Q10_Cs (r, rmse)','k_CS ') +'\n' )
Optimised_Q10s={}

for gcm in GCMs:
    Optimised_Q10s[gcm]={}
    print(gcm)
    outf.write( '%35s '%(gcm[0-min(34,len(gcm)):]) )
    FILE_C = DATA_DIR+gcm+'/'+RUNID+'_'+gcm+'_'+scenario+'.'+time_res+'_carbon.'+str(YEAR)+'.nc'
    FILE_HT = DATA_DIR+gcm+'/'+RUNID+'_'+gcm+'_'+scenario+'.'+time_res+'_h2o_t.'+str(YEAR)+'.nc'
   
    # Read Data:
    inf_c=nc.Dataset(FILE_C,'r')
    inf_ht=nc.Dataset(FILE_HT,'r')
    ntime=len(inf_c.variables['time'][:])
    lats=inf_c.variables['latitude'][:].squeeze()
    nland=len(lats)
    T_SOIL = inf_ht.variables['t_soil'][:].reshape(ntime,nz,nland)
    FWETL  = inf_ht.variables['fwetl'][:].reshape(ntime,nland)
    CS     = inf_c.variables['cs'][:].reshape(ntime,npools,nz,nland)
    FWETL_AREA =FWETL*AREA 
    inf_c.close()
    inf_ht.close()

    # Turetsky T is top 10 cm, i.e. to 2 layers
    T_turetsky = np.mean(T_SOIL[:,:2,:],axis=1)-273.15
    #T_sensitivity_turetsky = [ q10**(0.1*T_turetsky) for q10 in Turetsky_Q10s ]
    T_sensitivity_turetsky = (Turetsky_Q10**(0.1*T_turetsky))
    # FCH4 in kg per m^2 per year:
    FCH4_turetsky = Turetsky_A * T_sensitivity_turetsky * mgperday_to_kgperyear
    # Total kg per year for gridcell
    CH4_turetsky = FCH4_turetsky*FWETL_AREA   

    # Cs T sensitivity
    T_sensitivity_CS_ensemble= np.zeros([nQ10s,ntime,nland])

    #pdb.set_trace()
    # Calculate the T sensitive component of the JULES equation:
    for iQ10 in range(nQ10s):
        q10 = Q10_ensemble[iQ10]
        for ilayer in range(nz):
            depth_factor = CH4_depth_factor[ilayer]
            jules_q10 = JULES_Q10(q10,T_SOIL[:,ilayer,:])
            for ipool in range(npools):
                kappa_weight=Kappa_weights[ipool]
                cs=CS[:,ipool,ilayer,:]
                T_sensitivity_CS_ensemble[iQ10,:]+=(kappa_weight*cs*jules_q10*depth_factor)
    
    T_sensitivity_CS_correlation = \
            np.array([ stats.pearsonr((T_sensitivity_turetsky*FWETL_AREA).flatten(),
                                      (T_sensitivity_CS_ensemble[iQ10,:]*FWETL_AREA).flatten())[0]
                                             for iQ10 in range(nQ10s) ]).flatten() 

    CS_FACTORS     = np.zeros([nGlobTots,nQ10s])

    RMSE_CS = np.zeros([nGlobTots,nQ10s])

    RMSE_CS_scaled = np.zeros([nGlobTots,nQ10s])

    for iGlobTot in range(nGlobTots):
        # convert year tot to kg of CH4 per year
        YEAR_tot=GlobTot_ensemble[iGlobTot]*1e9

        Turetsky_Factor = YEAR_tot / np.sum( np.mean(CH4_turetsky,axis=0) )
        CH4_turetsky_scaled = CH4_turetsky*Turetsky_Factor
        FCH4_turetsky_scaled = FCH4_turetsky*Turetsky_Factor

        # Calculate the global k values from JULES equation based on the global total:
        # JULES works in kg C per m^2 per s
        # Multiply T sensitivty by wetland areai and convert from C to to CH4:
        CH4_CS_temporary      = T_sensitivity_CS_ensemble*FWETL_AREA*(16.04/12.01)
        
        # Convert to annual total by converting to per year then mean over time axis:
        CH4_CS_temporary      = np.mean(CH4_CS_temporary*secs_to_year,axis=1)
        
        # Now calculate the k by desired global total by sum of T_sensitve component:
        CS_FACTORS[iGlobTot,:]     = (YEAR_tot) / np.sum(CH4_CS_temporary,axis=1)
        
        # FCH4 in kgCH4 per m^2 per year
        FCH4_CS     = np.array([T_sensitivity_CS_ensemble[iQ10,:]    *CS_FACTORS[iGlobTot,iQ10]   \
                                                *secs_to_year*(16.04/12.01)
                                for iQ10 in range(nQ10s) ] )
        
        #CH4_CS = FCH4_CS * FWETL

        RMSE_CS[iGlobTot,:]     = np.array( [ np.mean( (((FCH4_turetsky-FCH4_CS[iQ10,:])*FWETL)**2.)**0.5 )
                                             for iQ10 in range(nQ10s) ] )
    
        RMSE_CS_scaled[iGlobTot,:] = np.array( [ np.mean( (((FCH4_turetsky_scaled-FCH4_CS[iQ10,:])*FWETL)**2.)**0.5 )
                                              for iQ10 in range(nQ10s) ] )
    

    # Optimised Q10 values minimum RMSE where correlation is greater than 90% of highest correlation
    CS_r_index= np.where(T_sensitivity_CS_correlation>(min_corre*T_sensitivity_CS_correlation.max()))[0]
    CS_Q10_index = CS_r_index[ np.argmin(RMSE_CS_scaled[glob_tot_index,CS_r_index]) ]
    CS_Q10=Q10_ensemble[CS_Q10_index]
    CS_Q10_rval=T_sensitivity_CS_correlation[CS_Q10_index] 
    CS_Q10_RMSE=RMSE_CS_scaled[glob_tot_index,CS_Q10_index]
    CS_kFACTOR=CS_FACTORS[glob_tot_index,CS_Q10_index]
    
    Optimised_Q10s[gcm] = { 'CS':    { 'Q10':CS_Q10,'r':CS_Q10_rval,'rmse':CS_Q10_RMSE,'k':CS_kFACTOR},
                            }

    outf.write('%4.1f (%4.2f, %6.2e), %9.3e'%tuple([ Optimised_Q10s[gcm]['CS'][key] 
                                                     for key in ['Q10','r','rmse','k']]) )
    outf.write('\n')
    
    # Plot Section:
    # Line plots
    fig=plt.figure(figsize=(8,5))
    ax=fig.add_subplot(111)
    ax.plot(Q10_ensemble,T_sensitivity_CS_correlation,label='Cs-r',c='b',lw=3)
    ax.plot([Optimised_Q10s[gcm]['CS']['Q10'],Optimised_Q10s[gcm]['CS']['Q10']],[0,1.0],ls='--',c='b',lw=2)
    ax.set_ylabel("Pearson's r")
    ax.set_ylim([0,ax.get_ylim()[-1]])
    ax.set_xlabel("Q10")
    ax.grid(True)
    ax2=ax.twinx()
    ax2.plot(Q10_ensemble,RMSE_CS_scaled[glob_tot_index,:],label='Cs-RMSE',ls='-.',c='b',lw=3)
    ax2.set_ylabel("RMSE (kgCH4 m$^{-2}$ year$^{-1}$)")
    handles1,labels1 = ax.get_legend_handles_labels()
    handles2,labels2 = ax2.get_legend_handles_labels()
    labels=labels1+labels2
    handles=handles1+handles2
    fig.legend(handles,labels,loc=8,ncol=2)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(OUT_DIR+'plots/'+gcm+'_pearsonsr_rmse180Tg_CS.png')
    plt.close()
    
    # Unscaled Turetsky 
    fig=plt.figure(figsize=(8,5))
    ax=fig.add_subplot(111)
    ax.plot(Q10_ensemble,T_sensitivity_CS_correlation,label='Cs-r',c='b',lw=3)
    ax.plot([Optimised_Q10s[gcm]['CS']['Q10'],Optimised_Q10s[gcm]['CS']['Q10']],[0,1.0],ls='--',c='b',lw=2)
    ax.set_ylabel("Pearson's r")
    ax.set_ylim([0,ax.get_ylim()[-1]])
    ax.set_xlabel("Q10")
    ax.grid(True)
    ax2=ax.twinx()
    ax2.plot(Q10_ensemble,RMSE_CS[glob_tot_index,:],label='Cs-RMSE',ls=':',c='b',lw=3)
    ax2.set_ylabel("RMSE (kgCH4 m$^{-2}$ year$^{-1}$)")
    handles1,labels1 = ax.get_legend_handles_labels()
    handles2,labels2 = ax2.get_legend_handles_labels()
    labels=labels1+labels2
    handles=handles1+handles2
    fig.legend(handles,labels, loc=8,ncol=2)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(OUT_DIR+'plots/'+gcm+'_pearsonsr_rmse180Tg_CS_unscaled.png')
    plt.close()

    
    # Contour plots - Scaled
    fig2,axes=plt.subplots(ncols=1,nrows=1,figsize=(12,10))
    fig2.subplots_adjust(wspace=0.5)
    #Plot CS
    ax=axes
    im=ax.contourf(Q10_ensemble,GlobTot_ensemble,RMSE_CS_scaled,cmap='gray')
    ax.set_ylabel('Global Total (Tg CH4 per year)')
    ax.set_xlabel('Q10 (CS)')
    cax=fig2.colorbar(im,ax=ax,orientation='horizontal')  #pad=0.1)
    cax.set_label('RMSE (Gg per year)')
    ax2=ax.twinx()
    ax2.plot(Q10_ensemble,T_sensitivity_CS_correlation,lw=3,c='r')
    ax2.plot([Optimised_Q10s[gcm]['CS']['Q10'],Optimised_Q10s[gcm]['CS']['Q10']],[0,1.0],ls=':',c='r',lw=3)
    ax2.set_xlim(Q10_ensemble[0],Q10_ensemble[-1])
    ax2.set_ylim(0,1.0)
    ax2.set_ylabel("Pearson's R")
    ax.set_title("Soil Carbon")
    fig2.tight_layout(w_pad=0.3,pad=4)
    fig2.savefig(OUT_DIR+'plots/'+gcm+'_RMSEcontour_CS.png')
    plt.close()

    # Plot Zonal Time-series
    if time_res=='Monthly':
        Boreal_mask = lats>50
        MidLat_mask = (lats>20)&(lats<50)
        Tropic_mask = (lats>-20)&(lats<20)
        plot_months=np.arange(1,13)

        FCH4_CS     =  T_sensitivity_CS_ensemble[CS_Q10_index,:] \
                     * CS_FACTORS[glob_tot_index,CS_Q10_index]   \
                        * secs_to_month*(16.04/12.01)
        CH4_CS=FCH4_CS*FWETL_AREA
        CH4_CS_timeseries = { 'Global': np.sum( CH4_CS, axis=1)*1e-9, 
                              'Boreal': np.sum( CH4_CS[:,Boreal_mask], axis=1)*1e-9,
                              'MidLat': np.sum( CH4_CS[:,MidLat_mask], axis=1)*1e-9,
                              'Tropic': np.sum( CH4_CS[:,Tropic_mask], axis=1)*1e-9,
                              }
        
        Turetsky_Factor = GLOBAL_TOTAL_ref / np.sum( np.mean(CH4_turetsky,axis=0) )
        CH4_turetsky_scaled = CH4_turetsky*Turetsky_Factor/12.
        CH4_turetsky_timeseries = { 'Global': np.sum( CH4_turetsky_scaled, axis=1), 
                                    'Boreal': np.sum( CH4_turetsky_scaled[:,Boreal_mask], axis=1),
                                    'MidLat': np.sum( CH4_turetsky_scaled[:,MidLat_mask], axis=1),
                                    'Tropic': np.sum( CH4_turetsky_scaled[:,Tropic_mask], axis=1),
                                   }
       
        fig=plt.figure(figsize=[12,18])
        ax=fig.add_subplot(411)
        ax.plot(plot_months,CH4_CS_timeseries['Global'],label='CH4 - CS',c='b')
        ax.plot(plot_months,CH4_turetsky_timeseries['Global'],label='CH4 - turetsky',c='k')
        ax.set_ylabel('Global Emissions (Tg CH4 per month)')
        ax.legend()
        ax=fig.add_subplot(412)
        ax.plot(plot_months,CH4_CS_timeseries['Boreal'],label='CH4 - CS',c='b')
        ax.plot(plot_months,CH4_turetsky_timeseries['Boreal'],label='CH4 - turetsky',c='k')
        ax.set_ylabel('Boreal Emissions (Tg CH4 per month)')
        ax=fig.add_subplot(413)
        ax.plot(plot_months,CH4_CS_timeseries['MidLat'],label='CH4 - CS',c='b')
        ax.plot(plot_months,CH4_turetsky_timeseries['MidLat'],label='CH4 - turetsky',c='k')
        ax.set_ylabel('MidLat Emissions (Tg CH4 per month)')
        ax=fig.add_subplot(414)
        ax.plot(plot_months,CH4_CS_timeseries['Tropic'],label='CH4 - CS',c='b')
        ax.plot(plot_months,CH4_turetsky_timeseries['Tropic'],label='CH4 - turetsky',c='k')
        ax.set_ylabel('Tropic Emissions (Tg CH4 per month)')
        ax.set_xlabel('Month')
        fig.suptitle('CH4 time series') 
        fig.savefig(OUT_DIR+'plots/'+gcm+'_ZonalTimeSeries_CS.png')


        plt.close()

outf.close()




