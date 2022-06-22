#!/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import parabolic, profile, ocean_co2, data_info

PLOT_TAG='Trpofiles_gcmECS'
out_dir='/users/eow/edwcom/CLIFFTOP/plots/toy_jules/'+PLOT_TAG+'/'
os.system('mkdir -p '+out_dir)

GCMs = data_info.GCMs()
nGCMs = len(GCMs)
cmip5_runs = data_info.cmip5_runs()

kappa = data_info.kappa(cmip5_runs)
lambda_o = data_info.lambda_o(cmip5_runs)
lambda_l = data_info.lambda_l(cmip5_runs)
f_ocean = data_info.ocean_frac(cmip5_runs)
nu = data_info.nu(cmip5_runs)
 
SCENARIOS={'1p5deg':{'dt_limit':1.5,'mu_zero':0.04,'mu_one':0.0,
                    'color':'#6495ED','label':'1.5$^o$C','ls':'-'},
        '1p81p5deg':{'dt_limit':1.5,'mu_zero':-0.024,'mu_one':0.00154,
                    'color':'#EDBC64','label':'1.5$^o$C (overshoot)','ls':'--'},
           #'1p81p5deg_CW':{'dt_limit':1.5,'mu_zero':-0.024,'mu_one':0.001540},
           '2deg':{'dt_limit':2.0,'mu_zero':0.04,'mu_one':0.0,
               'color':'#ED6495','label':'2.0$^o$C','ls':':'}, 
           '3deg':{'dt_limit':3.0,'mu_zero':0.065,'mu_one':0.0,
               'color':'#6F36C8','label':'3.0$^o$C','ls':'-.'} }
scen_plot_order = ['3deg','2deg','1p81p5deg','1p5deg']
l_outTprof=True   # Output temperature profile?
version='vn3p0'   # For saving output driving files
l_saveimg=True

#Tprofiles=np.loadtxt('CMIP5_tprofile_data.dat')
#np.save('CMIP5_tprofile_data.npy',Tprofiles)
Tprofiles=np.load('CMIP5_tprofile_data.npy')
Tprofiles=np.ma.masked_equal(Tprofiles,-9999.0)
# Also normalise the curves so that they end up at the current temperature and gradient estimate
beta_T=0.025               # K/yr
dTg_dt=beta_T

dQnonCO2_dt = 0.0122     # Current rate of change in nonCO2 radiative forcing with time (Wm-2/yr)
dQCO2_dt = 0.03          # Current rate of change in CO2 radiative forcing with time (Wm-2/yr)
dQ_dt = dQnonCO2_dt+dQCO2_dt 

k = ( ( (dQ_dt/dTg_dt)*(f_ocean+nu-(nu*f_ocean)) ) - (nu*lambda_l)
      / ( (f_ocean*nu)-nu )  )

# Current rate of change W m-2/yr
zeta = (f_ocean+nu-(nu*f_ocean)) / (lambda_l-(k*f_ocean)+(k*f_ocean*nu))
beta_Qco2_gcm = zeta*dQCO2_dt
beta_Q_gcm = zeta*dQ_dt

#dt_now=0.89
#yr_now=2015
dt_now=0.89
yr_now=2015
yr_beta=1995
end_year=2200
start_year=1850
n_yr=end_year-start_year+1

for scenario in SCENARIOS.keys():
    SCENARIOS[scenario]['Tprofile']=np.array([profile.profile(beta_Q_gcm[igcm], dt_now, 
    #SCENARIOS[scenario]['Tprofile']=np.array([profile.profile(0.025+(beta_Q_gcm[igcm]*0.), dt_now, 
                                                       SCENARIOS[scenario]['dt_limit'], 
                                                       SCENARIOS[scenario]['mu_zero'], 
                                                       SCENARIOS[scenario]['mu_one'], 
                                                       yr_now=yr_now, yr_beta=yr_beta)[:n_yr]
                                               for igcm in range(nGCMs) ])

yr_plot=np.arange(start_year,end_year+1)
now_index=np.where(yr_plot==yr_now)[0][0]
print(yr_plot[:251].shape)
print(Tprofiles.shape)

FIG,AXES=plt.subplots(figsize=[14,6],ncols=1,nrows=1)
FONTSIZE=35
# Plot total RF in top plot
ax=AXES
ax.plot(yr_plot[:251],Tprofiles.transpose(1,0),c='#e8e8e8',lw=2)
for scenario in scen_plot_order: #SCENARIOS.keys():
    leg_label= SCENARIOS[scenario]['label']                        \
             + r', $\mu_{0}$='+str(SCENARIOS[scenario]['mu_zero']) \
             + r', $\mu_{1}$='+str(SCENARIOS[scenario]['mu_one']) 
    
    ax.plot(yr_plot,SCENARIOS[scenario]['Tprofile'],
            #label=SCENARIOS[scenario]['label'],
            label=leg_label,
            color=SCENARIOS[scenario]['color'],lw=3.)
    # Output Tprofile if required
    if l_outTprof:
        out_Tprof_filename=out_dir+scenario+'_global_temp_anomaly_'+version+'.dat'
        outf=open(out_Tprof_filename,'w')
        for year,Temp in zip(yr_plot,SCENARIOS[scenario]['Tprofile']):
            line='%4i  %8.4f\n'%(year,Temp)
            outf.write(line)
        outf.close()
ax.set_ylabel('Global Temperature Change (K)',fontsize=FONTSIZE/2.)
ax.grid(True)
ax.legend(loc=2)
ax.set_ylim([-0.25,3.5])
ax.set_yticks(np.arange(-0.25,3.51,0.25))
ax.set_xlim([1850,2100])
ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=2.)
#ax.plot(ax.get_xlim(),[1.8,1.8],c='k',lw=1.5)
#ax.plot(ax.get_xlim(),[1.75,1.75],c='k',lw=.5,ls=':')
#FIG.suptitle(PLOT_TAG+', Imogen vs '+PLOT_TAG+' '+method,fontsize=FONTSIZE)
if l_saveimg:
    FIG.savefig(out_dir+'Tprofiles'+version+'.png',bbox_inches='tight')
    FIG.savefig(out_dir+'Tprofiles'+version+'.eps',bbox_inches='tight')
#plt.show()
plt.close()

FIG,AXES=plt.subplots(figsize=[14,6],ncols=1,nrows=1)
FONTSIZE=35
# Plot total RF in top plot
ax=AXES
ax.plot(yr_plot[:251],Tprofiles.transpose(1,0),c='#e8e8e8',lw=2)
for scenario in scen_plot_order: #SCENARIOS.keys():
    leg_label= SCENARIOS[scenario]['label']                        \
             + r', $\mu_{0}$='+str(SCENARIOS[scenario]['mu_zero']) \
             + r', $\mu_{1}$='+str(SCENARIOS[scenario]['mu_one']) 
    ax.plot(yr_plot,SCENARIOS[scenario]['Tprofile'],
            #label=SCENARIOS[scenario]['label'],
            label=leg_label,
            ls=SCENARIOS[scenario]['ls'],lw=4.,c='k')
    # Output Tprofile if required
    if l_outTprof:
        out_Tprof_filename=out_dir+scenario+'_global_temp_anomaly_'+version+'.dat'
        print(out_Tprof_filename)
        outf=open(out_Tprof_filename,'w')
        for year,Temp in zip(yr_plot,SCENARIOS[scenario]['Tprofile']):
            line='%4i  %8.4f\n'%(year,Temp)
            outf.write(line)
        outf.close()
ax.set_ylabel('Global Temperature Change (K)',fontsize=FONTSIZE/2.)
ax.grid(True)
ax.legend(loc=2)
ax.set_ylim([-0.25,3.5])
ax.set_yticks(np.arange(-0.25,3.51,0.25))
ax.set_xlim([1850,2100])
ax.plot([yr_now,yr_now],ax.get_ylim(),c='k',lw=2.)
#ax.plot(ax.get_xlim(),[1.8,1.8],c='k',lw=1.5)
#ax.plot(ax.get_xlim(),[1.75,1.75],c='k',lw=.5,ls=':')
#FIG.suptitle(PLOT_TAG+', Imogen vs '+PLOT_TAG+' '+method,fontsize=FONTSIZE)
if l_saveimg:
    FIG.savefig(out_dir+'Tprofiles_greyscale.png',bbox_inches='tight')
    FIG.savefig(out_dir+'Tprofiles_greyscale.eps',bbox_inches='tight')
plt.show()
#plt.close()




print(np.max(SCENARIOS['1p81p5deg']['Tprofile']),
      (SCENARIOS['1p5deg']['Tprofile']-SCENARIOS['1p81p5deg']['Tprofile'])[250])
print(SCENARIOS['2deg']['Tprofile'][250])
print(SCENARIOS['3deg']['Tprofile'][250])
#print(np.max(SCENARIOS['1p81p5deg_CW']['Tprofile']),
#      (SCENARIOS['1p5deg']['Tprofile']-SCENARIOS['1p81p5deg_CW']['Tprofile'])[250])


    
