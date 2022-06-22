#!/bin/env python3.5

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

plot_dir='/prj/GREENHOUSE/GREENHOUSE_sites/plots/Brattleby/crop_tuning_Season1/'
SITE_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
JULES_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/output/Brattleby/Calibration2/'
Site='Brattleby'
tag='_crop_Season'
iSEAS = 1

Season_Limits=('2011-10-01','2012-09-01')
Corre_Limits=('2012-06-01','2012-08-20')

Lc_H2O = 2.501e6
kgC_to_umolsCO2_factor = (1./12.)* 1e9

wk_fit_vars=['lai','canht','shf','lhf','gpp','nee','ter']
JU_fit_vars=['lai','canht','ftl_gb','fqw_gb','gpp_gb','npp_gb','resp_s_gb','resp_p_gb']
SI1_fit_vars=['SHF','LHF','GPP','NEE','TER']
SI2_fit_vars=['H','LE','GPP','NEE','TER']
pft_index=5

SI_veg_vars=['lai','canht','lai_sd']
##################
#Read in SITE VEG data:
SITE_veg_file=SITE_data_dir+'Brattleby_LAI_Canht_data.nc'
Sinf=nc.Dataset(SITE_veg_file,'r')
SITE_DICT={ var:Sinf.variables[var][:].squeeze() \
            for var in SI_veg_vars }
SITE_time=nc.num2date(Sinf.variables['date'][:],        \
                      units=Sinf.variables['date'].units)
SITE_VEG_panda=pd.DataFrame(SITE_DICT,index=SITE_time)
Sinf.close()

##################
#Read in SITE flux data:
SITE_flux_file1=SITE_data_dir+'ConCrop_FullFlux_data_albmar.nc'
Sinf=nc.Dataset(SITE_flux_file1,'r')
SITE_time=nc.num2date( Sinf.variables['date'][:], \
                       units=Sinf.variables['date'].units )
SITE_dict={ wkvar:Sinf.variables[invar][:].squeeze()  \
        for wkvar,invar in zip(wk_fit_vars[2:],SI1_fit_vars) }
Sinf.close()
SITE_panda=pd.DataFrame(SITE_dict,index=SITE_time).resample('D').mean()

###################################
# Concat fluxes and vegetation
SITE_panda=pd.concat([SITE_panda,SITE_VEG_panda],axis=1)

SITE_index=((SITE_panda.index>Season_Limits[0])&(SITE_panda.index<Season_Limits[1]))

SITE_panda=SITE_panda[SITE_index]
SITE_panda['nee']*=-1.0

##################################################
#Read in the JULES data
MU_opts=[ '0.015', '0.02', '0.025'  ]#, '0.015', '0.04' ]
#CRIT_opts= ['1.45', '1.48', '1.50', '1.52', '1.55']
CRIT_opts= ['1.48', '1.49', '1.50', '1.51', '1.52' ]
SEN_opts=[ '1.00', '1.10', '1.15', '1.20'] #, '1.00', '1.15', '1.20' ]

nMU=len(MU_opts)
nCRIT=len(CRIT_opts)
nSEN=len(SEN_opts)

template_array=np.zeros([nMU,nCRIT,nSEN])
stats=['mean','max','correlation','rmse','rmse_rel','stddev']

STAT_dict = { var: { stat:template_array.copy() for stat in stats } \
                 for var in wk_fit_vars }

fit_params=['mu','dvi_crit','sen_dvi']
PARAM_dict= { param:template_array.copy() for param in fit_params } 

#######################################################################
for imu in range(nMU):
  for isen in range(nSEN):
    for icrit in range(nCRIT):
        param_vals=[MU_opts[imu],CRIT_opts[icrit],SEN_opts[isen]]
        for param,val in zip(fit_params,param_vals):
            PARAM_dict[param][imu,icrit,isen] =float(val)

        Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                MU_opts[imu]+'_'+   \
                CRIT_opts[icrit]+'_'+ \
                SEN_opts[isen]+ \
                '.day.nc'

        print(Jinfile)
        Jinf=nc.Dataset(Jinfile,'r')
        temp_time=nc.num2date( Jinf.variables['time'][:],         \
                               units=Jinf.variables['time'].units )
        temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                for wkvar,invar in zip(wk_fit_vars[2:5],JU_fit_vars[2:5]) }
        temp_dict['lhf']*=Lc_H2O 
        temp_dict['gpp']*=kgC_to_umolsCO2_factor
        temp_dict['nee']=(Jinf.variables['npp_gb'][:].squeeze() - \
                      Jinf.variables['resp_s_gb'][:].squeeze())  \
                     * kgC_to_umolsCO2_factor 
        temp_dict['ter']=(Jinf.variables['resp_s_gb'][:].squeeze() + \
                      Jinf.variables['resp_p_gb'][:].squeeze())  \
                     * kgC_to_umolsCO2_factor 
        temp_dict['lai']=Jinf.variables['lai'][:,pft_index].squeeze()
        temp_dict['canht']=Jinf.variables['canht'][:,pft_index].squeeze()
        Jinf.close()
        
        JU_panda=pd.DataFrame(temp_dict,index=temp_time)
        JU_index=((JU_panda.index>Corre_Limits[0])&(JU_panda.index<Corre_Limits[1]))
        JU_panda=JU_panda[JU_index]

        for var in wk_fit_vars:
            STAT_dict[var]['mean'][imu,icrit,isen]=(JU_panda)[var].mean()
            STAT_dict[var]['max'][imu,icrit,isen]=(JU_panda)[var].max()
            STAT_dict[var]['rmse'][imu,icrit,isen]=(JU_panda-SITE_panda)[var].abs().mean()
            STAT_dict[var]['stddev'][imu,icrit,isen]=(JU_panda-SITE_panda)[var].abs().std()
            STAT_dict[var]['correlation'][imu,icrit,isen]=JU_panda[var].corr(SITE_panda[var])
            STAT_dict[var]['rmse_rel'][imu,icrit,isen] = \
                    ((JU_panda-SITE_panda)/(JU_panda+SITE_panda))[var].abs().mean()*2.

#######################################################################
for stat in ('correlation','rmse'):
    print(stat)
    for var in ['lai','canht',]: #'shf','lhf','nee']:
        print(var,STAT_dict[var][stat])

lai_corr_min=0.9
lai_rmse_max=0.2
canht_corr_min=0.5
canht_rmse_max=0.1  
shf_corr_min=0.70
shf_rmse_max=20.  

#lhf_corr_min=0.37
#lhf_rmse_max=35.0
#nee_corr_min=0.73
#nee_rmse_max=9.7

index=np.where(  (STAT_dict['lai']['correlation']>lai_corr_min)     \
#              & (STAT_dict['canht']['correlation']>canht_corr_min) \
               & (STAT_dict['shf']['correlation']>shf_corr_min)     \
#              & (STAT_dict['lhf']['correlation']>lhf_corr_min)     \
#              & (STAT_dict['nee']['correlation']>nee_corr_min)     \
               & (STAT_dict['lai']['rmse']<lai_rmse_max)            \
               & (STAT_dict['canht']['rmse']<canht_rmse_max)        \
               & (STAT_dict['shf']['rmse']<shf_rmse_max)            \
#              & (STAT_dict['lhf']['rmse']<lhf_rmse_max)            \
#              & (STAT_dict['nee']['rmse']<nee_rmse_max)            \
               )

nCONFIGS=len(index[0])
print(nCONFIGS)

out_fname=plot_dir+'/Calibration2_fit.txt'
outf=open(out_fname,'w')
outf.write('Fit criteria: \n')
outf.write('LAI minimum correlation: %6.3f\n'%lai_corr_min)
outf.write('LAI maximum rmse: %6.3f\n'%lai_rmse_max)
outf.write('Canht minimum correlation: %6.3f\n'%canht_corr_min)
outf.write('Canht maximum rmse: %6.3f\n'%canht_rmse_max)
#outf.write('SHF minimum correlation: %6.3f\n'%shf_corr_min)
#outf.write('SHF maximum rmse: %6.3f\n'%shf_rmse_max)
#outf.write('LHF minimum correlation: %6.3f\n'%lhf_corr_min)
#outf.write('LHF maximum rmse: %6.3f\n'%lhf_rmse_max)

outf.write('\n\n')
outf.write('Satisfying Configurations:\n')
for param in fit_params: 
    print(param)
    outf.write('%15a'%(param)+':'+nCONFIGS*' %10.3f ' % tuple(PARAM_dict[param][index]) +'\n')

for var in ['lai','canht','shf','lhf','nee']:
    for stat in ['correlation','rmse']:
        string=var+','+stat[:4]
        outf.write('%15a'%(string)+':'+nCONFIGS*' %10.4f '% tuple(STAT_dict[var][stat][index]) +'\n')

outf.close()


# Plot selected configs
configs_temp = [ PARAM_dict[param][index] for param in fit_params ]

configs = [ str(configs_temp[0][iconf])+'_'+  \
            str('%4.2f'%configs_temp[1][iconf])+'_'+  \
            str('%4.2f'%configs_temp[2][iconf])       \
              for iconf in range(nCONFIGS) ]
configs=[ config.replace(' ','') for config in configs ]

configs += ['0.05_2.00_1.50']
nCONFIGS =len(configs)

JU_panda_list=[]
for config in configs:
    Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
            config + \
            '.day.nc'
    Jinf=nc.Dataset(Jinfile,'r')
    temp_time=nc.num2date( Jinf.variables['time'][:],         \
                            units=Jinf.variables['time'].units )
        
    temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                for wkvar,invar in zip(wk_fit_vars[2:5],JU_fit_vars[2:5]) }
    temp_dict['lhf']*=Lc_H2O 
    temp_dict['gpp']*=kgC_to_umolsCO2_factor
    temp_dict['nee']=(Jinf.variables['npp_gb'][:].squeeze() - \
                      Jinf.variables['resp_s_gb'][:].squeeze())  \
                       * kgC_to_umolsCO2_factor 
    temp_dict['ter']=(Jinf.variables['resp_s_gb'][:].squeeze() + \
                      Jinf.variables['resp_p_gb'][:].squeeze())  \
                         * kgC_to_umolsCO2_factor 
    temp_dict['lai']=Jinf.variables['lai'][:,pft_index].squeeze()
    temp_dict['canht']=Jinf.variables['canht'][:,pft_index].squeeze()

    temp_dict['fsmc']=Jinf.variables['fsmc'][:,pft_index].squeeze() 

    Jinf.close()

    temp_panda=pd.DataFrame(temp_dict,index=temp_time)
    temp_index=((temp_panda.index>Season_Limits[0])&(temp_panda.index<Season_Limits[1]))
    temp_panda=temp_panda[temp_index]
    
    JU_panda_list.append(temp_panda.copy()) 


full_config_names= configs 

#long_color_list=['b','g','k','y','c','orange','aqua','khaki','gold','grey']
#plot_colors = long_color_list[:nCONFIGS] 

plot_colors = [ cm.get_cmap('Spectral')(fraction) for fraction in np.arange(0,1.00000001,1./nCONFIGS) ]

plot_var_list=['shf','lhf','nee','lai','canht','gpp']
units={'shf':'W m$^{-2}$','lhf':'W m$^{-2}$','nee':'umolCO2 m$^{-2}$ $s^{-1}$',\
        'lai':'m$^2$ m$^{-2}$','canht':'m','gpp':'umolCO2 m$^{-2}$ $s^{-1}$'}

fsmc_index=-1
print(fsmc_index)
FIG,AXES=plt.subplots(ncols=3,nrows=2, figsize=[25,12])

for var,ax in zip(plot_var_list,AXES.flatten()):
    DF=pd.concat( [JU_panda[var].copy() for JU_panda in JU_panda_list],axis=1 )
    DF.coloumns=full_config_names

    DF.plot(legend=False,ax=ax,color=plot_colors,lw=1.2)

    if var=='lai':
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,ls='',marker='',c='k',
                                        yerr=SITE_panda['lai_sd'])
    
    if var in ['canht','lai']:
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,ls='',marker='.',c='r')
    else:
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,lw=0.8,c='r')

    if var in ['nee','gpp']:
        #ax2=ax.twinx()
        (JU_panda_list[fsmc_index]['fsmc']*10.).plot(ax=ax,legend=False,grid=False,c='purple',lw=2)
        #ax2.set_ylim(0.1,1.1)
        
    ax.set_title(var.upper(),fontsize=25)
    ax.set_ylabel(var.upper()+' ('+units[var]+')')
    ax.set_xlim(Season_Limits)


handles,labels=AXES[0,0].get_legend_handles_labels()
lgd=FIG.legend(handles[:nCONFIGS],configs,loc=8,ncol=min(nCONFIGS,8))
title=FIG.suptitle('Brattleby Winter Wheat Calibration', fontsize=25)
FIG.subplots_adjust(bottom=0.1)

plt.savefig(plot_dir+'/Calibration2_fit_Season'+str(iSEAS)+'.png',\
             bbox_extra_artists=(lgd,title,), bbox_inches='tight')

plt.close()


quit()

