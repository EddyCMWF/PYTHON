#!/bin/env python3.5

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plot_dir='/prj/GREENHOUSE/GREENHOUSE_sites/plots/Brattleby/crop_tuning_Season1/'
SITE_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
JULES_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/output/Brattleby/Calibration1a/'
Site='Brattleby'
tag='_crop_Season'

Lc_H2O = 2.501e6
kgC_to_umolsCO2_factor = (1./12.)* 1e9

wk_fit_vars=['lai','canht'] 
JU_fit_vars=['lai','canht']
SI_fit_vars=['lai','canht','lai_sd']
pft_index=5

#Read in SITE data:
SITE_veg_file=SITE_data_dir+'Brattleby_LAI_Canht_data.nc'

Sinf=nc.Dataset(SITE_veg_file,'r')
SITE_DICT={ var:Sinf.variables[var][:].squeeze() \
            for var in SI_fit_vars }

SITE_time=nc.num2date(Sinf.variables['date'][:],        \
                      units=Sinf.variables['date'].units)
SITE_panda=pd.DataFrame(SITE_DICT,index=SITE_time)
Sinf.close()

SITE_index=((SITE_panda.index>'2011-10-01')&(SITE_panda.index<'2012-09-01'))# |\
SITE_panda=SITE_panda[SITE_index]
#
#LAI_index=((SITE_panda.index>'2011-10-01')&(SITE_panda.index<'2012-06-15')) #|\
#SITE_panda['lai']=SITE_panda['lai'][LAI_index]
#SITE_panda['lai_sd']=SITE_panda['lai_sd'][LAI_index]

fit_params=['ttveg','alpro','alpst','alple','gamma','kappa'] 
stats=['mean','max','correlation','rmse','rmse_rel','stddev']

configs=[]
config_names=[]
ttvegs=['1700','1800','1900']
#plottag='alro_alst_'+ttveg+'ttv'
#plottag='alro_alst_ttv'
colors  = []

#ttvegs = ('1700','1800','1900')
#ttvegs = (['1700'])
nttvegs=len(ttvegs)
nttvegs=0

alros= ('19.7','19.8','19.9','20.0', '20.1', '20.2')
alros= ('20.1','20.2','20.3')
#alros= ('20.2', '20.4')
alros= [('20.0')]
nalros=len(alros)
nalros=0

alsts= ('15.2', '15.3', '15.4', '15.5', '15.6', '15.7')
alsts= ('15.4', '15.5', '15.6','15.7','15.8')
#alsts= ('15.4', '15.6', '15.8')
alsts= [('15.5')]# , '15.8')
nalsts=len(alsts)
nalsts=0

for ittveg in range(nttvegs):
   for ialst in range(nalsts):
      for ialro in range(nalros):
        iGREEN=ittveg
        alro=alros[ialro]
        alst=alsts[ialst]
        ttveg=ttvegs[ittveg]
        configs.append(ttveg+'_'+alro+'_'+alst+'_18.0_27.3_1.40')
        #config_names.append('TT$_{VEG}$='+ttveg+r', $\alpha_{root}$='+alro+r', $\alpha_{stem}$='+alst)
        config_names.append('TT$_{VEG}$='+ttveg)
        colors.append( ( 0.2,    \
                         (float(ittveg)/(nttvegs+1))+0.3, \
                         0.2   ) )
#        colors.append( ( float(ialro+1)/(nalros+1.)+0.1,    \
#                         0.3+(0.2*iGREEN) , \
#                         float(ialst+1)/(nalsts+1.)+0.1   ) )


#configs += ['1900_21.6_15.7_18.0_27.3_1.40']
#configs += ['1900_21.6_15.7_18.0_27.3_1.40']
#configs += ['1900_21.6_15.7_18.0_27.3_1.40']
#configs += ['1900_21.6_15.7_18.0_27.3_1.40']

configs += ['1700_20.2_15.5_18.0_27.3_1.40']
config_names += [r'TT$_{VEG}$=1700, $\alpha_{root}$=20.2, $\alpha_{stem}$=15.5, $\alpha_{leaf}$=18.0']
colors.append( (0.5,0.8,0.5) )

configs += ['844_18.5_16.0_18.0_27.3_1.40']
config_names += ['Summer_Wheat']
colors.append( (0,0,0) )
#colors='k'

plottag='Summer_Wheat' 
plottag='Calibration1' 
plottag='Calibration1_wtfullLAI' 
print(colors)
nCONFIGS =len(configs)

JU_panda_list=[]
for iSEAS in [1]:
    for config in configs:
        Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                config + \
                '.day.nc'
        if iSEAS==1:
            print(config)
        print(Jinfile)
        Jinf=nc.Dataset(Jinfile,'r')
        temp_time=nc.num2date( Jinf.variables['time'][:],         \
                                units=Jinf.variables['time'].units )
        temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                     for wkvar,invar in zip(wk_fit_vars,JU_fit_vars) }
        temp_dict['nee']=(Jinf.variables['npp_gb'][:].squeeze() - \
                          Jinf.variables['resp_s_gb'][:].squeeze())  \
                           * kgC_to_umolsCO2_factor
        Jinf.close()
        for var in wk_fit_vars:
            temp_dict[var]=temp_dict[var][:,pft_index]
        
        temp_panda=pd.DataFrame(temp_dict,index=temp_time)
        temp_index=((temp_panda.index>'2011-10-01')&(temp_panda.index<'2012-10-01')) 
        temp_panda=temp_panda[temp_index]
        JU_panda_list.append(temp_panda.copy()) 
            

full_config_names= [ config+'_Seas1' for config in configs ] #+ \
#                   [  config+'_Seas2' for config in configs ] 

LAI_panda = pd.concat( [JU_panda['lai'].copy() for JU_panda in JU_panda_list],axis=1 )
LAI_panda.columns=full_config_names

canht_panda = pd.concat( [JU_panda['canht'].copy() for JU_panda in JU_panda_list],axis=1 )
canht_panda.columns=full_config_names

#long_color_list=['b','g','y','c','orange','aqua','khaki','gold','grey','lime','cadetblue','olive','coral','brown','teal','palegreen','wheat']
#long_color_list[nCONFIGS-1]='k'

plot_colors = colors  #long_color_list#[:nCONFIGS] #+ long_color_list[:nCONFIGS]

FIG,AXES=plt.subplots(ncols=2,nrows=1, figsize=[20,5])

#ttl=FIG.suptitle('Crop Model - TTVEG sensitivity',fontsize=25)

LAI_panda.plot(legend=False,ax=AXES[0],color=plot_colors,lw=2)
SITE_panda['lai'].plot(ax=AXES[0],legend=False,grid=True,ls='',marker='',c='k',
                        yerr=SITE_panda['lai_sd'])
SITE_panda['lai'].plot(ax=AXES[0],legend=False,grid=True,ls='',marker='.',c='r')
AXES[0].set_title('LAI',fontsize=25)
AXES[0].set_ylabel('LAI (m$^2$ m$^{-2}$)')
AXES[0].set_ylim([0,8])

canht_panda.plot(legend=False,ax=AXES[1],color=plot_colors,lw=2)
SITE_panda['canht'].plot(ax=AXES[1],legend=False,grid=True,ls='',marker='.',c='r')
AXES[1].set_title('Canopy Height',fontsize=25)
AXES[1].set_ylabel('Canopy Height (m)')
AXES[1].set_ylim([0,1.0])


handles,labels=AXES[1].get_legend_handles_labels()
lgd=FIG.legend(handles[:nCONFIGS],config_names,loc=8,ncol=nCONFIGS)
FIG.subplots_adjust(bottom=0.2) #,top=0.85)

plt.savefig(plot_dir+'SensitivityPlot_'+plottag+'.png',\
              bbox_extra_artists=(lgd,), bbox_inches='tight')
#              bbox_extra_artists=(ttl,lgd,), bbox_inches='tight')

plt.close()


