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

SITE_index=((SITE_panda.index>'2011-10-01')&(SITE_panda.index<'2012-07-01'))# |\
SITE_panda=SITE_panda[SITE_index]

LAI_index=((SITE_panda.index>'2011-10-01')&(SITE_panda.index<'2012-06-15')) #|\
SITE_panda['lai']=SITE_panda['lai'][LAI_index]
SITE_panda['lai_sd']=SITE_panda['lai_sd'][LAI_index]


fit_params=['ttveg','alpro','alpst','alple','gamma','kappa'] 
stats=['mean','max','correlation','rmse','rmse_rel','stddev']

TTVEG_opts=['1700'] #, '1850'] #, '1900', '1950']#  ,'2000']
ALPRO_opts=['20.1', '20.2', '20.3' ]
ALPST_opts=['15.4', '15.5', '15.6', '15.7', '15.8']
ALPLE_opts=['18.0'] #, '19.5', '20.0'] 
#ALPRO_opts=['20.5', '21.5', '22.5']
#ALPST_opts=['14.0', '15.0', '16.0']
#ALPLE_opts=['18.0', '19.0', '20.0']
GAMMA_opts=['27.3'] #, '28.0', '28.5' ] #26.5
KAPPA_opts=['1.40'] #1.90', '2.00'] # '1.70', '1.80', '1.90' ]

nTTVEG=len(TTVEG_opts)
nALPRO=len(ALPRO_opts)
nALPST=len(ALPST_opts)
nALPLE=len(ALPLE_opts)
nGAMMA=len(GAMMA_opts)
nKAPPA=len(KAPPA_opts)

template_array=np.zeros([nTTVEG,nALPRO,nALPST,nALPLE,nGAMMA,nKAPPA])

STAT_dict = { var: { stat:template_array.copy() for stat in stats } \
                 for var in wk_fit_vars }

PARAM_dict= { param:template_array.copy() for param in fit_params } 

for ittveg in range(nTTVEG):
 for ialpro in range(nALPRO):
  for ialpst in range(nALPST):
   for ialple in range(nALPLE):
    for igamma in range(nGAMMA):
     for ikappa in range(nKAPPA):
       
       param_vals=[TTVEG_opts[ittveg],\
                  ALPRO_opts[ialpro],ALPST_opts[ialpst],ALPLE_opts[ialple],\
                  GAMMA_opts[igamma],KAPPA_opts[ikappa]]
       for param,val in zip(fit_params,param_vals):
           PARAM_dict[param][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=float(val)

       pd_list=[]
    
       for iSEAS in [1]:
           Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                     TTVEG_opts[ittveg]+'_'+   \
                     ALPRO_opts[ialpro]+'_'+   \
                     ALPST_opts[ialpst]+'_'+   \
                     ALPLE_opts[ialple]+'_'+   \
                     GAMMA_opts[igamma]+'_'+   \
                     KAPPA_opts[ikappa]+ \
                     '.day.nc'
           print(Jinfile)
           Jinf=nc.Dataset(Jinfile,'r')
           temp_time=nc.num2date( Jinf.variables['time'][:],         \
                                  units=Jinf.variables['time'].units )
           temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                        for wkvar,invar in zip(wk_fit_vars,JU_fit_vars) }
           Jinf.close()
           for var in wk_fit_vars:
               temp_dict[var]=temp_dict[var][:,pft_index]
           temp_panda=pd.DataFrame(temp_dict,index=temp_time)
           pd_list.append(temp_panda.copy()) 

       JU_panda=pd.concat(pd_list)
        
       for var in wk_fit_vars:
           Spanda=SITE_panda[var].dropna()
           Jpanda=JU_panda[var].dropna()
           DIFFpanda=Jpanda-Spanda
           STAT_dict[var]['mean'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=Jpanda.mean()
           STAT_dict[var]['max'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=Jpanda.max()
           STAT_dict[var]['rmse'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=DIFFpanda.abs().mean()
           STAT_dict[var]['stddev'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=DIFFpanda.abs().std()
           STAT_dict[var]['correlation'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]=JU_panda[var].corr(Spanda)
           STAT_dict[var]['rmse_rel'][ittveg,ialpro,ialpst,ialple,igamma,ikappa] = \
                   (DIFFpanda/Jpanda).abs().mean()
       print(STAT_dict['lai']['correlation'][ittveg,ialpro,ialpst,ialple,igamma,ikappa],   \
             STAT_dict['canht']['correlation'][ittveg,ialpro,ialpst,ialple,igamma,ikappa], \
             STAT_dict['lai']['rmse'][ittveg,ialpro,ialpst,ialple,igamma,ikappa],          \
             STAT_dict['canht']['rmse'][ittveg,ialpro,ialpst,ialple,igamma,ikappa]         )

lai_corr_min=0.96
canht_corr_min=0.98
lai_rmse_max=0.6
canht_rmse_max=0.08
lai_max_min=5
canht_max_min=0.6
index=np.where(  (STAT_dict['lai']['correlation']>lai_corr_min)     \
               & (STAT_dict['canht']['correlation']>canht_corr_min) \
               & (STAT_dict['lai']['rmse']<lai_rmse_max)            \
               & (STAT_dict['canht']['rmse']<canht_rmse_max)        \
               & (STAT_dict['lai']['max']>lai_max_min)        \
               & (STAT_dict['canht']['max']>canht_max_min)        )

nCONFIGS=len(index[0])
print(nCONFIGS)
out_fname=plot_dir+'Calibration1a_fit.txt'
outf=open(out_fname,'w')
outf.write('Fit criteria: \n')
outf.write('LAI minimum correlation: %6.3f\n'%lai_corr_min)
outf.write('LAI maximum rmse: %6.3f\n'%lai_rmse_max)
outf.write('Canht minimum correlation: %6.3f\n'%canht_corr_min)
outf.write('Canht maximum rmse: %6.3f\n'%canht_rmse_max)

outf.write('\n\n')
outf.write('Satisfying Configurations:\n')
for param in fit_params: 
    print(param)
    outf.write('%15a'%(param)+':'+nCONFIGS*' %10.2f ' % tuple(PARAM_dict[param][index]) +'\n')

for stat in ['correlation','rmse','max']:
    for var in ['lai','canht']:
        string=var+','+stat[:4]
        outf.write('%15a'%(string)+':'+nCONFIGS*' %10.4f '% tuple(STAT_dict[var][stat][index]) +'\n')


outf.close()
if nCONFIGS>15:
    quit()

# Plot selected configs
configs_temp = [ PARAM_dict[param][index] for param in fit_params ]

configs = [ str('%4i'%configs_temp[0][iconf])+'_'+  \
            str('%4.1f'%configs_temp[1][iconf])+'_'+  \
            str('%4.1f'%configs_temp[2][iconf])+'_'+  \
            str('%4.1f'%configs_temp[3][iconf])+'_'+  \
            str('%4.1f'%configs_temp[4][iconf])+'_'+  \
            str('%4.2f'%(configs_temp[5][iconf])) \
              for iconf in range(nCONFIGS) ]
configs=[ config.replace(' ','') for config in configs ]

#configs += ['844_18.5_16.0_18.0_27.3_1.40']
#configs += ['1900_20.5_16.5_19.0_27.3_1.40']
#configs += ['1850_21.5_16.0_19.1_27.5_1.90']
#configs += ['1850_20.5_15.0_18.0_28.5_1.90']
nCONFIGS =len(configs)

JU_panda_list=[]
for iSEAS in [1]:
    for config in configs:
        Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                config + \
                '.day.nc'
        if iSEAS==1:
            print(config)
        #print(Jinfile)
        Jinf=nc.Dataset(Jinfile,'r')
        temp_time=nc.num2date( Jinf.variables['time'][:],         \
                                units=Jinf.variables['time'].units )
        temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                     for wkvar,invar in zip(wk_fit_vars,JU_fit_vars) }
        temp_dict['nee']=(Jinf.variables['npp_gb'][:].squeeze() - \
                          Jinf.variables['resp_s_gb'][:].squeeze())  \
                           * kgC_to_umolsCO2_factor
        #if config=='1850_21.5_16.0_19.0_27.5_1.80':
        if 'fsmc' in Jinf.variables:
            temp_dict['fsmc']=Jinf.variables['fsmc'][:,pft_index].squeeze()
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

long_color_list=['b','g','y','c','orange','aqua','khaki','gold','grey','lime','cadetblue','olive','coral','brown','teal','palegreen','wheat']
long_color_list[nCONFIGS-1]='k'

plot_colors = long_color_list#[:nCONFIGS] #+ long_color_list[:nCONFIGS]

FIG,AXES=plt.subplots(ncols=2,nrows=1, figsize=[20,8])

LAI_panda.plot(legend=False,ax=AXES[0],color=plot_colors,lw=2)
SITE_panda['lai'].plot(ax=AXES[0],legend=False,grid=True,ls='',marker='',c='k',
                        yerr=SITE_panda['lai_sd'])
SITE_panda['lai'].plot(ax=AXES[0],legend=False,grid=True,ls='',marker='.',c='r')
AXES[0].set_title('LAI',fontsize=25)
AXES[0].set_ylabel('LAI (m$^2$ m$^{-2}$)')
AXES[0].set_ylim([0,7])

canht_panda.plot(legend=False,ax=AXES[1],color=plot_colors,lw=2)
SITE_panda['canht'].plot(ax=AXES[1],legend=False,grid=True,ls='',marker='.',c='r')
AXES[1].set_title('Canopy Height',fontsize=25)
AXES[1].set_ylabel('Canopy Height (m)')

handles,labels=AXES[1].get_legend_handles_labels()
lgd=FIG.legend(handles[:nCONFIGS],configs,loc=8,ncol=min(nCONFIGS,4))
FIG.subplots_adjust(bottom=0.3)

plt.savefig(plot_dir+'Calibration1a_fit.png',\
              bbox_extra_artists=(lgd,), bbox_inches='tight')

plt.close()


NEE_panda = pd.concat( [JU_panda['nee'].copy() for JU_panda in JU_panda_list],axis=1 )
NEE_panda.columns=full_config_names
NEE_panda.plot(legend=False)
#JU_panda_list[0]['fsmc'].plot(legend=False)
#JU_panda_list[nCONFIGS]['fsmc'].plot(legend=False)
plt.savefig(plot_dir+'Calibration1a_NEE.png',\
              bbox_inches='tight')

plt.close()

