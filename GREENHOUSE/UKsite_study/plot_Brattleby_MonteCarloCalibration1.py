#!/bin/env python3.5

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SITE_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
JULES_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/output/Brattleby/Calibration1/'
Site='Brattleby'
tag='_crop_Season'

wk_fit_vars=['lai','canht'] 
JU_fit_vars=['lai','canht']
SI_fit_vars=['lai','canht']
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

SITE_index=((SITE_panda.index>'2011-10-01')&(SITE_panda.index<'2012-06-15')) |\
           ((SITE_panda.index>'2014-10-01')&(SITE_panda.index<'2015-08-01'))
SITE_panda=SITE_panda[SITE_index]

TTVEG_opts=['1800','1900','1950','2000','2050']
alpro_opts=['20.0','20.5','21.0','21.5','22.0']
alpst_opts=['16.0','16.5','17.0','17.5','18.0']
alple_opts=['18.5','19.0','19.5','20.0','20.5']

nTTVEG=len(TTVEG_opts)
nALPRO=len(alpro_opts)
nALPST=len(alpst_opts)
nALPLE=len(alple_opts)

template_array=np.zeros([nTTVEG,nALPRO,nALPST,nALPLE])
stats=['mean','max','correlation','rmse','rmse_rel','stddev']

STAT_dict = { var: { stat:template_array.copy() for stat in stats } \
                 for var in wk_fit_vars }

fit_params=['ttveg','alpro','alpst','alple'] 
PARAM_dict= { param:template_array.copy() for param in fit_params } 

for itvg in range(nTTVEG):
 for ialpro in range(nALPRO):
  for ialpst in range(nALPST):
   for ialple in range(nALPLE):
       
       param_vals=[TTVEG_opts[itvg],alpro_opts[ialpro],alpst_opts[ialpst],alple_opts[ialple]]
       for param,val in zip(fit_params,param_vals):
           PARAM_dict[param][itvg,ialpro,ialpst,ialple]=float(val)

       pd_list=[]
    
       for iSEAS in [1,2]:
           Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                     TTVEG_opts[itvg]+'_'+   \
                     alpro_opts[ialpro]+'_'+ \
                     alpst_opts[ialpst]+'_'+ \
                     alple_opts[ialple]+     \
                     '.day.nc'
           #print(Jinfile)
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
           STAT_dict[var]['mean'][itvg,ialpro,ialpst,ialple]=(JU_panda)[var].mean()
           STAT_dict[var]['max'][itvg,ialpro,ialpst,ialple]=(JU_panda)[var].max()
           STAT_dict[var]['rmse'][itvg,ialpro,ialpst,ialple]=(JU_panda-SITE_panda)[var].abs().mean()
           STAT_dict[var]['stddev'][itvg,ialpro,ialpst,ialple]=(JU_panda-SITE_panda)[var].abs().std()
           STAT_dict[var]['correlation'][itvg,ialpro,ialpst,ialple]=JU_panda[var].corr(SITE_panda[var])
           STAT_dict[var]['rmse_rel'][itvg,ialpro,ialpst,ialple] = \
                   ((JU_panda-SITE_panda)/(JU_panda+SITE_panda))[var].abs().mean()*2.


#fig,ax=plt.subplots(ncols=1,nrows=1)
#markers=['.','x','^','o','s','D']
#for itvg in range(nTTVEG):
#    for ialpro in range(nALPRO):
#        for ialpst in range(nALPST):
#            for ialple in range(nALPLE):
#                ax.plot(STAT_dict['lai']['correlation'][itvg,ialpro,ialpst,ialple],    \
#                         STAT_dict['lai']['rmse'][itvg,ialpro,ialpst,ialple],           \
#                         ls='',marker=markers[itvg],                                    \
#                         color=( float(ialpro+1)/6.,float(ialpst+1)/6.,float(ialple+1)/6.) )
                    

lai_corr_min=0.91
canht_corr_min=0.99
lai_rmse_max=0.75
canht_rmse_max=0.08
lai_max_min=5.0
canht_max_min=0.6
index=np.where(  (STAT_dict['lai']['correlation']>lai_corr_min)     \
               & (STAT_dict['canht']['correlation']>canht_corr_min) \
               & (STAT_dict['lai']['rmse']<lai_rmse_max)            \
               & (STAT_dict['canht']['rmse']<canht_rmse_max)        \
               & (STAT_dict['lai']['max']>lai_max_min)        \
               & (STAT_dict['canht']['max']>canht_max_min)        )

nCONFIGS=len(index[0])

out_fname='/prj/GREENHOUSE/GREENHOUSE_sites/plots/Brattleby/crop_tuning/Calibration1_fit.txt'
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


# Plot selected configs
configs_temp = [ PARAM_dict[param][index] for param in fit_params ]

configs = [ str('%4i'%configs_temp[0][iconf])+'_'+    \
            str('%3.1f'%(configs_temp[1][iconf]))+'_'+ \
            str('%3.1f'%(configs_temp[2][iconf]))+'_'+ \
            str('%3.1f'%(configs_temp[3][iconf]))      \
              for iconf in range(nCONFIGS) ]
configs += ['1950_19.0_16.0_18.0']
nCONFIGS+=1

JU_panda_list=[]
for iSEAS in [1,2]:
    for config in configs:
        Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
                config + \
                '.day.nc'
        Jinf=nc.Dataset(Jinfile,'r')
        temp_time=nc.num2date( Jinf.variables['time'][:],         \
                                units=Jinf.variables['time'].units )
        temp_dict={ wkvar:Jinf.variables[invar][:].squeeze() \
                     for wkvar,invar in zip(wk_fit_vars,JU_fit_vars) }
        Jinf.close()
        for var in wk_fit_vars:
            temp_dict[var]=temp_dict[var][:,pft_index]
        
        temp_panda=pd.DataFrame(temp_dict,index=temp_time)
        JU_panda_list.append(temp_panda.copy()) 
            

full_config_names= [ config+'_Seas1' for config in configs ] + \
                   [  config+'_Seas2' for config in configs ] 


LAI_panda = pd.concat( [JU_panda['lai'].copy() for JU_panda in JU_panda_list],axis=1 )
LAI_panda.columns=full_config_names

canht_panda = pd.concat( [JU_panda['canht'].copy() for JU_panda in JU_panda_list],axis=1 )
canht_panda.columns=full_config_names

long_color_list=['b','g','k','y','c','orange','aqua','khaki','gold','grey']

plot_colors = long_color_list[:nCONFIGS] + long_color_list[:nCONFIGS]

FIG,AXES=plt.subplots(ncols=2,nrows=1, figsize=[20,8])

LAI_panda.plot(legend=False,ax=AXES[0],color=plot_colors,lw=2)
SITE_panda['lai'].plot(ax=AXES[0],legend=False,grid=True,ls='',marker='.',c='r')
AXES[0].set_title('LAI',fontsize=25)
AXES[0].set_ylabel('LAI (m$^2$ m$^{-2}$)')

canht_panda.plot(legend=False,ax=AXES[1],color=plot_colors,lw=2)
SITE_panda['canht'].plot(ax=AXES[1],legend=False,grid=True,ls='',marker='.',c='r')
AXES[1].set_title('Canopy Height',fontsize=25)
AXES[1].set_ylabel('Canopy Height (m)')

handles,labels=AXES[1].get_legend_handles_labels()
lgd=FIG.legend(handles[:nCONFIGS],configs,loc=8,ncol=min(nCONFIGS,5))
FIG.subplots_adjust(bottom=0.2)

plt.savefig('/prj/GREENHOUSE/GREENHOUSE_sites/plots/Brattleby/crop_tuning/Calibration1_fit.png',\
              bbox_extra_artists=(lgd,), bbox_inches='tight')

plt.close()


