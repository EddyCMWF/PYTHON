#!/bin/env python3.5

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

plot_dir='/prj/GREENHOUSE/GREENHOUSE_sites/plots/Brattleby/crop_tuning_Season2/'
SITE_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/data/Brattleby/'
JULES_data_dir='/prj/GREENHOUSE/GREENHOUSE_sites/output/Brattleby/Calibration2/'
Site='Brattleby'
tag='_crop_Season'
iSEAS=2

Season_Limits=('2015-03-01','2015-09-01')

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
SITE_flux_file2=SITE_data_dir+'ConCrop_Flux_data.nc'
Sinf=nc.Dataset(SITE_flux_file2,'r')
SITE_time=nc.num2date( Sinf.variables['time'][:], \
                       units=Sinf.variables['time'].units )
SITE_dict={ wkvar:Sinf.variables[invar][:].squeeze()  \
        for wkvar,invar in zip(wk_fit_vars[2:],SI2_fit_vars) }
Sinf.close()
SITE_panda=pd.DataFrame(SITE_dict,index=SITE_time).resample('D').mean()
##################

SITE_panda=pd.concat([SITE_panda,SITE_VEG_panda],axis=1)

SITE_index=((SITE_panda.index>Season_Limits[0])&(SITE_panda.index<Season_Limits[1]))

SITE_panda=SITE_panda[SITE_index]
SITE_panda['nee']*=-1.0

##################################################


# Plot selected configs
configs = [ ]
colors = [ ]
config_names = [ ]

configs += ['1900_20.3_15.9_18.0_27.3_1.40']
#config_names += [r'TT$_{VEG}$=1700, $\alpha_{root}$=20.2, $\alpha_{stem}$=15.5, $\alpha_{leaf}$=18.0']
config_names += [r'$\mu$=0.05, DVI$_{senescence}$=2.00, DVI$_{critical}$=1.50']
colors.append( (0.5,0.8,0.5) )
#colors= (0.5,0.8,0.5) 

configs += ['0.025_1.55_1.45']
config_names += [r'$\mu$=0.025, DVI$_{senescence}$=1.45, DVI$_{critical}$=1.55']
colors.append( (0.4,0.4,0.6) )

#configs += ['844_18.5_16.0_18.0_27.3_1.40']
#config_names += ['Summer_Wheat']
#colors.append( (0,0,0) )

plottag='step1'
plottag='step2'


#configs += ['0.02_1.50_1.00']
nCONFIGS =len(configs)

JU_panda_list=[]
for config in configs:
    Jinfile=JULES_data_dir+Site+tag+str(iSEAS)+'_'+ \
            config + \
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
    
    if 'fsmc' in Jinf.variables:
        FSMC=Jinf.variables['fsmc'][:,pft_index].squeeze()
        temp_dict['gpp']/=FSMC
        temp_dict['nee']/=FSMC

    Jinf.close()

    temp_panda=pd.DataFrame(temp_dict,index=temp_time)
    temp_index=((temp_panda.index>Season_Limits[0])&(temp_panda.index<Season_Limits[1]))
    temp_panda=temp_panda[temp_index]
    
    JU_panda_list.append(temp_panda.copy()) 


full_config_names= config_names 

#long_color_list=['b','g','k','y','c','orange','aqua','khaki','gold','grey']
#plot_colors = long_color_list[:nCONFIGS] 
#plot_colors = [ cm.get_cmap('Spectral')(fraction) for fraction in np.arange(0,1.00000001,1./nCONFIGS) ]
plot_colors=colors

plot_var_list=['lai','canht','shf','lhf','nee','gpp']
units={'shf':'W m$^{-2}$','lhf':'W m$^{-2}$','nee':'umolCO2 m$^{-2}$ $s^{-1}$',\
        'lai':'m$^2$ m$^{-2}$','canht':'m','gpp':'umolCO2 m$^{-2}$ $s^{-1}$'}
titles={'shf':'Sensible Heat Flux','lhf':'Latent Heat Flux',\
        'nee':'Net Ecosystem Exchange','gpp':'Gross Primary Productivity',\
        'lai':'Leaf Area Index','canht':'Canopy Height'}
ylims={'shf':[-50,150],'lhf':[0,250],\
       'nee':[-20,30],'gpp':[0,30],\
       'lai':[0,8],'canht':[0,1.]}

FIG,AXES=plt.subplots(ncols=2,nrows=3, figsize=[16,13])

for var,ax in zip(plot_var_list,AXES.flatten()):
    DF=pd.concat( [JU_panda[var].copy() for JU_panda in JU_panda_list],axis=1 )
    DF.coloumns=full_config_names

    DF.plot(legend=False,ax=ax,color=plot_colors,lw=2)

    if var=='lai':
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,ls='',marker='',c='k',
                                        yerr=SITE_panda['lai_sd'])
    
    if var in ['canht','lai']:
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,ls='',marker='.',c='r')
    else:
        SITE_panda[var].plot(ax=ax,legend=False,grid=True,lw=1.2,c='r')

    #if var in ['nee','gpp']:
        #(JU_panda_list[-1]['fsmc']*10.).plot(ax=ax,legend=False,grid=False,c='purple',lw=2)

    ax.set_title(titles[var],fontsize=20)
    ax.set_ylabel(units[var])
    ax.set_xlim(Season_Limits)
    ax.set_ylim(ylims[var])


handles,labels=AXES[0,0].get_legend_handles_labels()
lgd=FIG.legend(handles[:nCONFIGS],full_config_names,loc=8,ncol=nCONFIGS)
#title=FIG.suptitle('Brattleby Winter Wheat Calibration', fontsize=25)
FIG.tight_layout(h_pad=1.0,w_pad=2.0,pad=1.0)
FIG.subplots_adjust(bottom=0.08)

plt.savefig(plot_dir+'/Calibration2_'+plottag+'_Season'+str(iSEAS)+'.png',\
             bbox_extra_artists=(lgd,)) #, bbox_inches='tight')
             #bbox_extra_artists=(lgd,title,), bbox_inches='tight')

plt.close()


quit()


