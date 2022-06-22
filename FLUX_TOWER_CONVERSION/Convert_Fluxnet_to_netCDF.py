#!/bin/env python 
import netCDF4 as nc
import glob
import datetime as dt
import pandas as pd
import numpy as np
import maths_tools.DateTimeTools as DTT
import maths_tools.FluxSiteTools as FST
import scipy.ndimage as im
from maths_tools import MetTools as MT

import matplotlib.pylab as plt


IN_DIR='/data/grp/fluxdata/temp_extract_dir/'
OUT_DIR=IN_DIR+'converted/'

#SITES = ['UKAMo','UKEBu','UKESa','UKGri','UKHam','UKHar','UKLBT','UKPL3','UKTad']
#site_longnames = ['Authencorth Moss','Easter Bush','East Saltoun','Griffin','Alice Holt',\
#                   'Harwood','London BT Tower','Pang Lambourne (Forest)','Tadham Moor']


SITES = ['UKHam']# 'UKAMo']#,'UKHam']
START_YEARS = [2005] #[2002,2004]
END_YEARS = [2005]  #[2014,2012]
SITE_LONGNAMES = ['Authencorth Moss','Alice Holt']
nSITES = len(SITES)

TIME_NAME_IN='TIMESTAMP_START'
TIME_NAME_OUT='time'

VARNAMES_IN = ['SW_IN','NETRAD','PA','TA','H2O','P','WS','CO2',\
                'FC','SH','H', 'PPFD_IN', \
                'TS','TS_2','TS_3',\
               ]

VARNAMES_OUT= ['sw_down','rad_net','pstar','t','q_molmr','precip','wind','co2_atm', \
                'co2_flux','sensible_heat','latent_heat', 'ppfd_in', \
                't_soil1','t_soil2','t_soil3',            \
               ]

fill_value=-9999.

for iSITE in range(nSITES):
    site=SITES[iSITE]
    site_longname=SITE_LONGNAMES[iSITE]
    start_year=START_YEARS[iSITE]
    end_year=END_YEARS[iSITE]

    DATA_DICT = { varname:[] for varname in VARNAMES_OUT }
    TIME = []
    for year in range(start_year,end_year+1):
        print(year)
        infile=glob.glob(IN_DIR+'EFDC_L2_Flx_'+site+'_'+str(year)+'*.csv')[0]
        
        inlines=open(infile).readlines()    
    
        headers=inlines.pop(0)
        headers=headers[:-1].split(',')
    
        for line in inlines:
            split=line[:-1].split(',')
            for invarname,outvarname in zip(VARNAMES_IN,VARNAMES_OUT):
                if invarname in headers:
                    DATA_DICT[outvarname].append(split[headers.index(invarname)])
                else:
                    DATA_DICT[outvarname].append(fill_value)
            TIME.append(split[headers.index(TIME_NAME_IN)])
        
    TIME =np.array([ dt.datetime( int(datestring[:4]),   \
                                  int(datestring[4:6]),  \
                                  int(datestring[6:8]),  \
                                  int(datestring[8:10]), \
                                  int(datestring[10:12]) ) \
                       for datestring in TIME              ])
                     
    TIME_secs=nc.date2num(TIME,units='seconds since '+str(TIME[0]))
    

    for varname in VARNAMES_OUT:
        DATA_DICT[varname]=np.array([float(val) for val in DATA_DICT[varname]])
        DATA_DICT[varname]=np.ma.masked_equal(DATA_DICT[varname],fill_value)


    DF=pd.DataFrame(DATA_DICT,index=TIME)
    DF['t_soil']=DF[['t_soil1','t_soil2','t_soil3']].mean(axis=1)
    DF['t_soil']=im.median_filter(DF['t_soil'],10)
    del DF['t_soil1']
    del DF['t_soil2']
    del DF['t_soil3']

    GPP_data, TER_data, GPPFitParams =    \
           FST.GPP_from_NEE_SW_T(np.ma.masked_invalid(DF['co2_flux'].values),\
                                 np.ma.masked_invalid(DF['sw_down'].values), \
                                 np.ma.masked_invalid(DF['t_soil'].values),  \
                                 spike_filter='median',\
                                 FIT_maxfev=2000,      \
                                 ReturnResp=True,      \
                                 ReturnFitParams=True, fill_value=-999 )
    
    temp_DF = pd.DataFrame( {'gpp':GPP_data,'ter':TER_data},index=TIME )

    DF = pd.concat([DF,temp_DF],axis=1)



    OUT_DATA_DICT={}

    # Fill missing Pressure data with 101.3 kPa
    temp_data=DATA_DICT['pstar'].data.copy()
    temp_data[DATA_DICT['pstar'].mask==True]=101300.
    OUT_DATA_DICT['pstar']={ 'data':temp_data, 'units':'Pa' }
    del temp_data


    # Temperature
    # Convert t to K
    for var in ['t','t_soil']:
        print(var)
        #temp_data=DATA_DICT[var].copy()+273.15
        temp_data=np.ma.masked_invalid(DF[var].values).copy()+273.15
        badex=np.where(temp_data.mask==True)[0]
        #print(badex.shape)
        # fill with previous year temperature
        loop_cnt=0
        max_loop=10
        while (len(badex)>0) and (loop_cnt<=max_loop):
            temp_data[badex]=temp_data[badex-(48*365)]
            badex=np.where(temp_data.mask==True)[0]
            loop_cnt+=1
    
        # fill with previous day temperature
        while len(badex)>0:
            temp_data[badex]=temp_data[badex-48]
            badex=np.where(temp_data.mask==True)[0]
    
        OUT_DATA_DICT[var]={'data':temp_data.copy(),'units':'K'}
    
        plt.plot(TIME,temp_data,label='filled')
        plt.plot(TIME,DF[var]+273.15,label='raw')
        plt.legend()
        plt.show()


    # fill RH with mean diurnal cycle values
    var='q_molmr'
    print(var)
    temp_data=DATA_DICT[var].reshape(-1,48).copy()
    diurnal_cycle = temp_data.mean(axis=0)#
    for ihh in range(48):
        temp_data_hh=temp_data[:,ihh]
        temp_data_hh[temp_data_hh.mask==True]=diurnal_cycle[ihh]
        temp_data[:,ihh]=temp_data_hh
    
    
    # Convert to Specific Humidity
    spec_hum=MT.mixr2sh(temp_data.flatten()*18e-3/29)
    
    OUT_DATA_DICT['q']={'data':spec_hum,'units':'kg kg^-1'}
    
    #print(diurnal_cycle.shape)
    plt.plot(TIME,spec_hum,label='filled')
    plt.plot(TIME,MT.mixr2sh(DATA_DICT[var]*18e-3/29),label='raw')
    plt.legend()
    plt.show()
    
    #Out_DataDict[var]=temp_RH.flatten()

    # Radiation Fields
    for var in ['rad_net','sw_down']:
        print(var)
        temp_data=DATA_DICT[var].copy()
        badex=np.where(temp_data.mask==True)[0]
        # First fill with previous year data, then next year, then first year etc.
        loop_cnt=0
        max_loop=10
        while (len(badex)>0) and (loop_cnt<=max_loop):
            loop_cnt+=1
            replace_index=badex-(48*365) 
            replace_index[replace_index<0]+=(48*365)
            temp_data[badex]=temp_data[badex-(48*365)]
            badex=np.where(temp_data.mask==True)[0]
            if len(badex)>0: 
                replace_index=badex+(48*365) 
                replace_index[replace_index>=len(temp_data)]-=(48*365)
                temp_data[badex]=temp_data[replace_index]
                badex=np.where(temp_data.mask==True)[0]
            
        # fill with previous day temperature
        while len(badex)>0:
            temp_data[badex]=temp_data[badex-48]
            badex=np.where(temp_data.mask==True)[0]

        if var=='sw_down':
            temp_data[temp_data<0]=0.
            
        OUT_DATA_DICT[var]={'data':temp_data.copy(),'units':'W m^-2'}
        
        plt.plot(TIME,temp_data,label='filled')
        plt.plot(TIME,DATA_DICT[var],label='raw')
        plt.legend()
        plt.show()
        

    # Fill precip with zeros
    var='precip'
    print(var)

    temp_data=DATA_DICT[var].copy()/1800.
    temp_data[temp_data.mask==True]=0

    plt.plot(TIME,temp_data,label='filled')
    plt.plot(TIME,DATA_DICT[var]/1800.,label='raw')
    plt.legend()
    plt.show()

    OUT_DATA_DICT[var]={'data':temp_data,'units':'kg m^-2 s^-1'}


    var='wind'
    print(var)
    temp_data=DATA_DICT[var].copy()
    badex=np.where(temp_data.mask==True)[0]
    goodex=np.where(temp_data.mask==False)[0]
    filled_data=np.interp(TIME_secs,TIME_secs[goodex],temp_data[goodex])
    OUT_DATA_DICT[var]={'data':filled_data,'units':'m s^-1'}

    plt.plot(TIME,filled_data,label='filled')
    plt.plot(TIME,DATA_DICT['wind'],label='raw')
    plt.legend()
    plt.show()

    
    var='co2_atm'
    print(var)
    temp_data=DATA_DICT[var].copy()
    # Filter out less than 330 ppm
    temp_data[temp_data<330]=temp_data.fill_value
    temp_data=np.ma.masked_equal(temp_data.data,temp_data.fill_value)*(44e-6/29)

    badex=np.where(temp_data.mask==True)[0]
    goodex=np.where(temp_data.mask==False)[0]
    filled_data=np.interp(TIME_secs,TIME_secs[goodex],temp_data[goodex])
    OUT_DATA_DICT[var]={'data':filled_data,'units':'kg kg^-1'}

    plt.plot(TIME,filled_data,label='filled')
    plt.plot(TIME,DATA_DICT[var]*(44e-6/29),label='raw')
    plt.legend()
    plt.show()

    DF[ ['co2_flux','gpp','ter'] ].plot()

    vars=['co2_flux','gpp','ter']
    for var in vars:
        print(var)
        OUT_DATA_DICT[var]={'data':DF[var].values,'units':'umolC m^-2 s^-1'}
    
    vars=['sensible_heat','latent_heat'] 
    for var in vars:
        OUT_DATA_DICT[var]={'data':DF[var].values,'units':'W m^-2'}

    print('IN:')
    for var in DATA_DICT:
        print(var)

    print('\nOUT:')
    for var in OUT_DATA_DICT:
        print(var,OUT_DATA_DICT[var]['units'])

    #Output data to netCDF
    if start_year==end_year:
        OUTFILE=OUT_DIR+site+'_'+str(start_year)+'_JULES_drivedata.nc'
    else:
        OUTFILE=OUT_DIR+site+'_'+str(start_year)+'_'+str(end_year)+'_JULES_drivedata.nc'
    print('Producing output file: '+OUTFILE)
    outf=nc.Dataset(OUTFILE,'w')

    outf.createDimension('time',len(TIME))
    outf.createDimension('land',1)

    outvar=outf.createVariable('time','float32',('time'))
    outvar.units='seconds since '+str(TIME[0])
    outvar[:]=TIME_secs

    for var in OUT_DATA_DICT:
        outvar=outf.createVariable(var,'float32',('time','land'))
        outvar.units=OUT_DATA_DICT[var]['units']
        outvar[:]=OUT_DATA_DICT[var]['data']
    
    outf.title='Gap filled met data for the Alice Holt flux tower'
    
    outf.title='Non gap filled flux data for the Alice Holt flux tower'
    
    outf.owner='Edward Comyn-Platt (edwcom@ceh.ac.uk)'

    outf.close()

