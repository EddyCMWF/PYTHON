# Modules required:
#   netCDF4
#   numpy
#   calendar


def BVOC_sources():
    JULES_data_dir = '/prj/wetlands_africa/jules/JASMIN/BVOCs/JULES_OUTPUT/'
    EMEP4UK_data_dir = '/prj/wetlands_africa/jules/JASMIN/BVOCs/EMEP4UK_output/'
    # sources is a dictionaries containing: 
    #      [ name, { NAME, DIR, FILETAG, \
    #                tstep, file_period, start_year, end_year, 
    #                iso_name, terp_name, 
    #                scale_factor  } ]
    sources = [ \
            { 'NAME':'Jv4.3.1_std_BC', \
              'DIR':JULES_data_dir+'EMEP4UK_BC/', \
              'FILETAG':'EMEP4UK_BC.hourly.YYYYMM.nc', \
              'tstep':3600,\
              'file_period':'monthly',\
              'time_axis':0,\
              'start_year':2001,\
              'end_year':2013,\
              'iso_name':'isoprene_gb',\
              'terp_name':'terpene_gb',\
              'scale_factor':3600.*24.*365.*1000.*(68./60), \
              'gridfile':['/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc','land_index'],\
              },
            { 'NAME':'Jv4.3.1_std_VG',\
              'DIR':JULES_data_dir+'EMEP4UK_VG/', \
              'FILETAG':'EMEP4UK_VG.hourly.YYYYMM.nc', \
              'tstep':3600,\
              'file_period':'monthly',\
              'time_axis':0,\
              'start_year':2001,\
              'end_year':2013,\
              'iso_name':'isoprene_gb',\
              'terp_name':'terpene_gb',\
              'scale_factor':3600.*24.*365.*1000.*(68./60), \
              'gridfile':['/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc','land_index'],\
              },
            {'NAME':'EMEP4UK_rv4.3',\
             'DIR':EMEP4UK_data_dir+'rv4.3/', \
             'FILETAG':'EMEP4UK_UK_webrun_emep_4.3_YYYY_day_BVOCS.nc', \
             'tstep':86400,\
             'file_period':'yearly',\
             'time_axis':0,\
             'start_year':2001,\
             'end_year':2014,\
             'iso_name':'Emis_mgm2_BioNatC5H8',\
             'terp_name':'Emis_mgm2_BioNatAPINENE', \
             'scale_factor':365./1000., \
              'gridfile':[None,None],\
             },
            { 'NAME':'Jv4.6_std_VG_monthly',\
              'DIR':JULES_data_dir+'J4.6_EMEP4UK_VG/', \
              'FILETAG':'J4.6_EMEP4UK_VG.monthly.YYYY.nc', \
              'tstep':'monthly',\
              'file_period':'yearly',\
              'time_axis':0,\
              'start_year':2001,\
              'end_year':2014,\
              'iso_name':'isoprene',\
              'terp_name':'terpene',\
              'scale_factor':3600.*24.*365.*1000.*(68./60), \
              'gridfile':['/users/eow/edwcom/EMEP/EMEP4UK/EMEP4UK_JULES_output_index.nc','land_index'],\
              'fracfile':[JULES_data_dir+'J4.6_EMEP4UK_VG/J4.6_EMEP4UK_VG.dump.20010101.0.nc','frac',range(10)],\
              },
            {'NAME':'EMEP4UK_rv4.3_monthly',\
             'DIR':EMEP4UK_data_dir+'rv4.3/', \
             'FILETAG':'EMEP4UK_UK_webrun_emep_4.3_YYYY_monthly_BVOCS.nc', \
             'tstep':'monthly',\
             'file_period':'yearly',\
             'time_axis':0,\
             'start_year':2001,\
             'end_year':2014,\
             'iso_name':'Emis_mgm2_BioNatC5H8',\
             'terp_name':'Emis_mgm2_BioNatAPINENE', \
             'scale_factor':365./1000., \
              'gridfile':[None,None],\
             },
              ]

    return sources 
                
def read_netCDF_BVOC_data_to_TRES(SOURCE, TRES=86400, \
                                    start_year=None, end_year=None, \
                                    return_time_vector=False,\
                                    fill_value=-999.):
    import netCDF4 as nc
    import numpy as np
    import calendar
    import datetime as dt
    from maths_tools import DateTimeTools as DTT

    if start_year==None:
        start_year=SOURCE['start_year']

    if end_year==None:
        end_year=SOURCE['end_year']
    
    # number of TSTEPS to go in to each averaging period
    if TRES>0:
        nTSTEPS = TRES/SOURCE['tstep']
    elif TRES==-1:
        if SOURCE['tstep']=='monthly':
            nTSTEPS=1
        else:
            nTSTEPS=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
            nTSTEPS*=86400/SOURCE['tstep']

    if 'fracfile' in SOURCE:
        pft_index=SOURCE['fracfile'][2]
        frinf=nc.Dataset(SOURCE['fracfile'][0],'r')
        frac=frinf.variables[SOURCE['fracfile'][1]][pft_index,:]
        frinf.close()
    else:
        frac=1.

    iso_data=None
    terp_data=None

    if SOURCE['file_period']!='full_run':
        for year in range(start_year,end_year+1):
            FILETAG_YEAR = SOURCE['FILETAG'].replace('YYYY',str(year))
            print(FILETAG_YEAR)
            
            if calendar.isleap(year):
                n_secs_in_year = 366.*86400.
                if isinstance(nTSTEPS,list):
                    if len(nTSTEPS)==12:
                        nTSTEPS[1]+=86400/SOURCE['tstep']
            else:
                n_secs_in_year = 365.*86400.

            if SOURCE['file_period']=='monthly':
                iso_data_temp=None
                terp_data_temp=None
                for mnth in range(12):
                    FILETAG = FILETAG_YEAR.replace('MM',"%02d"%(mnth+1))
                    FILE=SOURCE['DIR']+FILETAG
                    inf=nc.Dataset(FILE,'r')
                    if iso_data_temp==None:
                        iso_data_temp=inf.variables[SOURCE['iso_name']][:].squeeze()
                    else:
                        iso_data_temp=np.append(iso_data_temp,\
                                                inf.variables[SOURCE['iso_name']][:].squeeze(),\
                                                axis=SOURCE['time_axis'])
                    
                    if terp_data_temp==None:
                        terp_data_temp=inf.variables[SOURCE['terp_name']][:].squeeze()
                    else:
                        terp_data_temp=np.append(terp_data_temp,\
                                                 inf.variables[SOURCE['terp_name']][:].squeeze(),\
                                                 axis=SOURCE['time_axis'])
                    inf.close()

            elif SOURCE['file_period']=='yearly':
                FILE=SOURCE['DIR']+FILETAG_YEAR
                print(FILE)
                inf=nc.Dataset(FILE,'r')
                iso_data_temp=inf.variables[SOURCE['iso_name']][:].squeeze()
                terp_data_temp=inf.variables[SOURCE['terp_name']][:].squeeze()
                inf.close()
            
            # if fracfile provided sum over desired pfts
            if 'fracfile' in SOURCE:
                 iso_data_temp= np.sum(iso_data_temp*frac,axis=1)
                 terp_data_temp= np.sum(terp_data_temp*frac,axis=1)


            # Perform the average of data on a yearly basis to avoid over usage of memory
            if isinstance(nTSTEPS,int):
                nBINS=int(iso_data_temp.shape[0]/nTSTEPS)
                if nTSTEPS>1:
                    new_dims=[nBINS,nTSTEPS]+list(iso_data_temp.shape[1:])
                    iso_data_temp=np.mean(iso_data_temp.reshape(new_dims),axis=1)
                    terp_data_temp=np.mean(terp_data_temp.reshape(new_dims),axis=1)
                elif nTSTEPS<1:
                    print('No interpolation method installed yet')
                    quit()
            elif isinstance(nTSTEPS,list):
                if len(nTSTEPS)==12:
                    print('No monthly averaging installed yet')
                    quit()
            #else:
                # Nothing required here
            
            # correction for when datasets have missing data (EMEP4UK)
            # Fix by adding nans at end of year
            if (TRES>0) & (float(TRES)*nBINS!=n_secs_in_year):
                # new number of tsteps 
                new_nBINS=int(n_secs_in_year/TRES)
                n_fill_BINS = new_nBINS-nBINS
                fill_dims=[n_fill_BINS]+list(iso_data_temp.shape[1:])
                fill_data=np.zeros(fill_dims)*np.nan
                iso_data_temp=np.append(iso_data_temp,fill_data,axis=SOURCE['time_axis'])
                terp_data_temp=np.append(terp_data_temp,fill_data,axis=SOURCE['time_axis'])
                nBINS=new_nBINS
            
            # stick into np.array for full period:
            if iso_data==None:
                iso_data=iso_data_temp.copy()
            else:
                iso_data=np.append(iso_data,iso_data_temp.copy(),axis=SOURCE['time_axis'])

            if terp_data==None:
                terp_data=terp_data_temp.copy()
            else:
                terp_data=np.append(terp_data,terp_data_temp.copy(),axis=SOURCE['time_axis'])

            del iso_data_temp
            del terp_data_temp
    
    # convert to regular 2D grid if neccesary
    if SOURCE['gridfile'][0]!=None:
        print SOURCE['gridfile'][0],SOURCE['gridfile'][1]
        grinf=nc.Dataset(SOURCE['gridfile'][0],'r')
        grindex=grinf.variables[SOURCE['gridfile'][1]][:]
        grinf.close()
        grimask=np.ones_like(grindex)
        # apply index to last dimension
        iso_data=iso_data[...,grindex]*grimask
        terp_data=terp_data[...,grindex]*grimask
    
    # now mask any invalid data points
    try:
        iso_mask= (iso_data.mask==True)|np.isnan(iso_data)
    except:
        iso_mask= np.isnan(iso_data)
    iso_data=np.ma.masked_array(iso_data,mask=iso_mask,fill_value=fill_value)
    iso_data.data[iso_data.mask==True]=fill_value

    try:
        terp_mask= (terp_data.mask==True)|np.isnan(terp_data)
    except:
        terp_mask= np.isnan(terp_data)
    terp_data=np.ma.masked_array(terp_data,mask=terp_mask,fill_value=fill_value)
    terp_data.data[terp_data.mask==True]=fill_value
    
    # Finally multiply by the sacaling factor to ensure units = g m^-2 yr^-1
    iso_data=iso_data*SOURCE['scale_factor']
    terp_data=terp_data*SOURCE['scale_factor']
    
    
    if return_time_vector:
        if TRES>0:
            time_vector_secs = np.arange(0,iso_data.shape[0]*TRES,TRES)
            #print time_vector_secs.shape, iso_data.shape[0]
            time_vector = nc.num2date(time_vector_secs, \
                                      units='seconds since '+str(start_year)+'-01-01 00:00:00', \
                                      calendar='standard')
        elif TRES==-1:
            time_vector=np.array(DTT.DTarange_months(dt.datetime(start_year,1,1),dt.datetime(end_year,12,1)))
            print(time_vector.shape)
        return iso_data, terp_data, time_vector
            
    else:
    
        return iso_data, terp_data




# END OF FUNCTIONS
                    

                

