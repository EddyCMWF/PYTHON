# The purpose of this code is to extract temperature values from rcp2.6 for Ed. 
# C. Huntingford (16th February 2017)
import numpy as np
import pickle, os, sys

cmip5_runs = [['BCC','bcc-csm1-1','rcp26'],
              ['BCC','bcc-csm1-1-m','rcp26'],
              ['BNU','BNU-ESM','rcp26'],
              ['CCCma','CanESM2','rcp26'],
              ['CMCC','CMCC-CMS','rcp45'],
              ['CNRM-CERFACS','CNRM-CM5','rcp26'],
              ['CSIRO-BOM','ACCESS1-0','rcp45'],
              ['CSIRO-BOM','ACCESS1-3','rcp45'],
              ['CSIRO-QCCCE','CSIRO-Mk3-6-0','rcp26'],
              ['INM','inmcm4','rcp45'],
              ['IPSL','IPSL-CM5A-LR','rcp26'],
              ['IPSL','IPSL-CM5A-MR','rcp26'],
              ['IPSL','IPSL-CM5B-LR','rcp45'],
              ['MIROC','MIROC5','rcp26'],
              ['MIROC','MIROC-ESM','rcp26'],
              ['MIROC','MIROC-ESM-CHEM','rcp26'],
              ['MOHC','HadGEM2-CC','rcp45'],
              ['MOHC','HadGEM2-ES','rcp26'],
              ['MPI-M','MPI-ESM-LR','rcp26'],
              ['MPI-M','MPI-ESM-MR','rcp26'],
              ['MRI','MRI-CGCM3','rcp26'],
              ['NASA-GISS','GISS-E2-H','rcp26'],
              ['NASA-GISS','GISS-E2-H-CC','rcp45'],
              ['NASA-GISS','GISS-E2-R','rcp26'],
              ['NASA-GISS','GISS-E2-R-CC','rcp45'],
              ['NCAR','CCSM4','rcp26'],
              ['NCC','NorESM1-M','rcp26'],
              ['NCC','NorESM1-ME','rcp26'],
              ['NOAA-GFDL','GFDL-CM3','rcp26'],
              ['NOAA-GFDL','GFDL-ESM2G','rcp26'],
              ['NOAA-GFDL','GFDL-ESM2M','rcp26'],
              ['NSF-DOE-NCAR','CESM1-BGC','rcp45'],
              ['NSF-DOE-NCAR','CESM1-CAM5','rcp26'],
              ['NSF-DOE-NCAR','CESM1-WACCM','rcp26']]
n_files_pickled = len([name for name in os.listdir('/users/global/chg/imogen/build/pickled_yearly')])
n_cmip5 = len(cmip5_runs)
if(n_files_pickled != n_cmip5*3):
    print 'Not getting all the GCMs' ; sys.exit()

# Now loop over the GCMs
dir_pickle = '/users/global/chg/imogen/build/pickled_yearly/'
i_initial=0 ; i_plot=1
for i in range(0, n_cmip5):
    centre_int=cmip5_runs[i][0]
    model_int=cmip5_runs[i][1]
    rcp_lower=cmip5_runs[i][2]

    # Get model data
    file_nam_pickle_historical_lower = 'CEN_'+centre_int+'_MOD_'+model_int+'_historical_'+rcp_lower+'.pickle' 
    timeseries_historical_lower = pickle.load( open(dir_pickle+file_nam_pickle_historical_lower, 'rb'))

    # Only include GCMs that are operating under rcp26. Also, only include the first ensemble member, and GCM of interest
    if (rcp_lower =='rcp26' and centre_int=='MOHC' and model_int=='HadGEM2-ES'):

        for i_ens in range(0,len(timeseries_historical_lower)):
            if i_ens == 0:
                tas_control_est = np.mean(timeseries_historical_lower[i_ens][4][0:40])
                # Estimate initial state (and subtract off) based on first 4 decades. Avoids jumps.
                tas_delta = timeseries_historical_lower[i_ens][4] - tas_control_est

                # yr_offset gets year of RF forcing lined up with that of timeseries 
                yr_start = timeseries_historical_lower[i_ens][7]

                # Open and write to file. Only print numbers up to year 2100
                f=open('ed_rcp26_global_temp_MOHC_HadGEM2-ES.dat', 'w')
                for i_print in range(0,tas_delta.shape[0]):
                    if yr_start+i_print <= 2100:
                        output_line = '{0:d}{1:s}{2:.4f}'.format(yr_start+i_print, '   ', tas_delta[i_print])
                        f.write(output_line + '\n')
                f.close()
sys.exit()
