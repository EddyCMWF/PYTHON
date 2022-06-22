#!/usr/bin/env python
#
########################################################
#
# Module storing the metadata for the various datasets 
#  required for the SC module.
#
########################################################

def select_source(sources,Message='Available Sources:'):
    print Message 
    for i in range(len(sources)):
        print str(i)+': '+sources[i]

    iSOURCE = input('Select an option:\n')

    return iSOURCE
    

def LAI_source_info(JULES_SOURCE_INFO):
    
    if (len(JULES_SOURCE_INFO['LAI_sources'][:])==1):
        iLAI = 0
    else:
        iLAI = select_source(JULES_SOURCE_INFO['LAI_sources'][:], \
                              Message='Select an LAI source: ')
    
    LAI_source = JULES_SOURCE_INFO['LAI_sources'][iLAI]
    
    if (LAI_source=='JULES_LAI'):
        LAI_file   = JULES_SOURCE_INFO['file']
        LAI_tsteps = JULES_SOURCE_INFO['tsteps']
        LAI_t_res  = JULES_SOURCE_INFO['t_res']
        LAI_nc_name= 'lai'
        LAI_tag    = ''
    elif ('static' in LAI_source) and (JULES_SOURCE_INFO['grid'] == 'WFDEI'):
        LAI_file   = '/users/eow/edwcom/WFD_EI/qrparm.veg.funcNew3.3.nc'
        LAI_tsteps = 1 
        LAI_t_res  = 'yearly'
        LAI_nc_name= 'field1392'
        LAI_tag    = 'STATlai'
    elif ('constant' in LAI_source):
        LAI_file   = 'Not_Applicable'
        LAI_tsteps = 1
        LAI_t_res  = 'yearly'
        LAI_nc_name= 'Not_Applicable'
        LAI_tag    = 'CONSTlai'
    elif (LAI_source=='JULES_pheno_LAI')  and (JULES_SOURCE_INFO['grid'] == 'WFDEI'):
        LAI_file   = '/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' + \
                     'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno.monthly_mean.nc'
        LAI_tsteps = JULES_SOURCE_INFO['tsteps']
        LAI_t_res  = JULES_SOURCE_INFO['t_res']
        LAI_nc_name= 'lai'
        LAI_tag    = 'JPHENOlai'
    elif ('climatology' in LAI_source) and (JULES_SOURCE_INFO['grid'] == 'CHESS'):
        LAI_file   = '/prj/chess/data/1km/v1.0/ancil/chess_monthly_lai_UK.nc'
        LAI_tsteps = 12 
        LAI_t_res  = 'monthly'
        LAI_nc_name= 'lai'
        LAI_tag    = 'CLIMlai'

    
    LAI_source_info = { 'source': LAI_source, \
                        'file': LAI_file,     \
                        'tsteps': LAI_tsteps, \
                        't_res': LAI_t_res,   \
                        'nc_name':LAI_nc_name,\
                        'tag':LAI_tag         \
                        }
    
    return LAI_source_info


def PFT_source_info(JULES_SOURCE_INFO):
   
    if (len(JULES_SOURCE_INFO['PFT_sources'][:])==1):
        iPFT=0
    else:
        iPFT = select_source(JULES_SOURCE_INFO['PFT_sources'][:], \
                              Message='Select an PFT source: ')
    
    PFT_source = JULES_SOURCE_INFO['PFT_sources'][iPFT]
    
    if (PFT_source == 'JULES_PFT'):
        PFT_file   = JULES_SOURCE_INFO['file']
        PFT_tsteps = JULES_SOURCE_INFO['tsteps']
        PFT_t_res  = JULES_SOURCE_INFO['t_res']
        PFT_nc_name= 'frac'
        PFT_tag    = ''
    elif ('static' in PFT_source) and (JULES_SOURCE_INFO['grid'] == 'WFDEI'):
        PFT_file   = '/users/eow/edwcom/WFD_EI/qrparm.veg.fracNew.nc'
        PFT_tsteps = 1 
        PFT_t_res  = 'yearly'
        PFT_nc_name= 'field1391'
        PFT_tag    = 'STATpft'
    elif (PFT_source=='JULES_pheno_PFT')  and (JULES_SOURCE_INFO['grid'] == 'WFDEI'):
        PFT_file   = '/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' + \
                     'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno.monthly_mean.nc'
        PFT_tsteps = JULES_SOURCE_INFO['tsteps']
        PFT_t_res  = JULES_SOURCE_INFO['t_res']
        PFT_nc_name= 'frac'
        PFT_tag    = 'JPHENOpft'
    elif ('static' in PFT_source) and (JULES_SOURCE_INFO['grid'] == 'CHESS'):
        PFT_file   = '/prj/chess/data/1km/v1.0/ancil/chess_landcover_2000.nc'
        PFT_tsteps = 1 
        PFT_t_res  = 'yearly'
        PFT_nc_name= 'frac'
        PFT_tag    = 'STATpft'

    
    PFT_source_info = { 'source': PFT_source,  \
                        'file': PFT_file,      \
                        'tsteps': PFT_tsteps,  \
                        't_res': PFT_t_res,    \
                        'nc_name': PFT_nc_name,\
                        'tag': PFT_tag         \
                        }
    
    return PFT_source_info

def SOIL_source_info(JULES_SOURCE_INFO):

    if (len(JULES_SOURCE_INFO['SOIL_sources'][:])==1):
        iSOIL=0
    else:
        iSOIL=select_source(JULES_SOURCE_INFO['SOIL_sources'][:], \
                            Message='Select a SOIL source: ')
    
    SOIL_source = JULES_SOURCE_INFO['SOIL_sources'][iSOIL]
    
    if (SOIL_source=='JULES_SOIL'):
        SOIL_file  = JULES_SOURCE_INFO['file']
        SOIL_tsteps = JULES_SOURCE_INFO['tsteps']
        SOIL_t_res  = JULES_SOURCE_INFO['t_res']
        SOIL_wilt_nc_name = 'sm_wilt'
        SOIL_sat_nc_name  = 'sm_sat'
        SOIL_tag  = ''
    elif (SOIL_source=='CHESS_ancil_BC'):
        SOIL_file   = '/prj/chess/data/1km/v1.0/ancil/chess_soilparams_hwsd_bc.nc'
        SOIL_tsteps = 1
        SOIL_t_res  = 'yearly'
        SOIL_wilt_nc_name = 'vwilt'
        SOIL_sat_nc_name  = 'vsat'
        SOIL_tag = ''

    SOIL_source_info = {'source':SOIL_source,            \
                        'file':SOIL_file,                \
                        'tsteps':SOIL_tsteps,            \
                        't_res':SOIL_t_res,              \
                        'wilt_nc_name':SOIL_wilt_nc_name,\
                        'sat_nc_name':SOIL_sat_nc_name,  \
                        'tag':SOIL_tag                   \
                        }
    
    return SOIL_source_info


def jules_sources_info():
    
    JULES_SOURCES = [ 'JULES-WFDEI pheno QS Global',        \
                      'JULES-WFDEI pheno QS TrifLAI Global',\
                      'JULES-WFDEI pheno TS Global',        \
                      'JULES-WFDEI pheno TS TrifLAI Global',\
                      'JULES-WFDEI BigSpin Global',         \
                      'JULES-WFDEI QuickSpin nocomp Global',\
                      'JULES-WFDEI TinySpin nocomp Global', \
                      'JULES-WFDEI BigSpin nocomp Global',  \
                      'JULES-WFDEI Koven Method',           \
                      ]
        

    JULES_SOURCES_info ={ 'JULES-WFDEI pheno TS Global':                \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' + \
                                   'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno_TS.monthly_mean.nc',    \
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                      \
                            'tag':'JULESv43_WFDEI_pheno_Global',                                \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],           \
                            'PFT_sources': ['JULES_PFT','static_PFT'],                          \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                                \
                          'JULES-WFDEI pheno QS Global':             \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' + \
                                   'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno_QS.monthly_mean.nc',    \
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                      \
                            'tag':'JULESv43_WFDEI_pheno_QS_Global',                             \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],           \
                            'PFT_sources': ['JULES_PFT','static_PFT'],                          \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                           \
                          'JULES-WFDEI pheno QS TrifLAI Global':     \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' +     \
                                   'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno_QS_trifLAI.monthly_mean.nc',\
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                          \
                            'tag':'JULESv43_WFDEI_pheno_QS_TrifLAI_Global',                         \
                            'LAI_sources':['JULES_LAI','static_LAI','constant_LAI'],                \
                            'PFT_sources': ['JULES_PFT','static_PFT'],                              \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                         \
                          'JULES-WFDEI pheno TS TrifLAI Global':        \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_pheno/' +   \
                                   'JULES_v4.3_WFDEI_RsQ10_GLOBAL_pheno_TS_trifLAI.monthly_mean.nc', \
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                        \
                            'tag':'JULESv43_WFDEI_pheno_TrifLAI_Global',                          \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],             \
                            'PFT_sources': ['JULES_PFT','static_PFT'],                            \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                         \
                          'JULES-WFDEI BigSpin Global':              \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_TRIFFID/' +  \
                                   'JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc',\
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                         \
                            'tag':'JULESv43_WFDEI_BigSpin_Global',                                 \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI','JULES_pheno_LAI'], \
                            'PFT_sources': ['JULES_PFT','static_PFT','JULES_pheno_PFT'],                \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                         \
                          'JULES-WFDEI QuickSpin nocomp Global':     \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_nocomp/' +  \
                                   'J4.3_WFDEI_TRIF_Q10_GLOBAL_QS_nocomp.monthly_mean.nc',             \
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                             \
                            'tag':'J43_WFDEI_QS_nocomp_Global',                                        \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],\
                            'PFT_sources': ['JULES_PFT'],               \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                         \
                          'JULES-WFDEI TinySpin nocomp Global':     \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_nocomp/' + \
                                   'J4.3_WFDEI_TRIF_Q10_GLOBAL_TS_nocomp.monthly_mean.nc',             \
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',                             \
                            'tag':'J43_WFDEI_TS_nocomp_Global',                                        \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],\
                            'PFT_sources': ['JULES_PFT'],               \
                            'SOIL_sources':['JULES_SOIL']                                            \
                            },                         \
                          'JULES-WFDEI BigSpin nocomp Global':     \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/WFDEI_nocomp/' + \
                                   'JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin_nocomp.monthly_mean.nc',\
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',              \
                            'tag':'J43_WFDEI_BS_nocomp_Global',                         \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],   \
                            'PFT_sources': ['JULES_PFT'],                  \
                            'SOIL_sources':['JULES_SOIL']                               \
                            },                         \
                          'JULES-WFDEI Koven Method':     \
                          { 'file':'/users/eow/edwcom/SC_simulator/JULES_output/KovenMethod/' + \
                                   'JULES_v4.3_TRIFFID_RsQ10_WFDEI_GLOBAL_KovenMethod_spin2.monthly_mean.nc',\
                            'grid':'WFDEI','tsteps':396,'t_res':'monthly',              \
                            'tag':'J43_WFDEI_Koven_Method',                         \
                            'LAI_sources': ['JULES_LAI','static_LAI','constant_LAI'],   \
                            'PFT_sources': ['JULES_PFT'],                  \
                            'SOIL_sources':['JULES_SOIL']                               \
                            },                         \
                      }

    return JULES_SOURCES, JULES_SOURCES_info


def jules_WFDEI_sources():
    std_JULES_dir = '/prj/GREENHOUSE/SC_simulator/JULES_output/'
    SC_dir = '/prj/GREENHOUSE/SC_simulator/OPERATIONAL_output/'
    JULES_SOURCES_AND_FILENAMES = \
            [  [ 'WFDEI_nocomp',std_JULES_dir+'WFDEI_nocomp/' +                               \
                    'JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin_nocomp.monthly_mean.nc',\
                    'blue' ], \
               [ 'WFDEI_KovenMethod', std_JULES_dir+'KovenMethod/' +                             \
                    'JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_KovenMethod_spin2.monthly_mean.nc', 
                    'green' ], \
               [ 'WFDEI_TRIFFID', std_JULES_dir+'WFDEI_TRIFFID/' +                     \
                    'JULES_v4.3_WFDEI_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc', \
                    'seagreen' ], \
               [ 'WFDEI_SC_spin1', SC_dir + 'work_WFDEI_GLOBAL/' + \
                    'J4.3_WFDEI_GLOBAL_spin1.monthly_mean.nc', \
                    'gold' ],   \
               [ 'WFDEI_SC_spin2', SC_dir + 'work_WFDEI_GLOBAL/' + \
                    'J4.3_WFDEI_GLOBAL_spin2.monthly_mean.nc', \
                    'goldenrod' ],   \
               [ 'WFDEI_SC_spin3', SC_dir + 'work_WFDEI_GLOBAL/' + \
                    'J4.3_WFDEI_GLOBAL_spin3.monthly_mean.nc', \
                    'orange' ],   \
               [ 'WFDEI_SC_10y_spin1', SC_dir + 'work_WFDEI_GLOBAL_10yr/' + \
                    'J4.3_WFDEI_GLOBAL_10yr_spin1.monthly_mean.nc', \
                    'pink' ],   \
               [ 'WFDEI_SC_10y_spin2', SC_dir + 'work_WFDEI_GLOBAL_10yr/' + \
                    'J4.3_WFDEI_GLOBAL_10yr_spin2.monthly_mean.nc', \
                    'fuchsia' ],   \
               [ 'WFDEI_SC_10y_spin3', SC_dir + 'work_WFDEI_GLOBAL_10yr/' + \
                    'J4.3_WFDEI_GLOBAL_10yr_spin3.monthly_mean.nc', \
                    'red' ],   \
               [ 'WFDEI_SC_20y_spin1', SC_dir + 'work_WFDEI_GLOBAL_20yr/' + \
                    'J4.3_WFDEI_GLOBAL_20yr_spin1.monthly_mean.nc', \
                    'pink' ],   \
               [ 'WFDEI_SC_20y_spin2', SC_dir + 'work_WFDEI_GLOBAL_20yr/' + \
                    'J4.3_WFDEI_GLOBAL_20yr_spin2.monthly_mean.nc', \
                    'fuchsia' ],   \
               [ 'WFDEI_SC_20y_spin3', SC_dir + 'work_WFDEI_GLOBAL_20yr/' + \
                    'J4.3_WFDEI_GLOBAL_20yr_spin3.monthly_mean.nc', \
                    'red' ],   \
                    ]
    return JULES_SOURCES_AND_FILENAMES
     


def jules_CHESS_sources():
    std_JULES_dir = '/prj/GREENHOUSE/SC_simulator/JULES_output/'
    SC_dir        = '/prj/GREENHOUSE/SC_simulator/OPERATIONAL_output/'
    JULES_SOURCES_AND_FILENAMES =  [  \
               [ 'CHESS_KovenMethod', std_JULES_dir+'CHESS_KovenMethod/' +            \
                    'JULES_v4.3_TRIFFID_RsQ10_CHESS_KovenMethod_spin2.monthly_mean.nc', 
                    'green' ], \
               [ 'CHESS_SC_10yr_spin1', SC_dir + 'work_CHESS_10yr/' + \
                    'J4.3_CHESS_spin1.monthly_mean.nc', \
                    'pink' ],   \
               [ 'CHESS_SC_10yr_spin2', SC_dir + 'work_CHESS_10yr/' + \
                    'J4.3_CHESS_spin2.monthly_mean.nc', \
                    'fuchsia' ],   \
               [ 'CHESS_SC_10yr_spin3', SC_dir + 'work_CHESS_10yr/' + \
                    'J4.3_CHESS_spin3.monthly_mean.nc', \
                    'red' ],   \
               [ 'CHESS_SC_20yr_spin1', SC_dir + 'work_CHESS_20yr/' + \
                    'J4.3_CHESS_spin1.monthly_mean.nc', \
                    'pink' ],   \
               [ 'CHESS_SC_20yr_spin2', SC_dir + 'work_CHESS_20yr/' + \
                    'J4.3_CHESS_spin2.monthly_mean.nc', \
                    'fuchsia' ],   \
               [ 'CHESS_SC_20yr_spin3', SC_dir + 'work_CHESS_20yr/' + \
                    'J4.3_CHESS_spin3.monthly_mean.nc', \
                    'steelblue', 'BC-Composite' ],   \
               [ 'CHESS_SC_BT_spin3', SC_dir + 'work_CHESS_BT/' + \
                    'J4.3_CHESS_BT_spin3.monthly_mean.nc', \
                    'olive', 'Broadleaf Tree' ],   \
               [ 'CHESS_SC_NT_spin3', SC_dir + 'work_CHESS_NT/' + \
                    'J4.3_CHESS_NT_spin3.monthly_mean.nc', \
                    'darkgreen', 'Needleleaf Tree' ],   \
               [ 'CHESS_SC_grass_spin3', SC_dir + 'work_CHESS_grass/' + \
                    'J4.3_CHESS_grass_spin3.monthly_mean.nc', \
                    'goldenrod', 'Grass' ],   \
               [ 'CHESS_SC_shrub_spin3', SC_dir + 'work_CHESS_shrub/' + \
                    'J4.3_CHESS_shrub_spin3.monthly_mean.nc', \
                    'saddlebrown', 'Shrub' ],                 \
               [ 'CHESS_SC_VG_spin3', SC_dir + 'work_CHESS_VG/' + \
                    'J4.3_CHESS_VG_spin3.monthly_mean.nc', \
                    'violet', 'VG-Composite' ],     \
               [ 'CHESS_SC_BC_spin3', SC_dir + 'work_CHESS_BC/' + \
                    'J4.3_CHESS_BC_spin3.monthly_mean.nc', \
                    'blue', 'BC-Composite' ],     \
               [ 'CHESS_SC_ConstSoil_spin3', SC_dir + 'work_CHESS_ConstSoil/' + \
                    'J4.3_CHESS_ConstSoil_spin3.monthly_mean.nc', \
                    'grey', '1Soil-Composite' ],     \
               [ 'CHESS_SC_BC_tdelta1_spin3', SC_dir + 'work_CHESS_BC_tdelta1/' + \
                    'J4.3_CHESS_BC_tdelta1_spin3.monthly_mean.nc', \
                    'fuchsia', 'BC-Tdelta1' ],     \
                                ]

    return JULES_SOURCES_AND_FILENAMES

def jules_sources_info_chess():
    
    JULES_SOURCES = [ 'JULES-CHESS Alberto run (monthly)'   \
                      ]

    JULES_SOURCES_info = { 'JULES-CHESS Alberto run (monthly)':                             \
                           { 'file':'/users/global/albmar/CHESS/outputs/' +                 \
                                   'chess_4.1.month.nc',                                   \
                             'grid':'CHESS','tsteps':624,'t_res':'monthly',                 \
                             'tag':'J42_CHESS_albmar_nocomp',                               \
                             'LAI_sources': ['climatology_LAI','static_LAI','constant_LAI'],\
                             'PFT_sources': ['static_PFT'],                                 \
                             'SOIL_sources': ['CHESS_ancil_BC']                             \
                            }                          \
                          }

    return JULES_SOURCES, JULES_SOURCES_info

