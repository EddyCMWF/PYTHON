import numpy as np
import os

package_dir=os.path.dirname(os.path.realpath(__file__))

def cmip5_runs():

    return    [['BCC','bcc-csm1-1'],
              ['BCC','bcc-csm1-1-m'],
              ['BNU','BNU-ESM'],
              ['CCCma','CanESM2'],
              ['CMCC','CMCC-CMS'],
              ['CNRM-CERFACS','CNRM-CM5'],
              ['CSIRO-BOM','ACCESS1-0'],
              ['CSIRO-BOM','ACCESS1-3'],
              ['CSIRO-QCCCE','CSIRO-Mk3-6-0'],
              ['INM','inmcm4'],
              ['IPSL','IPSL-CM5A-LR'],
              ['IPSL','IPSL-CM5A-MR'],
              ['IPSL','IPSL-CM5B-LR'],
              ['MIROC','MIROC5'],
              ['MIROC','MIROC-ESM'],
              ['MIROC','MIROC-ESM-CHEM'],
              ['MOHC','HadGEM2-CC'],
              ['MOHC','HadGEM2-ES'],
              ['MPI-M','MPI-ESM-LR'],
              ['MPI-M','MPI-ESM-MR'],
              ['MRI','MRI-CGCM3'],
              ['NASA-GISS','GISS-E2-H'],
              ['NASA-GISS','GISS-E2-H-CC'],
              ['NASA-GISS','GISS-E2-R'],
              ['NASA-GISS','GISS-E2-R-CC'],
              ['NCAR','CCSM4'],
              ['NCC','NorESM1-M'],
              ['NCC','NorESM1-ME'],
              ['NOAA-GFDL','GFDL-CM3'],
              ['NOAA-GFDL','GFDL-ESM2G'],
              ['NOAA-GFDL','GFDL-ESM2M'],
              ['NSF-DOE-NCAR','CESM1-BGC'],
              ['NSF-DOE-NCAR','CESM1-CAM5'],
              ['NSF-DOE-NCAR','CESM1-WACCM']]


        
def read_imogen_varfile(cmip5_runs,filename):
    n_cmip5=len(cmip5_runs)
    outdata=np.zeros(n_cmip5)
    infile=open(filename,'r')
    indata={}
    for line in infile.readlines():
        split=line.split()
        centre=split[1]
        model=split[2]
        data=split[0]
        indata[centre+model]=data

    for i_gcm in range(n_cmip5):
        outdata[i_gcm]=indata[cmip5_runs[i_gcm][0]+cmip5_runs[i_gcm][1]]

    return outdata

#def kappa(cmip5_runs,filename='/users/global/chg/imogen/build/imogen_vals/kappa.dat'):
def kappa(cmip5_runs,filename=package_dir+'/kappa.dat'):
    return read_imogen_varfile(cmip5_runs,filename)

#def lambda_o(cmip5_runs,filename='/users/global/chg/imogen/build/imogen_vals/lambda_o.dat'):
def lambda_o(cmip5_runs,filename=package_dir+'/lambda_o.dat'):
    return read_imogen_varfile(cmip5_runs,filename)

#def lambda_l(cmip5_runs,filename='/users/global/chg/imogen/build/imogen_vals/lambda_l.dat'):
def lambda_l(cmip5_runs,filename=package_dir+'/lambda_l.dat'):
    return read_imogen_varfile(cmip5_runs,filename)

#def ocean_frac(cmip5_runs,filename='/users/global/chg/imogen/build/imogen_vals/ocean_frac.dat'):
def ocean_frac(cmip5_runs,filename=package_dir+'/ocean_frac.dat'):
    return read_imogen_varfile(cmip5_runs,filename)

#def nu(cmip5_runs,filename='/users/global/chg/imogen/build/imogen_vals/nu.dat'):
def nu(cmip5_runs,filename=package_dir+'/nu.dat'):
    return read_imogen_varfile(cmip5_runs,filename)

def SCENARIOs():
    """
        Return a List containing SCENARIO tag names
    """
    SCENARIOs = [ '2deg', '1p5deg', '1p81p5deg' ]

def TILE_full_names():
    """
        Return a List of full names for the 17 tiles
    """
    TILEs = [ 'Deciduous-Broadleaf','Tropical-Broadleaf','Temperate-Broadleaf',
              'Deciduous-Needleleaf','Evergreen-Needleleaf',
              'C3-grass','C3-crop','C3-pasture',
              'C4-grass','C4-crop','C4-pasture',
              'Deciduous-shrub','Evergreen-shrub',
              'Urban','Inland-Water','Bare-Soil','Land-Ice']
    return TILEs

def TILE_short_names():
    """
        Return a List of short names for the 17 tiles
    """
    TILEs = [ 'Dec-BL','Trop-BL','Temp-BL',
              'Dec-NL','Ever-NL',
              'C3-grass','C3-crop','C3-past',
              'C4-grass','C4-crop','C4-past',
              'Dec-shrub','Ever-shrub',
              'Urban','Lake','Soil','Ice']
    return TILEs

def TILE_colours():
    """
        Return a List colours for the 17 tile classes 
    """
    TILE_colours = ['#5ccc06','#51b305','#469a04',
                   '#166e2d','#115925',
                   '#d6e309','#e3df83','#a8ffeb',
                   '#abb507','#b5b268','#86ccbc',
                   '#9a6d02','#cc9103',
                   '#ff0101','#015cf6','#4e3801','#dfdfdf']
    return TILE_colours

def GCMs():
    """
        Return a List containing GCM tag names
    """
    GCMs = [
            'CEN_BCC_MOD_bcc-csm1-1',
            'CEN_BCC_MOD_bcc-csm1-1-m',
            'CEN_BNU_MOD_BNU-ESM',
            'CEN_CCCma_MOD_CanESM2',
            'CEN_CMCC_MOD_CMCC-CMS',
            'CEN_CNRM-CERFACS_MOD_CNRM-CM5',
            'CEN_CSIRO-BOM_MOD_ACCESS1-0',
            'CEN_CSIRO-BOM_MOD_ACCESS1-3',
            'CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0',
            'CEN_INM_MOD_inmcm4',
            'CEN_IPSL_MOD_IPSL-CM5A-LR',
            'CEN_IPSL_MOD_IPSL-CM5A-MR',
            'CEN_IPSL_MOD_IPSL-CM5B-LR',
            'CEN_MIROC_MOD_MIROC5',
            'CEN_MIROC_MOD_MIROC-ESM',
            'CEN_MIROC_MOD_MIROC-ESM-CHEM',
            'CEN_MOHC_MOD_HadGEM2-CC',
            'CEN_MOHC_MOD_HadGEM2-ES',
            'CEN_MPI-M_MOD_MPI-ESM-LR',
            'CEN_MPI-M_MOD_MPI-ESM-MR',
            'CEN_MRI_MOD_MRI-CGCM3',
            'CEN_NASA-GISS_MOD_GISS-E2-H',
            'CEN_NASA-GISS_MOD_GISS-E2-H-CC',
            'CEN_NASA-GISS_MOD_GISS-E2-R',
            'CEN_NASA-GISS_MOD_GISS-E2-R-CC',
            'CEN_NCAR_MOD_CCSM4',
            'CEN_NCC_MOD_NorESM1-M',
            'CEN_NCC_MOD_NorESM1-ME',
            'CEN_NOAA-GFDL_MOD_GFDL-CM3',
            'CEN_NOAA-GFDL_MOD_GFDL-ESM2G',
            'CEN_NOAA-GFDL_MOD_GFDL-ESM2M',
            'CEN_NSF-DOE-NCAR_MOD_CESM1-BGC',
            'CEN_NSF-DOE-NCAR_MOD_CESM1-CAM5',
            'CEN_NSF-DOE-NCAR_MOD_CESM1-WACCM',
            ]

    return GCMs



def GCM_INFO():
    """
        Return a Dictionary containing relevant info for GCMs
    """
    GCMs = {
            'CEN_BCC_MOD_bcc-csm1-1':{},
            'CEN_BCC_MOD_bcc-csm1-1-m':{},
            'CEN_BNU_MOD_BNU-ESM':{},
            'CEN_CCCma_MOD_CanESM2':{},
            'CEN_CMCC_MOD_CMCC-CMS':{},
            'CEN_CNRM-CERFACS_MOD_CNRM-CM5':{},
            'CEN_CSIRO-BOM_MOD_ACCESS1-0':{},
            'CEN_CSIRO-BOM_MOD_ACCESS1-3':{},
            'CEN_CSIRO-QCCCE_MOD_CSIRO-Mk3-6-0':{},
            'CEN_INM_MOD_inmcm4':{},
            'CEN_IPSL_MOD_IPSL-CM5A-LR':{},
            'CEN_IPSL_MOD_IPSL-CM5A-MR':{},
            'CEN_IPSL_MOD_IPSL-CM5B-LR':{},
            'CEN_MIROC_MOD_MIROC5':{},
            'CEN_MIROC_MOD_MIROC-ESM':{},
            'CEN_MIROC_MOD_MIROC-ESM-CHEM':{},
            'CEN_MOHC_MOD_HadGEM2-CC':{},
            'CEN_MOHC_MOD_HadGEM2-ES':{},
            'CEN_MPI-M_MOD_MPI-ESM-LR':{},
            'CEN_MPI-M_MOD_MPI-ESM-MR':{},
            'CEN_MRI_MOD_MRI-CGCM3':{},
            'CEN_NASA-GISS_MOD_GISS-E2-H':{},
            'CEN_NASA-GISS_MOD_GISS-E2-H-CC':{},
            'CEN_NASA-GISS_MOD_GISS-E2-R':{},
            'CEN_NASA-GISS_MOD_GISS-E2-R-CC':{},
            'CEN_NCAR_MOD_CCSM4':{},
            'CEN_NCC_MOD_NorESM1-M':{},
            'CEN_NCC_MOD_NorESM1-ME':{},
            'CEN_NOAA-GFDL_MOD_GFDL-CM3':{},
            'CEN_NOAA-GFDL_MOD_GFDL-ESM2G':{},
            'CEN_NOAA-GFDL_MOD_GFDL-ESM2M':{},
            'CEN_NSF-DOE-NCAR_MOD_CESM1-BGC':{},
            'CEN_NSF-DOE-NCAR_MOD_CESM1-CAM5':{},
            'CEN_NSF-DOE-NCAR_MOD_CESM1-WACCM':{},
            }

    return GCM_info
