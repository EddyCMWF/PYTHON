#!/bin/env python 
import shapefile
import netCDF4 as nc
import pyproj
import numpy as np

# ECP Modules:
from SoilTools import BrooksCorey as BC
from SoilTools import VanGenuchten as VG
from SoilTools import Thermal as TH

##################################################################################################
#Define Directories and Filenames
SOIL_DIR='/prj/GREENHOUSE/SOIL_PROPERTIES/datasets/England_Wales_Soil_Data/LDE16_12_SRUC_Tarsitano/'
shp_file=SOIL_DIR+'Spatial Soil Data/NATMAP1000.shp'
Hyd_LUT_file=SOIL_DIR+'Tabular Attribute Data/HORIZONhydraulics_V2.csv'
Cst_LUT_file=SOIL_DIR+'Tabular Attribute Data/HORIZONfundamentals_V2.csv'

LatLon_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'

OUTFILE=SOIL_DIR+'EnW_Soil_WeightedCompositionProperties_CHESSgrid_ECOSSElayers_PG.nc'
##################################################################################################
# Options
fill_value=-9999.0

# Max number of Series to sum over
nSERIES=10
# LU Groups (Correspond to CHESS classes in  CHESS_LU_idx
LU_groups=   ['PG']  #,'AR','OT'] # ['LE']  #
CHESS_LU_idx= [range(8)]  # [[2],[4],[0,1,3,5,6,7]] # [range(8)]  #
nLU=len(LU_groups)

# Define Soil Depths
# JULES
#CH_soil_depths=np.array([0.1,0.35,1.0,3.0])
#CH_mean_soil_depths=np.array([0.05,0.225,0.675,2.0])

# ECOSEE
CH_soil_depths=np.arange(0.05,3.0001,0.05)
CH_mean_soil_depths=CH_soil_depths-0.025

nSD=len(CH_soil_depths)

# Define input/output Parameter names to extract
Comp_in_names=['SAND_TOTAL','SILT', 'CLAY', 'OC', 'PH','BULK_DENSITY'] 
Comp_out_names=['sand','silt','clay','org_carb','ph','Bulk_Density'] 

BC_Pnames=['bexp','sathh','satcon','sm_sat','sm_wilt','sm_crit','hcap','hcon']
VG_Pnames=['oneovernminusone','oneoveralpha','ksat','sm_sat','sm_wilt','sm_crit','hcap','hcon']

##################################################################################################
#read in latlon/xy data from chess_landcover
LLinf=nc.Dataset(LatLon_file,'r')
lats=LLinf.variables['lat'][:]
lons=LLinf.variables['lon'][:]
landcover=LLinf.variables['frac'][:]
x=LLinf.variables['x'][:]
y=LLinf.variables['y'][:]
nx,ny=len(x),len(y)
#LLinf.close()
LAND_MASK=landcover.mask[0,:]

# Create Land Use Map correpsonding to the soil properties land uses:
CHESS_LU_map = np.zeros_like(landcover[:nLU,:])
for iLU in range(nLU):
    CHESS_LU_map[iLU,:,:]=np.sum(landcover[CHESS_LU_idx[iLU],:],axis=0)


#################################################################################################
# Read in the Look-up tables
Cst_LUT_lines=[line[:-1] for line in open(Cst_LUT_file).readlines()]
Hyd_LUT_lines=[line[:-1] for line in open(Hyd_LUT_file).readlines()]
Cst_headers=Cst_LUT_lines.pop(0).split(',')
Hyd_headers=Hyd_LUT_lines.pop(0).split(',')

LUT_Dict={hdr:[] for hdr in Cst_headers+Hyd_headers}
for Cst_line,Hyd_line in zip(Cst_LUT_lines,Hyd_LUT_lines):
    Cst_vals=Cst_line.split(',')
    Hyd_vals=Hyd_line.split(',')
    # Unless there are values for sand, silt, clay, OC and ph, throw line away
    if (Cst_vals[7]!='')&(Cst_vals[11]!='')&(Cst_vals[12]!='')&(Cst_vals[13]!='')&(Cst_vals[14]!='')\
                        &(Hyd_vals[7]!=''):
        # Append Components to relevant list
        for hdr,val in zip(Cst_headers,Cst_vals):
            if val=='':
                val='-1'
            LUT_Dict[hdr].append(val)
        # Append Hydraulics to list (missing out the series names and ids etc.)
        for hdr,val in zip(Hyd_headers[7:],Hyd_vals[7:]):
            if val=='':
                val='-1'
            LUT_Dict[hdr].append(val)

headers_to_float=['UPPER_DEPTH', 'LOWER_DEPTH', 'SAND_TOTAL', \
                  'SAND_FINE', 'SAND_MED', 'SAND_COARSE', 'SILT', 'CLAY', 'OC', 'PH',\
                  'BULK_DENSITY', 'PARTICLE_DENSITY', 'TOTAL_POROSITY',\
                  'THV0', 'THV1', 'THV5', 'THV10', 'THV40', 'THV200', 'THV1500', \
                  'KSAT_SUBVERT', 'KSAT_LAT', \
                  'VG_TH_S', 'VG_TH_R', 'VG_ALPHA', 'VG_N', 'VG_M', \
                  'BC_TH_S', 'BC_ALPHA', 'BC_BETA',                  ]

# convert to numpy arrays
for hdr in headers_to_float:
    LUT_Dict[hdr]=[float(val) for val in LUT_Dict[hdr]]
for hdr in LUT_Dict:
    LUT_Dict[hdr]=np.array(LUT_Dict[hdr])
    
# Calculate Mean Depth
LUT_Dict['MEAN_DEPTH']=(LUT_Dict['LOWER_DEPTH']+LUT_Dict['UPPER_DEPTH'])/200.
#################################################################################################

#################################################################################################
# Read in Shape file
SF=shapefile.Reader(shp_file)
field_names=[field[0] for field in SF.fields]
field_names.pop(0)
Eindex=field_names.index('EAST_1K')
Nindex=field_names.index('NORTH_1K')

ShRcs=SF.shapeRecords()
#################################################################################################
####################################################################################
## Loop around each shape, extracting the Soil Composition and 
# computing the BC and VG parameters
# All calculated values are weighted accoring to the soil type weights in shapefile
####################################################################################

# Set up dictionaries for stroing the data:
CHESS_SoilComp={ varname: np.zeros([nSD,ny,nx])+fill_value \
                    for varname in Comp_out_names  }
BC_Properties={ varname: np.zeros([nSD,ny,nx])+fill_value \
                    for varname in BC_Pnames              }
VG_Properties={ varname: np.zeros([nSD,ny,nx])+fill_value \
                    for varname in VG_Pnames              }

for iShRc in range(len(ShRcs)):
    ShRc=ShRcs[iShRc]
    x_index=int(float(ShRc.record[Eindex])/1000.)
    y_index=int(float(ShRc.record[Nindex])/1000.)
    
    #Check to see if location is a CHESS Land Point
    if (LAND_MASK[y_index,x_index]==True):
        continue
    
    #Create Dictionary to store point data in
    Point_Data={ varname: np.zeros([nSD]) for varname in Comp_out_names  }
    BC_Properties_P={ varname: np.zeros([nSD]) for varname in BC_Pnames   }
    VG_Properties_P={ varname: np.zeros([nSD]) for varname in VG_Pnames   }
    
    # Collate Soil Series names and percetage cover
    SERIES_dict={'names':[],'pc':[]}
    for iS in range(nSERIES):
        Sname_index=int(iS*2.)+2
        Spc_index=Sname_index+1
        if (ShRc.record[Spc_index]!=None)&\
            (ShRc.record[Sname_index] in LUT_Dict['SERIES_NAME']):
            SERIES_dict['names'].append(ShRc.record[Sname_index])
            SERIES_dict['pc'].append(ShRc.record[Spc_index])
    # Convert to numpy arrays for searching
    SERIES_dict['names']=np.array(SERIES_dict['names'])
    SERIES_dict['pc']=np.array(SERIES_dict['pc'])
    # Normalise percentages (i.e. ignoring the small and 'other' components)
    SERIES_dict['pc']=SERIES_dict['pc']/np.sum(SERIES_dict['pc'])
    
    # Count for Soil Types in this gridcell
    SOIL_SUM_COUNT=0
    # Loop over soil series for location
    for iS in range(len(SERIES_dict['names'])):
        # Loop over LU groups:
        LU_fraction=np.zeros(nSD)
        
        #Create Dictionary to store series data in
        Series_Data={ varname: np.zeros([nSD]) for varname in Comp_out_names  }
        BC_Properties_S={ varname: np.zeros([nSD]) for varname in BC_Pnames  }
        VG_Properties_S={ varname: np.zeros([nSD]) for varname in VG_Pnames  }
    
        for iLU in range(len(LU_groups)):
            # if Landuse not relevent then move on
            if CHESS_LU_map[iLU,y_index,x_index]==0:
                continue
            
            LU=LU_groups[iLU]
            # index locations which have correct series
            #print(SERIES_dict['names'][iS],LU)
            LU_S_index= np.where((LUT_Dict['LU_GROUP']==LU)      & \
                                 (LUT_Dict['SERIES_NAME']==SERIES_dict['names'][iS]))[0]
            #if there is no match LU and Soil Series then continue:
            #print(LU_S_index)
            if len(LU_S_index)==0:
                continue
            
            # loop over output soil depths
            for iSD in range(nSD):
                # find closest input SD
                if (len(LU_S_index)==0):
                    print(LU, SERIES_dict['names'][iS])
                    print(ShRc.record)
                SD_index=np.argmin(np.abs( LUT_Dict['MEAN_DEPTH'][LU_S_index] \
                                          -CH_mean_soil_depths[iSD]) )
                # INDEX of the closest soil depth for the series and LU
                INDEX=LU_S_index[SD_index]
                #print(iSD,INDEX)
                LU_fraction[iSD]+=CHESS_LU_map[iLU,y_index,x_index]
                # Loop over output parameters
                for inname,outname in zip(Comp_in_names,Comp_out_names):
                    #for each param, at each soil depth, sum the weighted components
                    Series_Data[outname][iSD]+= \
                      (LUT_Dict[inname][INDEX]*SERIES_dict['pc'][iS]) \
                    * CHESS_LU_map[iLU,y_index,x_index]
                 
                
                # get values required for BC and VG calculations
                sand=LUT_Dict[Comp_in_names[0]][INDEX]
                silt=LUT_Dict[Comp_in_names[1]][INDEX]
                clay=LUT_Dict[Comp_in_names[2]][INDEX]
                org_carb=LUT_Dict[Comp_in_names[3]][INDEX]
                #print(INDEX,sand,silt,clay)
                # Normalise percentages to 100%
                ratio=(sand+silt+clay)/100.
                sand/=ratio
                silt/=ratio
                clay/=ratio

                # Compute the Brooks and Corey Soil Properties
                BC_Properties_raw=BC.get_BC_soil_properties(clay,sand,silt)
                BC_Properties_raw['hcap']=TH.hcap(clay,sand,silt,sm_sat=BC_Properties_raw['sm_sat'])
                BC_Properties_raw['hcon']=TH.hcon_Farouki(clay,sand,silt,sm_sat=BC_Properties_raw['sm_sat'])
                #for prop in BC_Properties_raw:
                #    print(prop, BC_Properties_raw[prop])
                
                #print('BC_name: ',BC_name)
                #print('iS: ',iS)
                #print('iLU: ',iLU)
                #print('y_index,x_index: ',(y_index,x_index) )
                # Sum using the percentage weighting and land use cover:
                for BC_name in BC_Pnames:
                    BC_Properties_S[BC_name][iSD]+= \
                    (BC_Properties_raw[BC_name]*SERIES_dict['pc'][iS]) \
                                    * CHESS_LU_map[iLU,y_index,x_index]
                
                # Compute the Van Genuchten Soil Properties
                VG_Properties_raw=VG.get_VG_soil_props_from_comp(clay,sand,silt,\
                                                                 OC_PC=org_carb,\
                                                                 SUB_SOIL=bool(iSD) ) 
                
                VG_Properties_raw['hcap']=TH.hcap(clay,sand,silt,                    \
                                                  sm_sat=VG_Properties_raw['sm_sat'],\
                                                  l_vg=True                          )
                
                VG_Properties_raw['hcon']=TH.hcon_Farouki(clay,sand,silt,                    \
                                                          sm_sat=VG_Properties_raw['sm_sat'],\
                                                          l_vg=True                          )
                
                # Sum using the percentage weighting and land use cover:
                for VG_name in VG_Pnames:
                    VG_Properties_S[VG_name][iSD]+= \
                    (VG_Properties_raw[VG_name]*SERIES_dict['pc'][iS]) \
                    * CHESS_LU_map[iLU,y_index,x_index]
                        
        # Normalise for any absent Land Uses if required
        if (LU_fraction.min()>0.1):
            SOIL_SUM_COUNT+=1
            for outname in Comp_out_names:
                #for each param, at each soil depth, sum the weighted components
                Point_Data[outname] += Series_Data[outname]/LU_fraction
            for BC_name in BC_Pnames:
                BC_Properties_P[BC_name] += BC_Properties_S[BC_name]/LU_fraction
            for VG_name in VG_Pnames:
                VG_Properties_P[VG_name] += VG_Properties_S[VG_name]/LU_fraction
        else:
            continue
    
    # If no Soil Found, reset to fill value
    if SOIL_SUM_COUNT!=0:
        for outname in Comp_out_names:
            #for each param, at each soil depth, sum the weighted components
            CHESS_SoilComp[outname][:,y_index,x_index]=Point_Data[outname]
        for BC_name in BC_Pnames:
            BC_Properties[BC_name][:,y_index,x_index]=BC_Properties_P[BC_name]
        for VG_name in VG_Pnames:
            VG_Properties[VG_name][:,y_index,x_index]=VG_Properties_P[VG_name]


###############################################################
# Output data to netCDF file:
###############################
outf=nc.Dataset(OUTFILE,'w')

# Create Dimensions
outf.createDimension('x',nx)
outf.createDimension('y',ny)
outf.createDimension('z',nSD)

# write out dimension variables
for var in ['x','y']:
    outvar=outf.createVariable(var,'float32',(var))
    for att in LLinf.variables[var].ncattrs():
        outvar.setncattr(str(att),LLinf.variables[var].getncattr(str(att)))
        outvar[:]=LLinf.variables[var][:]
    
outvar=outf.createVariable('z','float32',('z'))
outvar.units='m'
outvar.long_name='lower soil depth'
outvar[:]=CH_soil_depths

# Loop around soil Composition parameters and 
for var in CHESS_SoilComp:
    outvar=outf.createVariable(var,'float32',('z','y','x'),fill_value=fill_value)
    if var=='ph':
        outvar.units='pH'
    elif var=='Bulk_Density':
        outvar.units='g cm^-3'
    else:
        outvar.units='percent'
    outvar[:]=CHESS_SoilComp[var]

# Loop around soil BC parameters and 
for var in BC_Properties:
    outvar=outf.createVariable('BC_'+var,'float32',('z','y','x'),fill_value=fill_value)
    outvar[:]=BC_Properties[var]

# Loop around soil VG parameters and 
for var in VG_Properties:
    outvar=outf.createVariable('VG_'+var,'float32',('z','y','x'),fill_value=fill_value)
    outvar[:]=VG_Properties[var]
    
outf.Title='Soil composition and properties map for England and Wales based on the Cranfield Soils Database'
outf.owner='Edward Comyn-Platt'
outf.note='VG and BC soil properties calculated at the soil type scale then aggregated up to the 1km resolution'

outf.close()

LLinf.close()

