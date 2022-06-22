#!/bin/env python
import shapefile
import netCDF4 as nc
import pyproj
import numpy as np
from shapely.geometry import Polygon

# Plotting Modules:
#import matplotlib.pyplot as plt
######################################################################################
def Square_Polygon(i,j,spacing):
    return Polygon( ( (i*spacing,j*spacing),        \
                      (i*spacing,j*spacing+spacing),    \
                      (i*spacing+spacing,j*spacing+spacing),\
                      (i*spacing+spacing,j*spacing),    \
                      (i*spacing,j*spacing),        \
                      )  )
    
######################################################################################
# Filenames, directories and options
SOIL_DIR='/prj/GREENHOUSE/SOIL_PROPERTIES/datasets/Scotland_Soil_Data/'
shpvar='soil_250k'
shp_file=SOIL_DIR+shpvar+'/'+shpvar+'.shp'

LatLon_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'

LUT_file=SOIL_DIR+'SSKIB_OCT13_orgcarbGapFilled.csv'

#OUTFILE=SOIL_DIR+'Scot_Soil_WeightedComposition_CHESSgrid.nc'
OUTFILE=SOIL_DIR+'Scot_Soil_WeightedComposition_CHESSgrid_ECOSSElayers.nc'
#####################################################################################

##################################################################################################
# Options
fill_value=-9999.0

SEMI=False
CULT=True 

if SEMI:
    OUTFILE=OUTFILE[:-3]+'_SEMI.nc'
elif CULT:
    OUTFILE=OUTFILE[:-3]+'_CULT.nc'

# Max number of Series to sum over
nSERIES=10

# Define Soil Depths
# JULES
#CH_soil_depths=np.array([0.1,0.35,1.0,3.0])
#CH_mean_soil_depths=np.array([0.05,0.225,0.675,2.0])

# ECOSEE
CH_soil_depths=np.arange(0.05,3.0001,0.05)
CH_mean_soil_depths=CH_soil_depths-0.025

nSD=len(CH_soil_depths)

# Define input/output Parameter names to extract
Comp_in_names=['SAND_MED', 'SILT_MED', 'CLAY_MED', 'CARBON_MED', \
               'organic_matter_med', 'PH_W_MED', 'NITROGEN_MED', \
               'CA_MED','MG_MED','NA_MED','K_MED','H_MED','TOTAL_P_MED', ]
                
Comp_out_names=['sand','silt','clay','org_carb',  \
                'org_matter','ph','nitrogen',     \
                'calcium','magnesium','sodium','potassium','hydrogen','phospherous']
                 
##################################################################################################

##################################################################################################
#read in latlon/xy data from chess_landcover
LLinf=nc.Dataset(LatLon_file,'r')
lats=LLinf.variables['lat'][:]
lons=LLinf.variables['lon'][:]
landcover=LLinf.variables['frac'][:]
x=LLinf.variables['x'][:]
y=LLinf.variables['y'][:]
nx,ny=len(x),len(y)

CHESS_cropFrac=landcover[4,:]
if SEMI:
    CHESS_cropFrac.data[CHESS_cropFrac.mask==False]=0.0
elif CULT:
    CHESS_cropFrac.data[CHESS_cropFrac.mask==False]=1.0
##################################################################################################

##################################################################################################
#Read in LUT
LUT_lines=open(LUT_file).readlines()
LUT_lines=[line[:-1] for line in LUT_lines]
headers=LUT_lines.pop(0).split(',')
sand_index=headers.index('SAND_MED')
silt_index=headers.index('SILT_MED')
clay_index=headers.index('CLAY_MED')
temp_LUT_Dict={hdr:[] for hdr in headers}
for line in LUT_lines:
    vals=line.split(',')
    # Unless there are values for sand, silt and clay, throw line away
    if (vals[sand_index]!='')&(vals[silt_index]!='')&(vals[clay_index]!=''):
        for hdr,val in zip(headers,vals):
            if val=='':
                val='-1'
            temp_LUT_Dict[hdr].append(val)

Req_Headers=['LAND_USE', 'SERIES_NAME', 'SERIES_CODE',              \
             'HORZ_TOP', 'HORZ_BOTTOM',                             \
             'CA_MED', 'MG_MED','NA_MED','K_MED','H_MED',           \
             'PH_W_MED','CARBON_MED','NITROGEN_MED','TOTAL_P_MED',  \
             'SAND_MED','SILT_MED','CLAY_MED', 'organic_matter_med' ]

LUT_Dict={hdr:temp_LUT_Dict[hdr] for hdr in Req_Headers}

for hdr in Req_Headers[2:]:
    LUT_Dict[hdr]=[float(val) for val in LUT_Dict[hdr]]
for hdr in Req_Headers:
    LUT_Dict[hdr]=np.array(LUT_Dict[hdr])

LUT_Dict['MEAN_DEPTH']=(LUT_Dict['HORZ_TOP']+LUT_Dict['HORZ_BOTTOM'])/200.

##################################################################################################

##################################################################################################
# Read in Shape file
SF=shapefile.Reader(shp_file)
field_names=[field[0] for field in SF.fields]
field_names.pop(0)

ShRcs=SF.shapeRecords()
##################################################################################################

##################################################################################################
# Loop round each shape and, calculate the composition for the shape, 
#   then attribute to each CHESS grid cell according to overlap fraction

# Set up dictionaries for stroing the data:
CHESS_SoilComp={ varname: np.zeros([nSD,ny,nx]) for varname in Comp_out_names  }

# Array for storing the Total Fraction accounted for in each cell
TOTAL_FRACTION=np.zeros([ny,nx])
MAX_SERIES=8

for iShRc in range(len(ShRcs)):
    ShRc=ShRcs[iShRc]
    
    # Get Soil Series ID and percent lists for shape
    series_num=1
    Soil_Series_ID=[]
    Soil_Series_PC=[]
    SIDindex=field_names.index('sercde'+str(series_num))
    SPCindex=field_names.index('pcent'+str(series_num))
    temp_seriesID=ShRc.record[SIDindex]
    temp_seriesPC=ShRc.record[SPCindex]
    while (temp_seriesID[:1]!='\x00')&(series_num<MAX_SERIES):
        Soil_Series_ID.append(int(temp_seriesID))
        Soil_Series_PC.append(float(temp_seriesPC)/100.)
        series_num+=1
        SIDindex=field_names.index('sercde'+str(series_num))
        SPCindex=field_names.index('pcent'+str(series_num))
        temp_seriesID=ShRc.record[SIDindex]
        temp_seriesPC=ShRc.record[SPCindex]
    
    nSS=len(Soil_Series_ID)
    if nSS==0:
        continue
        
    # Create profiles arrays for each output parameter for SEMI and CULT Land uses
    SEMI_profiles={ varname: np.zeros([nSD]) for varname in Comp_out_names }
    CULT_profiles={ varname: np.zeros([nSD]) for varname in Comp_out_names }
    PROFILE_PERCENT=0
    for iSS in range(nSS):
        full_LUT_index=np.where(LUT_Dict['SERIES_CODE']==Soil_Series_ID[iSS])[0]
        
        # Calculate Cultivated and Semi land-use independently, if and where available
        # Semi first
        if 'SEMI'in LUT_Dict['LAND_USE'][full_LUT_index]:
            LUT_index = full_LUT_index[LUT_Dict['LAND_USE'][full_LUT_index]=='SEMI']
            # loop over output soil depths
            for iSD in range(nSD):
                # find closest input SD
                SD_index=np.argmin(np.abs(  LUT_Dict['MEAN_DEPTH'][LUT_index] \
                                          - CH_mean_soil_depths[iSD]) )
                # INDEX of the closest soil depth for the series
                INDEX=LUT_index[SD_index]
                # Loop over output parameters
                for inname,outname in zip(Comp_in_names,Comp_out_names):
                    #for each param, at each soil depth, sum the weighted components
                    SEMI_profiles[outname][iSD]+=(LUT_Dict[inname][INDEX]*Soil_Series_PC[iSS])
                
            COPY_FROM_CULT=0
        else:
            COPY_FROM_CULT=1
    
        # Next 'CULT
        if 'CULT'in LUT_Dict['LAND_USE'][full_LUT_index]:
            LUT_index = full_LUT_index[LUT_Dict['LAND_USE'][full_LUT_index]=='CULT']
            # loop over output soil depths
            for iSD in range(nSD):
                # find closest input SD
                SD_index=np.argmin(np.abs(  LUT_Dict['MEAN_DEPTH'][LUT_index]   \
                                          - CH_mean_soil_depths[iSD]) )
                # INDEX of the closest soil depth for the series
                INDEX=LUT_index[SD_index]
                # Loop over output parameters
                for inname,outname in zip(Comp_in_names,Comp_out_names):
                    #for each param, at each soil depth, sum the weighted components
                    CULT_profiles[outname][iSD] += (LUT_Dict[inname][INDEX]*Soil_Series_PC[iSS])
                                
                sand=LUT_Dict[Comp_in_names[0]][INDEX]
                silt=LUT_Dict[Comp_in_names[1]][INDEX]
                clay=LUT_Dict[Comp_in_names[2]][INDEX]
                org_carb=LUT_Dict[Comp_in_names[3]][INDEX]
                    
            COPY_FROM_SEMI=0
        else:
            COPY_FROM_SEMI=1
    
        # Copy Profile Dictionaries from other land-use type where possible Where neccessary
        if (COPY_FROM_CULT==1)&(COPY_FROM_SEMI==0):
            if 'CULT'in LUT_Dict['LAND_USE'][full_LUT_index]:
                LUT_index = full_LUT_index[LUT_Dict['LAND_USE'][full_LUT_index]=='CULT']
                for iSD in range(nSD):
                    SD_index=np.argmin(np.abs(  LUT_Dict['MEAN_DEPTH'][LUT_index]   \
                                              - CH_mean_soil_depths[iSD] ) )
                    INDEX=LUT_index[SD_index]
                    for inname,outname in zip(Comp_in_names,Comp_out_names):
                        SEMI_profiles[outname][iSD] += (LUT_Dict[inname][INDEX]*Soil_Series_PC[iSS])
                     
                    sand=LUT_Dict[Comp_in_names[0]][INDEX]
                    silt=LUT_Dict[Comp_in_names[1]][INDEX]
                    clay=LUT_Dict[Comp_in_names[2]][INDEX]
                    org_carb=LUT_Dict[Comp_in_names[3]][INDEX]
                    
                COPY_FROM_CULT=0
            else:
                COPY_FROM_CULT=1
        elif (COPY_FROM_CULT==0)&(COPY_FROM_SEMI==1):
            if 'SEMI'in LUT_Dict['LAND_USE'][full_LUT_index]:
                LUT_index = full_LUT_index[LUT_Dict['LAND_USE'][full_LUT_index]=='SEMI']
                for iSD in range(nSD):
                    SD_index=np.argmin(np.abs(  LUT_Dict['MEAN_DEPTH'][LUT_index]  \
                                              - CH_mean_soil_depths[iSD]) )
                    INDEX=LUT_index[SD_index]
                    for inname,outname in zip(Comp_in_names,Comp_out_names):
                        CULT_profiles[outname][iSD] += (LUT_Dict[inname][INDEX]*Soil_Series_PC[iSS])
                    
                    sand=LUT_Dict[Comp_in_names[0]][INDEX]
                    silt=LUT_Dict[Comp_in_names[1]][INDEX]
                    clay=LUT_Dict[Comp_in_names[2]][INDEX]
                    org_carb=LUT_Dict[Comp_in_names[3]][INDEX]
                    
                COPY_FROM_CULT=0
            else:
                COPY_FROM_CULT=1
            
        if (COPY_FROM_CULT==0)|(COPY_FROM_SEMI==0):
            PROFILE_PERCENT+=Soil_Series_PC[iSS]
    
    # if Profile Percent = 0 there is no relevent information for the shape so move on
    # This will be closest neighbour gap filled later
    if PROFILE_PERCENT==0:
        continue
        
    # Now Normalise the Profiles based on the percentage complete
    for outname in Comp_out_names:
        CULT_profiles[outname]=CULT_profiles[outname]/PROFILE_PERCENT
        SEMI_profiles[outname]=SEMI_profiles[outname]/PROFILE_PERCENT
    
    # Calculate the fraction of each overlapping CHESS chess grid square
    limits=[int(bbox/1e3) for bbox in ShRc.shape.bbox]
    x_limits=[limits[0],limits[2]]
    y_limits=[limits[1],limits[3]]
    if min(y_limits)>ny:
        continue
    if min(x_limits)>nx:
        continue
    INpoly=Polygon(ShRc.shape.points)
    INpoly=INpoly.buffer(0)
    
    for i in range(x_limits[0],x_limits[1]+1):
        for j in range(y_limits[0],y_limits[1]+1):
            if not CHESS_cropFrac.mask[j,i]:
                # CHESS polygons are 1km squares with a LL corner = j,i
                OUTpoly = Square_Polygon(i,j,1e3)
                # Check for intersecting polygons
                if OUTpoly.intersects(INpoly):
                    # Get aarea of overlap
                    OVERLAP=OUTpoly.intersection(INpoly)
                    FRACTION=OVERLAP.area/1e6
                    TOTAL_FRACTION[j,i]+=FRACTION
                    for outname in Comp_out_names:
                        #for each param, at each soil depth, sum the weighted components
                        CHESS_SoilComp[outname][:,j,i]+=       \
                                ( CULT_profiles[outname] * CHESS_cropFrac[j,i] +     \
                                  SEMI_profiles[outname] * (1-CHESS_cropFrac[j,i]) ) \
                                * FRACTION

################################################################################################

#################################################################################################
# Filter any locations with a fraction cover less than 20%
ZEROindex=np.where(TOTAL_FRACTION<0.2)
TOTAL_FRACTION[ZEROindex]=0

# Mask missing values and Normalise to a TOTAL_Fraction of 1
TFindex=np.where(TOTAL_FRACTION>=0.2)
# Soil Composition
for outname in Comp_out_names:    
    CHESS_SoilComp[outname][:,ZEROindex[0],ZEROindex[1]]=0.
    CHESS_SoilComp[outname][:,TFindex[0],TFindex[1]]/= TOTAL_FRACTION[TFindex]
    CHESS_SoilComp[outname]=np.ma.masked_equal(CHESS_SoilComp[outname],0.)
    CHESS_SoilComp[outname].data[CHESS_SoilComp[outname].mask==True]=fill_value
    CHESS_SoilComp[outname].fill_value=fill_value
    
# Normalise sand/silt/clay fractions
TOTAL_SSC=(CHESS_SoilComp['sand']+CHESS_SoilComp['silt']+CHESS_SoilComp['clay'])/100.
for soil_comp in ['sand','silt','clay']:
    CHESS_SoilComp[soil_comp]/=TOTAL_SSC

#################################################################################################

#################################################################################################
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

# Loop around soil parameters and 
for var in CHESS_SoilComp:
    outvar=outf.createVariable(var,'float32',('z','y','x'),fill_value=fill_value)
    if var=='ph':
        outvar.units='pH'
    else:
        outvar.units='percent'
    outvar[:]=CHESS_SoilComp[var]
    

outf.Title='Soil composition and properties map for Scotland based on the Scotland Soils Database'
outf.owner='Edward Comyn-Platt'

outf.close()


LLinf.close()

