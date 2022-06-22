#!/bin/python
import netCDF4 as nc
import numpy as np

in_site_dump='/users/eow/edwcom/GREENHOUSE/PALS_comparison/Jvn4.3.1-E-F_fullsetup/output/Loobos.dump.19970101.0.nc'

in_CHESS_dump='/users/eow/edwcom/CHESS/dumps/J4.3_CHESS_ConstSoil_spin2_startdump.nc'

landindex_file='/users/eow/edwcom/CHESS/chess_jules_land_index.nc'
landcover_file='/users/eow/edwcom/CHESS/chess_landcover_2000.nc'

out_CHESS_dump='/users/eow/edwcom/CHESS/dumps/J4.3-E-F_CHESS_Loobos_startdump.nc'

grinf=nc.Dataset(landindex_file,'r')
grindex=grinf.variables['index_to1D'][:]
grinf.close()

LCinf=nc.Dataset(landcover_file,'r')
frac=LCinf.variables['frac'][:]
LCinf.close()

infCH=nc.Dataset(in_CHESS_dump,'r')
#lats=infCH.variables['latitude'][:]
#lons=infCH.variables['longitude'][:]
land_dim_len=len(infCH.dimensions['land'])
tile_dim_len=len(infCH.dimensions['tile'])

#infCH.close()


outf=nc.Dataset(out_CHESS_dump,'w')

inf=nc.Dataset(in_site_dump,'r')

# Copy dimensions, change land dim to new land dim
for dim in inf.dimensions:
    if str(dim)=='land':
        outf.createDimension(str(dim),land_dim_len)
    elif str(dim)=='tile':
        outf.createDimension(str(dim),tile_dim_len)
    elif str(dim)=='type':
        outf.createDimension(str(dim),tile_dim_len)
    else:
        outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    if str(var) in ['b','sathh','satcon','sm_sat','sm_crit','sm_wilt','hcap','hcon','albsoil']:
        continue

    print var
    outvar=outf.createVariable(str(var),'float32',inf.variables[str(var)].dimensions)
    if (str(var)=='latitude') | (str(var)=='longitude') |  \
       ('land' not in inf.variables[str(var)].dimensions):
        outvar[:]=infCH.variables[str(var)][:]
    elif ('tile' in inf.variables[str(var)].dimensions):
        outdata=np.array( [ inf.variables[str(var)][:8,:].squeeze() \
                            for i in range(land_dim_len) ] )
        outvar[:]=outdata
    else:
        outdata=np.array( [ inf.variables[str(var)][:].squeeze() \
                              for i in range(land_dim_len) ] )
        outvar[:]=outdata

inf.close()
outf.close()
infCH.close()


