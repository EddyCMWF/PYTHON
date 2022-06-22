#!//bin/env python

import netCDF4 as nc
import sys

runid=sys.argv[1]

START_YEAR=1990
END_YEAR=2015
Nyears=END_YEAR-START_YEAR

IN_DIR='/prj/wetlands_africa/jules/JASMIN/BVOCs/JULES_OUTPUT/WFDEI_global/'
OUT_DIR=IN_DIR+'yearly/'


infile=IN_DIR+runid+'.monthly.nc'
print(infile)
inf=nc.Dataset(infile,'r')

for iyear in range(Nyears+1):
    year=str(START_YEAR+iyear)
    outfile=OUT_DIR+runid+'.monthly.'+year+'.nc'
    print(outfile)
    outf=nc.Dataset(outfile,'w')
    
    for dim in inf.dimensions:
        if str(dim)!='time':
            outf.createDimension(str(dim),len(inf.dimensions[dim]))
        else:
            outf.createDimension('time',12)

    for var in inf.variables:
        print(var)
        invar=inf.variables[var]
        outvar=outf.createVariable(str(var),'float32',invar.dimensions)
        for att in invar.ncattrs():
            outvar.setncattr(att,invar.getncattr(att))

        if str(var) == 'time':
            outvar[:]=invar[iyear*12:(iyear*12)+12]
        elif (str(var)=='latitude')|(str(var)=='longitude'):
            outvar[:]=invar[:]
        else:
            outvar[:]=invar[iyear*12:(iyear*12)+12]
    
    outf.close()

inf.close()






