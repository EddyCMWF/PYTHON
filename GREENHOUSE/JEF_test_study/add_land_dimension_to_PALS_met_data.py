#! /bin/python
#
# Script to add a Dummy land dimension to the PALS met data
#
import netCDF4 as nc
import glob

IN_DIR = '/data/grp/fluxdata/PALS_sites/'

OUT_DIR = '/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/PALS_comparison/met_data/'

files=glob.glob(IN_DIR+'*Fluxnet.1.4_met.nc')


for file in files:
    
    inf=nc.Dataset(file,'r')
    outf=nc.Dataset(OUT_DIR+file.split('/')[-1],'w')
    
    for dim in inf.dimensions:
        outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))
        
    outf.createDimension('land',1)
    
    for var in inf.variables:
        if str(var) in ['lat','lon','alt']:
            outvar=outf.createVariable(str(var),\
                                       inf.variables[var][:].dtype, \
                                       ('land',)                    )
        else:
            outvar=outf.createVariable(str(var), \
                                       inf.variables[var].dtype, \
                                       ('time','land') )

        for attr in inf.variables[var].ncattrs():
            outvar.setncattr( str(attr), \
                             inf.variables[var].getncattr(attr) )
        outvar[:]=inf.variables[var][:]

    outf.close()
    inf.close()

