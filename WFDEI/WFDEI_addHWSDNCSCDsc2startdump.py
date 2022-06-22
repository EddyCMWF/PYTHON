#!/usr/bin/python2.7
#
# Program: PLOT_STARTDUMPS.py
# Purpose: Plot data from start dumps to check quality
# 
# Edward Comyn-Platt, 2015
#
######################################################
import numpy as np
import netCDF4 as nc

in_SD_oti='/work/scratch/ecomynplatt/Wetland_Diagnostics_Output/GLOBAL_0.5degree_WFD-EI_DougDiag/spindumps/JULES_v42_WFD-EI_MPI_global_DougDiag.dump.spin5.19800101.0.nc'
in_SD_nti='/work/scratch/ecomynplatt/Wetland_Diagnostics_Output/GLOBAL_0.5degree_WFD-EI_DougDiag/spindumps/MPI_intel_WFD-EI_v4_global_DougDiag_newtopo.dump.spin5.19800101.0.nc'

SC_file='/group_workspaces/jasmin/jules/data/WFD-EI-Forcing/ancils/qrparm.soil_merge_HWSD_NCSCD_cont_cosbyWFDEI.nc'
SC_param_name='field1397'

outfile_oti='/home/users/ecomynplatt/JULES/MPI_tstep_v4.2_global_WFD_EI/dumps/WFDEI_startdump_withHWSDNCSCD-SC_ECP20150909.nc'
outfile_nti='/home/users/ecomynplatt/JULES/MPI_tstep_v4.2_global_WFD_EI/dumps/WFDEI_newtopo_startdump_withHWSDNCSCD-SC_ECP20150909.nc'

#get SC data to substitute in
inf=nc.Dataset(SC_file)
new_SC=inf.variables[SC_param_name][:].flatten()
inf.close()

# first oti
inf=nc.Dataset(in_SD_oti,'r')
outf=nc.Dataset(outfile_oti,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    #for ncattr in inf.variables[str(var)].ncattrs():
    #    outvar.setncattr(str(ncattr),inf.variables[str(var)].getncattr(str(ncattr)))
    #
    if not (str(var) in ['cs']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        print 'here'
        outvar[:]=new_SC

outf.close()
inf.close()


# first oti
inf=nc.Dataset(in_SD_nti,'r')
outf=nc.Dataset(outfile_nti,'w')

for dim in inf.dimensions:
    outf.createDimension(str(dim),len(inf.dimensions[str(dim)]))

for var in inf.variables:
    in_dims=list(inf.variables[str(var)].dimensions)
    in_dtype=inf.variables[str(var)].dtype
    #
    outvar=outf.createVariable(str(var),in_dtype,in_dims)
    #
    #for ncattr in inf.variables[str(var)].ncattrs():
    #    outvar.setncattr(str(ncattr),inf.variables[str(var)].getncattr(str(ncattr)))
    #
    if not (str(var) in ['cs']):
        outvar[:]=inf.variables[str(var)][:]
    else:
        outvar[:]=new_SC

outf.close()
inf.close()





