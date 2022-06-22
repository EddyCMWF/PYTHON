# ECP netCDF tools

def copy_dimension(inf,outf,dimension):
    outf.createDimension(dimension,len(inf.dimensions[dimension]))

def copy_all_dimensions(inf,outf):
    for dim in inf.dimensions:
        copy_dimension(inf,outf,str(dim))

def copy_attribute(invar,outvar,attr):
    outvar.setncattr( str(attr),            \
                      invar.getncattr(attr) )

def copy_all_attributes(invar,outvar):
    for attr in invar.ncattrs():
        # Don't copy '_FillValue', 
        # this must be defined when variable is created
        if str(attr)!='_FillValue':
            copy_attribute(invar,outvar,str(attr))

def copy_variable(inf,outf,var,               \
        dtype=None,dimensions=None,fill_value=None):
    var=str(var)
    invar=inf.variables[var]
    if dtype==None:
        dtype=invar.dtype
    if dimensions==None:
        dimensions=invar.dimensions
    if (fill_value==None)&('_FillVallue' in invar.ncattrs()):
        fill_value=invar._FillValue
    ##print var, dtype, dimensions, fill_value
    if fill_value==None:
        outvar=outf.createVariable(var,dtype,dimensions)
    else:
        outvar=outf.createVariable(var,dtype,dimensions,fill_value=fill_value)

    copy_all_attributes(invar,outvar)
    outvar[:]=inf.variables[var][:]



