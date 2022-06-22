############################################
# Series of tools for jules post processing
############################################

#############################################
# nee_from_jules_output:
# Calculate the NEE from jules output:
#   (npp * frac) - sum(resp_s)
# or:
#   npp_gb - resp_s
# or:
#   any combination of the above
#
#############################################

def nee_from_jules_output(filename,\
                          l_npp_gb=False,l_resp_s_gb=False,\
                          frac=None):
    import netCDF4 as nc
    import numpy as np
    
    inf=nc.Dataset(filename,'r')
    tsteps=len(inf.dimensions['time'])
    
    if (l_npp_gb==False):
        # load frac if in file
        if (frac==None):
            frac=inf.variables['frac'][:].squeeze()
        # else populate from static input vector
        else:
            frac = np.array([ frac for i in range(tsteps) ])
        npp=inf.variables['npp'][:]. squeeze()
        npp_gb = np.sum( npp * frac, axis=1 )
    else:
        npp_gb = inf.variables['npp_gb'][:]. squeeze()
        
    if (l_resp_s_gb==False):
        resp_s=inf.variables['resp_s'][:]. squeeze()
        resp_s_gb = np.sum( resp_s, axis=1 )
    else:
        resp_s_gb = inf.variables['resp_s_gb'][:]. squeeze()
        
    

 

