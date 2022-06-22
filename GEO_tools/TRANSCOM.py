########################################
# Tool box for TRANSCOM
########################################

def region_names():
    return ['North American Boreal',     \
            'North American Temperate',  \
            'South American Tropical',   \
            'South American Temperate',  \
            'Northern Africa',           \
            'Southern Africa',           \
            'Asia Boreal',               \
            'Asia Temperate',            \
            'Asia Tropical',             \
            'Australia',                 \
            'Europe',                    \
            'North Pacific Temperate',   \
            'West Pacific Tropics',      \
            'East Pacific Tropics',      \
            'South Pacific Temperate',   \
            'Northern Ocean',            \
            'North Atlantic Temperate',  \
            'Atlantic Tropics',          \
            'South Atlantic Temperate',  \
            'Southern Ocean',            \
            'Indian Tropics',            \
            'South Indian Temperate',    \
            'Antarctic and Greenland',   \
            ]


def region_Snames():
    return ['N. Amer. Bor.',  \
            'N. Amer. Temp.', \
            'S. Amer. Trop.', \
            'S. Amer. Temp.', \
            'N. Africa',      \
            'S. Africa',      \
            'Asia Bor.',      \
            'Asia Temp.',     \
            'Asia Trop.',     \
            'Australia',      \
            'Europe',         \
            'N. Pac. Temp.',  \
            'W. Pac. Trop.',  \
            'E. Pac. Trop.',  \
            'S. Pac. Temp.',  \
            'Nthn Ocean',     \
            'N. Atl. Temp.',  \
            'Atl. Trop.',     \
            'S. Atl. Temp.',  \
            'Sthn Ocean',     \
            'Ind. Trop.',     \
            'S.Ind. Temp.',   \
            'Ant. and GL',    \
            ]


def land_region_names():
    return region_names()[:11]


def land_region_Snames():
    return region_Snames()[:11]

def load_data(Tncfile='/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5.nc', \
        varname='transcom_regions', SHIFT=True):
    import netCDF4 as nc
    import numpy as np
    inf=nc.Dataset(Tncfile,'r')
    Tdata=inf.variables[varname][:]
    inf.close()
    if SHIFT==True:
        Tdata=np.append(Tdata[:,360:],Tdata[:,:360],axis=1)
    return Tdata

#def get_transcomindexes(Tncfile='/prj/ALANIS/UM_Modelling/TRANSCOM_Regions_0.5.nc', \
#                            varname='transcom_regions')
#    load





