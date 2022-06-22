#!/bin/env python

from mapper import Mapper

MAPPING_FILE = "/users/global/ppha/swelter/um_amip/prep/amip_map.asc"

NC_METUM = 192
NR_METUM = 144
# NL_METUM = 10664

DLON_METUM = 1.8750
DLAT_METUM = 1.2500
LON0_METUM = 0.9375
LAT0_METUM = -89.375

NC_EURO = 23
NR_EURO = 20
# NP_EURO = NC_EURO*NR_EURO
# NL_EURO = 350

# DLON_EURO = 1.8750
# DLAT_EURO = 1.2500
LON0_EURO = 349.6875 - 360.0
LAT0_EURO = 35.6250


metum = Mapper(NC_METUM,NR_METUM,
               DLON_METUM,DLAT_METUM,LON0_METUM,LAT0_METUM,
               MAPPING_FILE)

import netCDF4 as nc
import matplotlib.pyplot as plt

FILE_NC = "/prj/SWELTER/data/amip/nc/amxvg/amxvga.pa19960101.nc"
ncid = nc.Dataset(FILE_NC,"r")
tstar = ncid.variables["tstar"][:,:,:].mean(axis=0)
ncid.close()

F = plt.figure()
A = F.add_subplot(111)

I = metum.plotMap(A,tstar)

plt.show()
