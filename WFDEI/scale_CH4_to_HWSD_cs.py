#!/usr/bin/python
#
# rescale JULES fch4_wetl from the Zinke soil carbon to the 
#  merged HWSD and NCSCD soil carbon
# 
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# April 2015
#
# Contains
import numpy as np
import netCDF4 as nc
import netcdftime as nctime
#
fill_value=-999.0
# Files and directories:
WFDEI_dir='/users/eow/edwcom/WFD_EI/'

Zinke_SC_file='qrparm.soil_HWSD_class3_van_genuchtenNew.nc'
HWSD_NCSCD_SC_file='qrparm.soil_merge_HWSD_NCSCD_cont_cosbyWFDEI.nc'
HWSD_SC_file='soil_HWSD_soilcarbon_landpoints_0p5.safe.nc'


# Names and parameters:
Zinke_SC_name='field1397'
HWSD_SC_name='field1397'





# Open SC files and calculate ratio:
# Zincke:
inf=nc.Dataset(WFDEI_dir+Zinke_SC_file,'r')
Zinke_SC=inf.variables[Zinke_SC_name][:].squeeze()
inf.close()
# HWSD:
inf=nc.Dataset(WFDEI_dir+HWSD_SC_file,'r')
HWSD_SC=inf.variables[HWSD_SC_name][:].squeeze()
inf.close()
# HWSD/NCSCD:
inf=nc.Dataset(WFDEI_dir+HWSD_NCSCD_SC_file,'r')
HWSD_NCSCD_SC=inf.variables[HWSD_SC_name][:].squeeze()
inf.close()


C_ratio=np.zeros_like(Zinke_SC)

C_ratio=HWSD_SC/Zinke_SC
C_ratio_merged=HWSD_NCSCD_SC/Zinke_SC
































quit()

WFD_gridfile='/users/eow/edwcom/WFD_EI/wfdei-land-mask.nc'
inf_grid = nc.Dataset(WFD_gridfile,'r')
LAND_FRAC= inf_grid.variables['land_fraction'][:]
LAND_IND = inf_grid.variables['land_index'][:]
inf_grid.close()


Zinke_2D   = Zinke_SC[LAND_IND-1]
HWSD_2D    = HWSD_SC[LAND_IND-1]
C_ratio_2D = C_ratio[LAND_IND-1]

Zinke_2D[np.where(LAND_IND.mask==True)]=fill_value
HWSD_2D[np.where(LAND_IND.mask==True)]=fill_value
C_ratio_2D[np.where(LAND_IND.mask==True)]=fill_value


Zinke_2D   = np.ma.masked_equal(Zinke_2D ,fill_value)
HWSD_2D    = np.ma.masked_equal(HWSD_2D ,fill_value)
C_ratio_2D = np.ma.masked_equal(C_ratio_2D ,fill_value)



plt.subplot(1,3,1)
plt.imshow(Zinke_2D,origin='bottom',vmin=0,vmax=50)
plt.colorbar(fraction=0.03)
plt.title('Zinke Soil Carbon')
plt.subplot(1,3,2)
plt.imshow(HWSD_2D,origin='bottom',vmin=0,vmax=50)
plt.colorbar(fraction=0.03)
plt.title('HWSD Soil Carbon')
plt.subplot(1,3,3)
plt.imshow(C_ratio_2D,origin='bottom',vmax=5)
plt.colorbar(fraction=0.03)
plt.title('SC Ratio')

plt.show()











