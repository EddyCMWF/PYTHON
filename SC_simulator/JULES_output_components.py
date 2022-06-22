 #!/usr/bin/env python

################################################################################
# 
# Program: SC_simulator.py
# 
# Python Script to simulate soil carbon estimate for a long-term equilibrium
#   condition, i.e. dCs/dt=0 .'. litterfall=microbial_soil_respiraiton       
#
# Author: Edward Comyn-Platt, edwcom@ceh.ac.uk
#
################################################################################
#
# 
#  
################################################################################
#
#
import pylab as plt
import numpy as np
import netCDF4 as nc
import JULES_SC_functions as J_SC_F
#

INTERVAL='SANITY_CHECK'
fill_value=-999.

#####################################################################################
# 1. Set dirs and filenames
###########################
print 'Setting Filenames and Directories'
# Input files:
JULES_output_dir  = '/users/eow/edwcom/CRUNCEP/n96/JULES_output/'
JULES_file = 'JULES_v4.3_TRIFFID_RsQ10_GLOBAL.monthly_mean.nc'
JULES_file_BS = 'JULES_v4.3_TRIFFID_RsQ10_GLOBAL_BigSpin.monthly_mean.nc'
#
# Output Dir
OUT_DIR = '/users/eow/edwcom/SC_simulator/output/'
out_tag = 'CRUNCEPn96_simulation_JULESoutput_BigSpin_'+INTERVAL
#

# 2.4 Read in BigSpin CS from JULES output
print 'Reading in JULES output data'
inf = nc.Dataset(JULES_output_dir+JULES_file_BS,'r')
lats=inf.variables['latitude'][:].squeeze()
cs_JULES_BS      = inf.variables['cs'][:].squeeze()
lit_C_JULES_BS   = inf.variables['lit_c'][:].squeeze()
resp_S_JULES_BS  = inf.variables['resp_s'][:].squeeze()
frac_JULES_BS  = inf.variables['frac'][:].squeeze()
c_veg_JULES_BS  = inf.variables['c_veg'][:].squeeze()
lai_JULES_BS  = inf.variables['lai'][:].squeeze()
inf.close()
#

cs_JULES_BS=np.ma.masked_invalid(cs_JULES_BS)
lit_C_JULES_BS=np.ma.masked_invalid(lit_C_JULES_BS)
resp_S_JULES_BS=np.ma.masked_invalid(resp_S_JULES_BS)
frac_JULES_BS=np.ma.masked_invalid(frac_JULES_BS)

cs_JULES_BS=cs_JULES_BS[-12:,:,:]
lit_C_JULES_BS=lit_C_JULES_BS[-12:,:,:]
resp_S_JULES_BS=resp_S_JULES_BS[-12:,:,:]*3600.*24.*360.
frac_JULES_BS=frac_JULES_BS[-12:,:,:]

lit_C_JULES_cum=np.sum(lit_C_JULES_BS,axis=0)
resp_S_JULES_cum=np.sum(resp_S_JULES_BS,axis=0)
lit_C_JULES_mean=np.mean(lit_C_JULES_BS,axis=0)
resp_S_JULES_mean=np.mean(resp_S_JULES_BS,axis=0)

frac_JULES_mean=np.mean(frac_JULES_BS,axis=0)

cs_JULES_cum =np.sum(cs_JULES_BS,axis=0)
cs_JULES_mean=np.mean(cs_JULES_BS,axis=0)


kappa_s_RothC  = np.array([ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10 ])

clay_content=0.23   # as in JULES
beta_r  = 1.0 / (4.0895 + 2.672*  np.exp(-0.0786 * 100.0 * clay_content))
alpha_dr= np.array([0.25,0.25,0.67,0.67,0.33])
f_dpm = alpha_dr/(1.+alpha_dr)

Lit_DPMpool=np.zeros_like(lit_C_JULES_cum[0,:],dtype='float64')
Lit_RPMpool=np.zeros_like(lit_C_JULES_cum[0,:],dtype='float64')
for i in range(5):
    Lit_DPMpool+= (f_dpm[i]*lit_C_JULES_mean[i,:]*frac_JULES_mean[i,:])
    Lit_RPMpool+= ((1-f_dpm[i])*lit_C_JULES_mean[i,:]*frac_JULES_mean[i,:])
    

C_dpm_equi = (Lit_DPMpool)/(resp_S_JULES_cum[0,:])
C_rpm_equi = (Lit_RPMpool)/(resp_S_JULES_cum[1,:])
C_bio_equi = (0.46*beta_r)*(np.sum(resp_S_JULES_cum,axis=0)/resp_S_JULES_cum[2,:])
C_hum_equi = (0.54*beta_r)*(np.sum(resp_S_JULES_cum,axis=0)/resp_S_JULES_cum[3,:])


n96_landindexfile='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

# Open and read the land index file
inf=nc.Dataset(n96_landindexfile,'r')
index=inf.variables['index_2D'][:]
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
inf.close()


cs_JULES_mean_2D=cs_JULES_mean[:,index]
for i in range(4):
    cs_JULES_mean_2D[i,:,:][np.where(index.mask==True)]=fill_value
cs_JULES_mean_2D=np.ma.masked_equal(cs_JULES_mean_2D,-999.)

resp_S_JULES_mean_2D=resp_S_JULES_mean[:,index]
for i in range(4):
    resp_S_JULES_mean_2D[i,:,:][np.where(index.mask==True)]=fill_value
resp_S_JULES_mean_2D=np.ma.masked_equal(resp_S_JULES_mean_2D,-999.)


Lit_DPMpool_2D=Lit_DPMpool[index]
Lit_DPMpool_2D[np.where(index.mask==True)]=fill_value
Lit_DPMpool_2D=np.ma.masked_equal(Lit_DPMpool_2D,-999.)

Lit_RPMpool_2D=Lit_RPMpool[index]
Lit_RPMpool_2D[np.where(index.mask==True)]=fill_value
Lit_RPMpool_2D=np.ma.masked_equal(Lit_RPMpool_2D,-999.)

C_dpm_equi_2D=C_dpm_equi[index]
C_rpm_equi_2D=C_rpm_equi[index]
C_bio_equi_2D=C_bio_equi[index]
C_hum_equi_2D=C_hum_equi[index]

C_dpm_equi_2D[:,:][np.where(index.mask==True)]=fill_value
C_rpm_equi_2D[:,:][np.where(index.mask==True)]=fill_value
C_bio_equi_2D[:,:][np.where(index.mask==True)]=fill_value
C_hum_equi_2D[:,:][np.where(index.mask==True)]=fill_value

C_dpm_equi_2D=np.ma.masked_equal(C_dpm_equi_2D,-999.)
C_rpm_equi_2D=np.ma.masked_equal(C_rpm_equi_2D,-999.)
C_bio_equi_2D=np.ma.masked_equal(C_bio_equi_2D,-999.)
C_hum_equi_2D=np.ma.masked_equal(C_hum_equi_2D,-999.)


plt.subplot(2,2,1)
plt.imshow(C_dpm_equi_2D,origin='bottom',vmin=-5,vmax=5,cmap='RdBu')
plt.colorbar(fraction=0.04)
plt.title('DPM')
plt.subplot(2,2,2)
plt.imshow(C_rpm_equi_2D,origin='bottom',vmin=-5,vmax=5,cmap='RdBu')
plt.colorbar(fraction=0.04)
plt.title('RPM')
plt.subplot(2,2,3)
plt.imshow(C_bio_equi_2D,origin='bottom',vmin=-5,vmax=5,cmap='RdBu')
plt.colorbar(fraction=0.04)
plt.title('Bio')
plt.subplot(2,2,4)
plt.imshow(C_hum_equi_2D,origin='bottom',vmin=0,vmax=2,cmap='RdBu')
plt.colorbar(fraction=0.04)
plt.title('Hummus')

plt.show()



plt.subplot(1,2,1)
plt.imshow(Lit_DPMpool_2D,origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Litter in DPM pool')
plt.subplot(1,2,2)
plt.imshow(Lit_RPMpool_2D,origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Litter in RPM pool')
plt.show()




plt.subplot(1,2,1)
plt.imshow(cs_JULES_mean_2D[0,:,:],origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Carbon in DPM pool')
plt.subplot(1,2,2)
plt.imshow(cs_JULES_mean_2D[1,:,:],origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Carbon in RPM pool')
plt.show()




plt.subplot(1,2,1)
plt.imshow(resp_S_JULES_mean_2D[0,:,:]*(3600.*24.*360.),origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Soil Resp. of DPM pool')
plt.subplot(1,2,2)
plt.imshow(resp_S_JULES_mean_2D[1,:,:]*(3600.*24.*360.),origin='bottom')#,vmin=0,vmax=20)
plt.colorbar(fraction=0.04)
plt.title('Soil Resp. of RPM pool')
plt.show()


