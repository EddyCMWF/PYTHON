def ocean_co2(iyr, co2_atmos_ppmv, co2_atmos_ppmv_init, dt_ocean,
   fa_ocean,ocean_area,co2_atmos_change_ppmv,nyr,t_ocean_init,
   nfarray):

   import response
   import numpy as np
   import sys, os

# This subroutine describes a simple uptake 
# by the oceans of atmospheric CO2
# The work is based upon the paper by Joos. 

# The parameters chosen replicate those of the 3-D model, see 
# Table 2 of Joos et al.

# Ocean CO2 not called for last year of run 
   year_run=nyr-1

# Parameters from the Joos model
   k_g = 0.1306            # Gass exchange coefficient (/m2/yr)
   h = 40.0                # Mixed-layer depth (m)
   c = 1.722E17            # Units conversion parameter
   t_mixed_init = t_ocean_init - 273.15

# Set up the number of calls per year
   ncallyr = 20

   if(year_run*ncallyr > nfarray):
      print 'Array size too small for FA_OCEAN'

# Get the Green's function for use "under integral". (This is badly coded as doing each call)
   rs=response.response(ncallyr,nyr)

   d_ocean_atmos = 0.0
# Assume initial mixed-layer temperature identical to initial temperature (units change needed)
   timestep_co2 = 1.0/np.float(ncallyr)

   for j in range(1,ncallyr+1):
      istep = (iyr*ncallyr)+j
      co2_ocean_init = co2_atmos_ppmv_init

# Introduce linear correction in atmospheric CO2 concentration 
# down to timescale 1/NCALLYR (and things centred mid-year).
      co2_atmos_short_timescale = co2_atmos_ppmv + \
           co2_atmos_change_ppmv * ((np.float(j)/np.float(ncallyr))-0.5) 

      dco2_atmos = co2_atmos_short_timescale-co2_atmos_ppmv_init

# Calculate perturbation in dissolved inorganic carbon (Eqn (3) of Joos)
# for timestep indexed by istep.

      dco2_ocean_mol = 0.0		#Initialised for start of summation 
      if(istep >= 2):
         for i in range(1,istep):
            dco2_ocean_mol = dco2_ocean_mol + \
           (c/h)*fa_ocean[i-1]*rs[istep-i-1]*timestep_co2

# Relation between DCO2_OCEAN_MOL and DCO2_OCEAN:   Joos et al.
# ie Convert from umol/kg to ppm
      dco2_ocean = \
         (1.5568-(1.3993*1.0E-2*t_mixed_init))*dco2_ocean_mol \
       + (7.4706-(0.20207*t_mixed_init))*1.0E-3*(dco2_ocean_mol**2) \
       - (1.2748-(0.12015*t_mixed_init))*1.0E-5*(dco2_ocean_mol**3) \
       + (2.4491-(0.12639*t_mixed_init))*1.0E-7*(dco2_ocean_mol**4) \
       - (1.5468-(0.15326*t_mixed_init))*1.0E-10*(dco2_ocean_mol**5) 

# Now incorporate correction suggested by Joos 
      co2_ocean = co2_ocean_init + dco2_ocean
      co2_ocean = co2_ocean * np.exp(0.0423*dt_ocean)

      fa_ocean[istep-1] = \
         (k_g/ocean_area)*(co2_atmos_short_timescale-co2_ocean)

# Now calculate d_ocean_atmos (positive when flux is upwards)
      d_ocean_atmos = d_ocean_atmos -  \
         fa_ocean[istep-1]*ocean_area*timestep_co2

# Temporary Code
#   print 'iyr = ',iyr 
#   print 'co2_atmos_ppmv = ',co2_atmos_ppmv
#   print 'co2_atmos_ppmv_init = ',co2_atmos_ppmv_init
#   print 'dt_ocean = ',dt_ocean
#   print 'fa_ocean 1,2,3 = ',fa_ocean[0:3]
#   print 'ocean_area = ', ocean_area
#   print 'co2_atmos_change_ppmv = ', co2_atmos_change_ppmv
#   print 'year_run = ',year_run 
#   print 't_mixed_init = ',t_mixed_init
#   print 'nfarray = ',nfarray
#   print 'd_ocean_atmos = ',d_ocean_atmos 

   return d_ocean_atmos #, fa_ocean 
