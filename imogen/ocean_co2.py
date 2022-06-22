from imogen import response
import numpy as np

def ocean_co2(iyr, co2_atmos_ppmv, co2_atmos_ppmv_init, dt_ocean,\
   fa_ocean,ocean_area,co2_atmos_change_ppmv,nyr,t_ocean_init   ):


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


# Get the Green's function for use "under integral". (This is badly coded as doing each call)
   rs=response.response(ncallyr,nyr)

   d_ocean_atmos = 0.0
# Assume initial mixed-layer temperature identical to initial temperature (units change needed)
   timestep_co2 = 1.0/np.float(ncallyr)

   co2_ocean_init = co2_atmos_ppmv_init
   for j in range(1,ncallyr+1):
      istep = (iyr*ncallyr)+j

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


   return d_ocean_atmos 

def ocean_co2_ECP(iyr, co2_atmos_ppmv, co2_ocean_ppmv, 
                    co2_atmos_ppmv_init, dt_ocean, dco2_atmos_ppmv,
                    fa_ocean,rs,ncallyr=20,
                    t_ocean_init=289.28,ocean_area=3.627e14 ):

# This subroutine describes a simple uptake 
# by the oceans of atmospheric CO2
# The work is based upon the paper by Joos. 

# The parameters chosen replicate those of the 3-D model, see 
# Table 2 of Joos et al.

# Parameters from the Joos model
   k_g = 0.1306            # Gass exchange coefficient (/m2/yr)
   h = 40.0                # Mixed-layer depth (m)
   c = 1.722E17            # Units conversion parameter
   t_mixed_init = t_ocean_init - 273.15

   d_ocean_atmos = 0.0
# Assume initial mixed-layer temperature identical to initial temperature (units change needed)
   timestep_co2 = 1.0/np.float(ncallyr)

   #co2_ocean_init = co2_ocean_ppmv
   co2_ocean_init =  co2_atmos_ppmv_init
   for j in range(1,ncallyr+1):
      istep = (iyr*ncallyr)+j

# Introduce linear correction in atmospheric CO2 concentration 
# down to timescale 1/NCALLYR (and things centred mid-year).
      co2_atmos_short_timescale = co2_atmos_ppmv + \
           dco2_atmos_ppmv * ((np.float(j)/np.float(ncallyr))-0.5) 

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


   return d_ocean_atmos  ,fa_ocean 

def ocean_co2_parallel(co2_atmos_ppmv, dt_ocean,\
                        ocean_area=3.627e14,t_ocean_init=289.28   ):


    # This subroutine describes a simple uptake 
    # by the oceans of atmospheric CO2
    # The work is based upon the paper by Joos. 
    
    # Calculates for a number of gcms in parallel.
    # Input array dimensions must be [ YEARs, GCM ]
    
    # Variables and options
    # Parameters from the Joos model
    k_g = 0.1306            # Gass exchange coefficient (/m2/yr)
    h = 40.0                # Mixed-layer depth (m)
    c = 1.722E17            # Units conversion parameter
    t_mixed_init = t_ocean_init - 273.15
    
    # Set up the number of calls per year
    ncallyr = 20
    nfarray = 20000
    timestep_co2 = 1.0/np.float(ncallyr)
    
    # Get values and arrays from input data:
    nGCM = co2_atmos_ppmv.shape[1]
    nyr  = co2_atmos_ppmv.shape[0]
    co2_atmos_ppmv_init = np.copy(co2_atmos_ppmv[0,:]) 
    co2_ocean_init = np.copy(co2_atmos_ppmv_init)
    d_ocean_atmos = np.zeros_like(co2_atmos_ppmv)
    co2_atmos_change_ppmv = np.zeros_like(co2_atmos_ppmv)
    co2_atmos_change_ppmv[1:,:] = co2_atmos_ppmv[1:,:]-co2_atmos_ppmv[:-1,:]
    
    fa_ocean=np.zeros([nfarray,nGCM])
    # Ocean CO2 not called for last year of run 
    year_run=nyr-1

    # Get the Green's function for use "under integral". 
    rs=response.response(ncallyr,nyr)
    
    for iyr in range(nyr):
        for j in range(1,ncallyr+1):
            istep = (iyr*ncallyr)+j
            
            # Introduce linear correction in atmospheric CO2 concentration 
            # down to timescale 1/NCALLYR (and things centred mid-year).
            co2_atmos_short_timescale = co2_atmos_ppmv[iyr,:] + \
                    co2_atmos_change_ppmv[iyr,:] * ((np.float(j)/np.float(ncallyr))-0.5) 
            
            ##dco2_atmos = co2_atmos_short_timescale-co2_atmos_ppmv_init
            
            # Calculate perturbation in dissolved inorganic carbon (Eqn (3) of Joos)
            # for timestep indexed by istep.
            dco2_ocean_mol = np.zeros(nGCM)
            if(istep >= 2):
                for i in range(1,istep):
                    dco2_ocean_mol = dco2_ocean_mol + \
                            (c/h)*fa_ocean[i-1,:]*rs[istep-i-1]*timestep_co2
                            
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
            co2_ocean = co2_ocean * np.exp(0.0423*dt_ocean[iyr,:])
            
            fa_ocean[istep-1,:] = \
                    (k_g/ocean_area)*(co2_atmos_short_timescale-co2_ocean)
                    
            # Now calculate d_ocean_atmos (positive when flux is upwards)
            d_ocean_atmos[iyr,:] = d_ocean_atmos[iyr,:] -  \
                    fa_ocean[istep-1,:]*ocean_area*timestep_co2
                                

    return d_ocean_atmos 

