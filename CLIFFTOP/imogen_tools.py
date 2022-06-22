import numpy as np


def response(ncallyr, nyr):
    rs=np.zeros(ncallyr*nyr)
    for i in range(0,ncallyr*nyr):
        time_rs=np.float(i+1)/np.float(ncallyr)
        
        if (time_rs >= 1.0):
            rs[i] = 0.014819 + 0.70367*np.exp(-time_rs/0.70177) \
                    +0.24966*np.exp(-time_rs/2.3488) + 0.066485*np.exp(-time_rs/15.281) \
                    +0.038344*np.exp(-time_rs/65.359)+0.019439*np.exp(-time_rs/347.55)
        else:
            rs[i] = 1.0-2.2617*time_rs + 14.002*(time_rs**2)  \
                    -48.770*(time_rs**3) + 82.986*(time_rs**4) \
                    -67.527*(time_rs**5) + 21.037*(time_rs**6)
                    
    return rs


def ocean_co2_uptake(co2_atmos_ppmv, dt_ocean, fa_ocean):
    
    # co2_atmos_ppmv, time-series of atmospheric co2 concentration
    #                 Dims = [nyears]
    # dt_ocean, time-series of ocean surface delta temperature
    #                 Dims = [nyears]
    # fa_ocean, time-series of ocean memory
    #                 Dims = [nyear,nfarry]
    
    # This subroutine describes a simple uptake 
    # by the oceans of atmospheric CO2
    # The work is based upon the paper by Joos. 

    # The parameters chosen replicate those of the 3-D model, see 
    # Table 2 of Joos et al.
    
    # Parameters:
    k_g = 0.1306            # Gass exchange coefficient (/m2/yr)
    h = 40.0                # Mixed-layer depth (m)
    c = 1.722E17            # Units conversion parameter
    t_ocean_init = 289.28   # initial oean temperature
    t_mixed_init = t_ocean_init - 273.15
    ocean_area=3.627E14

    # Set up the number of calls per year
    ncallyr = 20
    nyr=co2_atmos_ppmv.shape[0]
    nfarray=fa_ocean.shape[1]
    # Ocean CO2 not called for last year of run 
    year_run=nyr-1
   
    co2_atmos_change_ppmv= np.append(0,co2_atmos_ppmv[1:]-co2_atmos_ppmv[:1])

    if(year_run*ncallyr > nfarray):
        print('Array size too small for FA_OCEAN')
    
    # Get the Green's function for use "under integral". (This is badly coded as doing each call)
    rs=response(ncallyr,nyr)

    d_ocean_atmos = np.zeros(nyr)

    # Assume initial mixed-layer temperature identical to initial temperature (units change needed)
    timestep_co2 = 1.0/np.float(ncallyr)
    co2_ocean_init = co2_atmos_ppmv[0]
    
    for iyr in range(nyr):
        for j in range(1,ncallyr+1):
            istep = (iyr*ncallyr)+j
            
            # Introduce linear correction in atmospheric CO2 concentration 
            # down to timescale 1/NCALLYR (and things centred mid-year).
            co2_atmos_short_timescale = co2_atmos_ppmv[iyr] + \
                    co2_atmos_change_ppmv[iyr] * ((np.float(j)/np.float(ncallyr))-0.5)
    
            dco2_atmos = co2_atmos_short_timescale-co2_atmos_ppmv[0]
            
            # Calculate perturbation in dissolved inorganic carbon (Eqn (3) of Joos)
            # for timestep indexed by istep.
            dco2_ocean_mol = 0.0      #Initialised for start of summation 
            if(istep >= 2):
                for i in range(1,istep):
                    dco2_ocean_mol = dco2_ocean_mol + \
                            (c/h)*fa_ocean[iyr,i-1]*rs[istep-i-1]*timestep_co2
                            
            # Relation between DCO2_OCEAN_MOL and DCO2_OCEAN:   Joos et al.
            # ie Convert from umol/kg to ppm
            dco2_ocean =(1.5568-(1.3993*1.0E-2*t_mixed_init))*dco2_ocean_mol \
                    + (7.4706-(0.20207*t_mixed_init))*1.0E-3*(dco2_ocean_mol**2) \
                    - (1.2748-(0.12015*t_mixed_init))*1.0E-5*(dco2_ocean_mol**3) \
                    + (2.4491-(0.12639*t_mixed_init))*1.0E-7*(dco2_ocean_mol**4) \
                    - (1.5468-(0.15326*t_mixed_init))*1.0E-10*(dco2_ocean_mol**5)
                    
            # Now incorporate correction suggested by Joos 
            co2_ocean = co2_ocean_init + dco2_ocean
            co2_ocean = co2_ocean * np.exp(0.0423*dt_ocean[iyr])
            
            fa_ocean[iyr,istep-1] = (k_g/ocean_area)*(co2_atmos_short_timescale-co2_ocean)
                    
            # Now calculate d_ocean_atmos (positive when flux is upwards)
            d_ocean_atmos[iyr] = d_ocean_atmos[iyr] - fa_ocean[iyr,istep-1]*ocean_area*timestep_co2
                    
    return d_ocean_atmos
    




