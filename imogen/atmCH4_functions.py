
def ch4_ppbv_adjust(ch4_ppbv_tm1, atmos_gain_ch4, ch4_ppbv_init, ch4_ppbv_init_tm1, 
                    atmos_gain_ch4_init=(180e9)*(12.01/16.04),
                    tau_ch4_ref=8.4, ch4_ppbv_ref=1751.02,CONV=0.471 ):
    # Required Inputs:
    #  ch4_ppbv_tm1 = previous timestep ch4_ppbv
    #  atmos_gain_ch4 = -land_gain_ch4, global flux of ch4 to atmosphere from wetlands, 
    #                                   positive into atmos. Units = kg per year
    #  ch4_ppbv_init = ch4_ppbv of the initial projection
    #  ch4_ppbv_init_tm1 = previous timestep ch4_ppbv of the initial projection
    # 
    # Optional Inputs:
    #  atmos_gain_ch4_ref = global flux of ch4 to atmosphere from wetlands for the the initial projection
    #                        fixed at 180 TgCH4 for deviation from RCP/SSP scenarios
    #                        Units = kg per year
    #  tau_ch4_ref = lifetime of atmospherice ch4 in years for ch4_ppbv_ref
    #                  Units = years
    #  ch4_ppbv_ref = reference atmospheric ch4 concentration for tau_ch4_ref, see Climate Change 2001
    #  CONV = 0.471, conversion of GtC to ppmv
    
    # OLD CODE:
    #land_gain_ch4 = -fCH4_noperma[:,:,iyear]    # GtC per year *1e12
    #d_land_atmos_ch4 =-( land_gain_ch4 + fch4_ref_Gtc ) * CONV * 1e3
    # Calculate difference in global ch4 flux, converted to ppbv
    d_land_atmos_ch4 = ( atmos_gain_ch4 - atmos_gain_ch4_init )*1e-12 * CONV * 1e3 
    
    # atmos ch4 lifetime, see CC2001
    tau_ch4 = tau_ch4_ref * np.exp( 0.28 * np.log(ch4_ppbv/ch4_ppbv_ref) )
    ch4_ppbv_new = ch4_ppbv_tm1 + d_land_atmos_ch4           \
                - ((ch4_ppbv_tm1 - ch4_ppbv_init_tm1) / tau_ch4) \
                + (ch4_ppbv_init - ch4_ppbv_init_tm1) 
    return ch4_ppbv_new

