# module containing Radiativie forcing Calculations

import numpy as np

a1=-2.4e-7
b1=7.2e-4
c1=-2.1e-4

a2 = -8.0e-6
b2 =  4.2e-6
c2 = -4.9e-6

a3 = -1.3e-6
b3 = -8.2e-6

def etminan_CO2(co2_ppm,n2o_ppb,co2_ppm_0=285.,n2o_ppb_0=272.):
    co2_diff=co2_ppm-co2_ppm_0
    N2Obar = (n2o_ppb+n2o_ppb_0)/2.
    delQ_co2 = (  (a1 * (co2_diff**2))  \
               + (b1 * co2_diff)        \
               + (c1 * N2Obar) + 5.36 ) \
             * np.log(co2_ppm/co2_ppm_0)

    return delQ_co2


def etminan_CO2_inverse(delQ,co2_init_ppm,n2o_ppb,
                    co2_ppm_0=285., n2o_ppb_0=272., maxitercnt=1000):
    Nbar = (n2o_ppb+n2o_ppb_0)*0.5
    co2_iter=999.
    co2=np.copy(co2_init_ppm)
    itercnt=0
    while ( (np.abs(co2_iter-co2)>0.001) & (itercnt<=maxitercnt) ):
        co2_iter=co2.copy()
        co2_diff = co2_iter-co2_ppm_0
        denom = ( (a1 * (co2_diff**2.))     
                + (b1 * co2_diff )           
                + (c1 * Nbar) + 5.36 )
        co2 = co2_ppm_0 * np.exp( delQ/denom )
        itercnt+=1
    return co2

def etminan_CO2_inverse_series(delQ,n2o_ppb,
                    co2_ppm_0=285., n2o_ppb_0=272., maxitercnt=1000):
    ntsteps = delQ.shape[0]
    co2_ppm_out = np.zeros(ntsteps)
    co2_init = co2_ppm_0
    for itstep in range(ntsteps):
        co2_ppm_out[itstep] = etminan_CO2_inverse(
                delQ[itstep], co2_init, n2o_ppb[itstep],
                co2_ppm_0=co2_ppm_0,n2o_ppb_0=n2o_ppb_0, maxitercnt=maxitercnt)
        co2_init = co2_ppm_out[itstep]*1.0
    return co2_ppm_out


def etminan_CH4(ch4_ppb,n2o_ppb,ch4_ppb_0=720.,n2o_ppb_0=272.):
    N2Obar = (n2o_ppb+n2o_ppb_0)*0.5
    CH4bar = (ch4_ppb+ch4_ppb_0)*0.5
    delQ_ch4 = ( (a3*CH4bar)  \
                +(b3*N2Obar)  \
                + 0.043 )                          \
              * ( np.sqrt(ch4_ppb)-np.sqrt(ch4_ppb_0) )

    return delQ_ch4

def etminan_N2O(n2o_ppb,co2_ppm,ch4_ppb,n2o_ppb_0=272.,co2_ppm_0=285.,ch4_ppb_0=720.):
    CO2bar = (co2_ppm+co2_ppm_0)/2.
    N2Obar = (n2o_ppb+n2o_ppb_0)/2.
    CH4bar = (ch4_ppb+ch4_ppb_0)/2.
    delQ_n2o = ( (a2*CO2bar)  \
                +(b2*N2Obar)  \
                +(c2*CH4bar)  \
                + 0.117 )                          \
              * ( np.sqrt(n2o_ppb)-np.sqrt(n2o_ppb_0) )

    return delQ_n2o

def collins_CH4StratoOzone(ch4_ppb,ch4_ppb_0=720.):
    delQ_Ozone_ch4=2.36e-4*(ch4_ppb-ch4_ppb_0)
    return delQ_Ozone_ch4

def SSP_CH4StratoOzone(ch4_ppb,ch4_ppb_0=720.):
    delQ_Ozone_ch4=0.043*5*(np.log(ch4_ppb)-np.log(ch4_ppb_0))
    return delQ_Ozone_ch4

