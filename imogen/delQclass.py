# module containing Radiativie forcing Calculations

import numpy as np

class delQ:
    def __init__(self):
        # Define etminan parameters
        self.a1=-2.4e-7
        self.b1=7.2e-4
        self.c1=-2.1e-4
        self.a2 = -8.0e-6
        self.b2 =  4.2e-6
        self.c2 = -4.9e-6
        self.a3 = -1.3e-6
        self.b3 = -8.2e-6
        # Define IPCC parameters
        self.q2co2 = 3.74
        
        #Addtional Parameters
        self.Collins_O3_c = 2.36e-4

        # Pre industrial parameters, these are modifiable can be bypassed 
        #  when calling the delQ equations
        self.co2_ppm_0=285.
        self.n2o_ppb_0=272.
        self.ch4_ppb_0=720.

    def etminan_CO2(self, co2_ppm, n2o_ppb, co2_ppm_0=None, n2o_ppb_0=None):
        if co2_ppm_0==None: co2_ppm_0=self.co2_ppm_0 
        if n2p_ppb_0==None: n2o_ppb_0=self.n2o_ppb_0
        
        co2_diff = co2_ppm-co2_ppm_0
        N2Obar = (n2o_ppb+n2o_ppb_0)/2.
        delQ_co2 =(  (self.a1 * (co2_diff**2))  
                   + (self.b1 * co2_diff)        
                   + (self.c1 * N2Obar) + 5.36 ) * np.log(co2_ppm/co2_ppm_0)
        return delQ_co2

    def etminan_CO2_inverse(self, delQ,co2_init_ppm,n2o_ppb,
                        co2_ppm_0=None, n2o_ppb_0=None, maxitercnt=1000):
        if co2_ppm_0==None: co2_ppm_0=self.co2_ppm_0 
        if n2p_ppb_0==None: n2o_ppb_0=self.n2o_ppb_0
        Nbar = (n2o_ppb+n2o_ppb_0)*0.5
        co2_iter=999.
        co2=np.copy(co2_init_ppm)
        itercnt=0
        while ( (np.abs(co2_iter-co2)>0.001) & (itercnt<=maxitercnt) ):
            co2_iter=co2.copy()
            co2_diff = co2_iter-co2_ppm_0
            denom = ( (self.a1 * (co2_diff**2.))
                    + (self.b1 *  co2_diff )           
                    + (self.c1 *  Nbar) + 5.36 )
            co2 = co2_ppm_0 * np.exp( delQ/denom )
            itercnt+=1
        return co2

    def etminan_CH4(self, ch4_ppb,n2o_ppb,ch4_ppb_0=None,n2o_ppb_0=None):
        if ch4_ppb_0==None: ch4_ppb_0=self.ch4_ppb_0 
        if n2p_ppb_0==None: n2o_ppb_0=self.n2o_ppb_0
        N2Obar = (n2o_ppb+n2o_ppb_0)*0.5
        CH4bar = (ch4_ppb+ch4_ppb_0)*0.5
        delQ_ch4 = ( (self.a3*CH4bar)  
                    +(self.b3*N2Obar)  
                    + 0.043 ) * ( np.sqrt(ch4_ppb)-np.sqrt(ch4_ppb_0) )
        return delQ_ch4

    def etminan_N2O(n2o_ppb, co2_ppm, ch4_ppb,
                    n2o_ppb_0=None,co2_ppm_0=None,ch4_ppb_0=None):
        if co2_ppm_0==None: co2_ppm_0=self.co2_ppm_0 
        if ch4_ppb_0==None: ch4_ppb_0=self.ch4_ppb_0 
        if n2p_ppb_0==None: n2o_ppb_0=self.n2o_ppb_0

        CO2bar = (co2_ppm+co2_ppm_0)/2.
        N2Obar = (n2o_ppb+n2o_ppb_0)/2.
        CH4bar = (ch4_ppb+ch4_ppb_0)/2.
        delQ_n2o = ( (self.a2*CO2bar)  
                    +(self.b2*N2Obar)  
                    +(self.c2*CH4bar)  
                    + 0.117 ) * ( np.sqrt(n2o_ppb)-np.sqrt(n2o_ppb_0) )
        return delQ_n2o

    def collins_CH4StratoOzone(ch4_ppb,ch4_ppb_0=None):
        if ch4_ppb_0==None: ch4_ppb_0=self.ch4_ppb_0 
        delQ_Ozone_ch4=self.Collins_O3_c * (ch4_ppb-ch4_ppb_0)
        return delQ_Ozone_ch4

    def SSP_CH4StratoOzone(ch4_ppb,ch4_ppb_0=None):
        if ch4_ppb_0==None: ch4_ppb_0=self.ch4_ppb_0 
        delQ_Ozone_ch4=0.043*5*(np.log(ch4_ppb)-np.log(ch4_ppb_0))
        return delQ_Ozone_ch4








