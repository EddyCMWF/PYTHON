#### Module replicating relationships in Piere Friedlingstein's C4MIP Feedback analysis
#### Friedlingstein, P., et al, 2006: Climate-Carbon Cycle Feedback Analysis: Results 
#### from the C4MIP Model Intercomparison. J. Climate, 19, 3337-3353, https://doi.org/10.1175/JCLI3800.1

import numpy as np

class CarbonClimateFeedbacks:
    def __init__(self, delCO2ppmv, delTglob, cmip='CMIP5'): # , cmip='C4MIP'):  #
        self.delCO2ppmv   = delCO2ppmv
        self.delTglob     = delTglob    
        self.cmip         = cmip.upper()
        if self.cmip in ['C4MIP','FRIEDLINGSTEIN','CMIP3']:
            print( 'Using parameters from C4MIP models' )
            self.CMIP='C4MIP'   # class recognised CMIP name
            self.params=C4MIP_params()
        elif self.cmip in ['CMIP5','ARORA']:
            print( 'Using parameters from CMIP5 models' )
            self.CMIP='CMIP5'    # class recognised CMIP name
            self.params=CMIP5_params()
        else:
            print( 'Unrecognised CMIP name' )

        #self.delCland_Accum = self.delCland_Accum_Ensemble()
    
    # Function to return the change in carbon pool for given atmospheric CO2 (ppmv), 
    #  global warming (K) and beta and gamma values
    def delC_Accum(self, beta, gamma):
        delC_Accum = (beta*self.delCO2ppmv) + (gamma*self.delTglob)
        return delC_Accum

    def delCland_Accum(self, model):
        beta_l = self.params.model_dict[model]['beta_l']
        gamma_l= self.params.model_dict[model]['gamma_l']
        return self.delC_Accum(beta_l, gamma_l)

    def delCland_Accum_Ensemble(self):#, delCO2ppmv, delTglob):
        delCland_Accum_Ensemble = { self.params.models[imod]:
                                    self.delCland_Accum(self.params.models[imod])
                                        for imod in range(self.params.nmodels) }
        return delCland_Accum_Ensemble 

    def delCocean_Accum(self, model):
        beta_o = self.params.model_dict[model]['beta_o']
        gamma_o= self.params.model_dict[model]['gamma_o']
        return self.delC_Accum(beta_o, gamma_o)

    def delCocean_Accum_Ensemble(self):#, delCO2ppmv, delTglob):
        delCland_Accum_Ensemble = { self.params.models[imod]:
                                    self.delCocean_Accum(self.params.models[imod])
                                        for imod in range(self.params.nmodels) }
        return delCland_Accum_Ensemble 


# Friedlingstein et al. (2006)
class C4MIP_params:
    def __init__(self):
        # Set parameters and model names:
        self.models = ['HadCM3LC','IPSL-CM2C','IPSL-CM4-LOOP','CSM-1','MPI','LLNL',
                'FRCGC','UMD','UVic-2.7','CLIMBER','BERN-CC']
        self.nmodels=len(self.models)
        self.alpha=[0.0066,0.0065,0.0072,0.0038,0.0082,0.0068,0.0059,0.0063,0.0053,0.0046]
        self.beta_l=[1.3,1.6,1.3,1.1,1.4,2.8,1.2,0.2,1.2,1.1,1.6]
        self.gamma_l=[-177.,-98.,-20.,-23.,-65.,-70.,-112.,-40.,-98.,-57.,-105.]
        self.beta_o=[0.8,1.6,1.1,0.9,1.1,0.9,1.2,1.5,1.1,0.9,1.3]
        self.gamma_o=[-24.,-30.,-16.,-17.,-22.,-14.,-46.,-67.,-43.,-22.,-39.]
        # store above params into a dictionary by model 
        self.model_dict = { self.models[imodel]: { 'gamma_l':self.gamma_l[imodel],
                                         'beta_l':self.beta_l[imodel],
                                         'gamma_o':self.gamma_o[imodel],
                                         'beta_o':self.beta_o[imodel], }
                       for imodel in range(self.nmodels) }

# Arora (2013) parameters
class CMIP5_params:
    def __init__(self):
        # Set parameters and model names:
        self.models = ['MPI-ESM-LR','IPSL-CM5A-LR','BCC-CSM1','HadGEM2','UVic ESCM',
                        'CanESM2','NorESM-ME','CESM1-BGC','MIROC ESM']
        self.nmodels=len(self.models)
        self.gamma_l=[-83.2,-58.6,-77.8,-30.1,-78.5,-71.9,-15.6,-21.3,-88.6]
        self.beta_l=[1.46,1.14,1.36,1.16,0.96,0.97,0.22,0.24,0.74]
        self.gamma_o=[-9.0,-6.2,-9.8,-10.0,-7.3,-7.8,-5.7,-2.4,-12.1]
        self.beta_o=[0.83,0.91,0.83,0.79,0.78,0.69,0.85,0.72,0.82]
        self.model_dict = { self.models[imodel]: { 'gamma_l':self.gamma_l[imodel],
                                                 'beta_l':self.beta_l[imodel],
                                                 'gamma_o':self.gamma_o[imodel],
                                                 'beta_o':self.beta_o[imodel], }
                               for imodel in range(self.nmodels) }






