#!/bin/env python2.7

from imogen import CMIP_CCFBs as CCFB
import numpy as np

delCO2=np.arange(0,100,1)
delT=np.arange(0,1,0.01)

test=CCFB.CarbonClimateFeedbacks(delCO2,delT)

#print(test.delCland_Accum_Ensemble())
print(test.delCocean_Accum_Ensemble())

#print(test.delCocean_Accum('HadGEM2'))


