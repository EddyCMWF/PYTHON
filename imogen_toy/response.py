def response(ncallyr, nyr):
   import numpy as np

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
