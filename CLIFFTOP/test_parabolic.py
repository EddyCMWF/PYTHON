#!/bin/env python


import numpy as np
from parabolic import parabolic


kappa=350
DTEMP_O = np.array( [0.0,0.028262032 ] ) #, 0.06025] )

para=parabolic(kappa,DTEMP_O)


