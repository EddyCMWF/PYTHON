#!/usr/bin/python
#
# Python module to plot the Wetladn Dignostics from the EMEP4UK JULES output
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# January 2015
#
# Contains
#
import os, sys
import numpy as np
import argparse
import netCDF4 as nc
import matplotlib.pyplot as plot
import plot_map_ECP as PM
#
