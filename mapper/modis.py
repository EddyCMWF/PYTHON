#!/bin/env python

"""Module of MODIS info a bit like modis.f90 really."""

class rec:
    def __init__(self,recData,recNums,recClim,name):
        self.recData = recData
        self.recNums = recNums
        self.recClim = recClim
        self.name = name

allt = rec(0,1,10,"allt")
forest = rec(2,3,11,"tree")
urban = rec(4,5,12,"urbn")
crop = rec(6,7,13,"crop")
misc = rec(8,9,14,"misc")

landTypes = [allt,forest,urban,crop,misc]

NPIX_MIN = 100.0
FMDI = -1.0E20
DTYPE = ">f4"
NFIELDS = 15

##
## End of script.
##
