#!/usr/bin/env python

import numpy as np
import sys
import string
import argparse
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt
import mpl_toolkits.basemap.pyproj as pyproj 



################################################################################
################################################################################
#
# Define the class
#
################################################################################
################################################################################

class llbng:

	def __init__(self):

		# projection from:
		# http://spatialreference.org/ref/epsg/27700/proj4/
		self.osgb36=pyproj.Proj('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' ) # UK Ordnance Survey, 1936 datum 

		# http://spatialreference.org/ref/epsg/4258/
		self.etrs89=pyproj.Proj('+proj=longlat +ellps=GRS80 +no_defs')

	def ll_to_bng(self,lon,lat):
		x, y = pyproj.transform(self.etrs89, self.osgb36, lon, lat )

		return x, y



	def read_lat_lon(self,fname):

		f=open(fname,'r')

		sites=[]
		for lin in f.readlines():
			if lin.strip()[0] != '#' :
				els=lin.strip().split(',')
				if len(els)>0:
					sites.append({})
					sites[-1]['sitename']=els[0]
					sites[-1]['lat']=float(els[1])
					sites[-1]['lon']=float(els[2])

		f.close()
		return sites

################################################################################
# Parse the input
################################################################################
	def parse_input(self):
		parser=argparse.ArgumentParser(description='Get BNG square(s) containing pair(s) of lat/lon coords')

		# positional
		parser.add_argument('infname',help='File containing sitename, lat, lon elev')
		parser.add_argument('outfname',help='Output file')

		args=parser.parse_args()

		return args.infname, \
				args.outfname

################################################################################
################################################################################
#
# Start the main routine
#
################################################################################
################################################################################

if __name__=='__main__':

	l=llbng()

	infname, outfname = l.parse_input()

	sites=l.read_lat_lon(infname)

	for site in sites:
		site['x'],site['y']=l.ll_to_bng(site['lon'],site['lat'])
		

	of=open(outfname,'w')
	for site in sites:
		of.write('%s   %f    %f    %f    xx\n'%(site['sitename'].replace(' ','_').replace("'",''),site['x'],site['y'],0.0))

	of.close()


