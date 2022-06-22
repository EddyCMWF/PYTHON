#!/bin/env python2.7

import ipdb
import numpy as np
import sys, glob, os
import zipfile as zf
from maths_tools import MetTools as MTs

site=sys.argv[1]

site_in = site.replace('_','-')
site_out= site.replace('-','_')

indir='/data/grp/fluxdata/FLUXNET_2015/raw_data/'
outdir='/data/grp/fluxdata/FLUXNET_2015/wetland_suite/'

infiletag=indir+'FLX_'+site_in+'_FLUXNET2015_FULLSET_'
inzipfile = glob.glob(infiletag+'*.zip')[0]

print(infiletag,inzipfile)

start_year=int(inzipfile[len(infiletag):len(infiletag)+4])
end_year=int(inzipfile[len(infiletag)+5:len(infiletag)+9])

zipdir = outdir+site_out+'_tempzip/'
os.system('mkdir -p '+zipdir)
zf_ref = zf.ZipFile(inzipfile,'r')
zf_ref.extractall(zipdir)
zf_ref.close()

in_csvfile = glob.glob(zipdir+'FLX_'+site_in+'_FLUXNET2015_FULLSET_HH_'
                        +str(start_year)+'-'+str(end_year)+'*.csv')[0]
out_datfile = outdir+'fluxnet/'+site_out+'-met.dat'

incsvf = open(in_csvfile,'r')
inlines = incsvf.readlines()
header = inlines.pop(0).replace('\n','')
headers = header.split(',')

tstep_len = 1800.
jules_varnames = [ 't',  'sw_down', 'lw_down', 'q', 'pstar', 'precip', 'wind', 'co2_ppmv' ]
var_formats    = [ '%6.2f,', '%7.3f,', '%7.3f,', '%8.5f,', '%8.1f,', '%10.5f,', '%6.3f,', '%7.3f\n'  ]
site_varnames  = [ 'TA_F', 'SW_IN_F', 'LW_IN_F', 'VPD_F', 'PA_F', 'P_F', 'WS_F', 'CO2_F_MDS'   ]
era_varnames   = [ 'TA_ERA', 'SW_IN_ERA', 'LW_IN_ERA', 'VPD_ERA', 'PA_ERA', 'P_ERA', 'WS_ERA', 'CO2_F_MDS'   ]
convert_functions = [ lambda x:x+273.15, lambda x:x, lambda x:x, lambda x:MTs.vpd2sh(x[0]*100.,T=x[1],P=x[2]), 
        lambda x:x*1e3, lambda x:x/tstep_len, lambda x:x, lambda x:x ]
nvars=len(jules_varnames)

var_dict = { jules_varnames[ivar]:{   #'varname':varnames[ivar],
                                   'convert':convert_functions[ivar],
                                   'format':var_formats[ivar],
                                   'ERA_index':headers.index(era_varnames[ivar]),
                                   'F_index':headers.index(site_varnames[ivar]),
                                   }
                for ivar in range(nvars) }

outf = open(out_datfile,'w')

for line in inlines:
    split=line.split(',')
    outstring = ''
    for var in jules_varnames:
        outvar = float(split[var_dict[var]['F_index']])
        if var != 'q':
            outvar = var_dict[var]['convert'](outvar)
        else:
            T = float(split[var_dict['t']['F_index']])+273.15
            P = float(split[var_dict['pstar']['F_index']])*1e3
            X = [outvar,T,P]
            outvar = var_dict[var]['convert'](X)
        outstring = outstring+var_dict[var]['format']%(outvar)
    
    outf.write(outstring)

outf.close()



