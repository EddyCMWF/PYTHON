#!/bin/env python3

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys,os

from imogen import delQ
from maths_tools import SmoothTools as STs

start_year,end_year=1850,2100
n_yrs=end_year-start_year

#SCENARIO='1p5deg'
SCENARIO='1p81p5deg'
#SCENARIO='2deg'
SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'

if SCENARIO=='1p5deg':
    SSP_scenario='IMAGE-SSP2-26'
    maxT=1.5
elif SCENARIO=='1p81p5deg':
    SSP_scenario='MESSAGE-GLOBIOM-SSP2-26'
    maxT=1.8
elif SCENARIO=='2deg':
    SSP_scenario='IMAGE-SSP2-34'
    maxT=2.0
else: 
    print('Unrecognised Scenario')

HadCRUT_fname=SCENARIO_DIR+'HadCRUT.4.5.0.0.annual_ns_avg_smooth.rtf'

RCP_in_CO2file=SCENARIO_DIR+'rcp2.6_concs_co2_vn1p0.txt'
RCP_in_CH4N2Ofile=SCENARIO_DIR+'rcp2.6_concs_ch4_n2o_vn1p0.txt'
RCP_in_Aerosolfile=SCENARIO_DIR+'rcp2.6_dq_aerosol.txt'
RCP_in_NonCO2file=SCENARIO_DIR+'rcp_qnonco2_vn1p0_2.6.txt'
RCP_in_Qco2ch4n2ofile=SCENARIO_DIR+'rcp2.6_q_co2_ch4_n2o_vn1p0.txt'
RCP_in_QstratoOzonefile=SCENARIO_DIR+'rcp2.6_qOzone_qstratoWV.txt'

SSP_in_file=SCENARIO_DIR+'SSP_data/'+SSP_scenario+'.csv'

SSP_out_Tfile=SCENARIO_DIR+SSP_scenario+'_temp_anomaly' #.dat'
SSP_out_CO2file=SCENARIO_DIR+SSP_scenario+'_concs_co2'# .txt'
SSP_out_CH4N2Ofile=SCENARIO_DIR+SSP_scenario+'_concs_ch4_n2o'#.txt'
SSP_out_Aerosolfile=SCENARIO_DIR+SSP_scenario+'_qaerosol'#.txt'
SSP_out_NonCO2file=SCENARIO_DIR+SSP_scenario+'_qnonco2'#.txt'

#Read in Hadcrut data
HadCRUT_lines=open(HadCRUT_fname,'r').readlines()
HadCRUT_T=np.array([float(line.split()[1]) for line in HadCRUT_lines])
HadCRUT_years=np.array([float(line.split()[0]) for line in HadCRUT_lines])
# Normalise T to median of first 51 years
median_51 = np.median(HadCRUT_T[:51])
HadCRUT_T-=median_51




#Read RCP data in:
rcp_start,rcp_end=1850,2100
rcp_nyr=rcp_end-rcp_start+1
RCP_co2lines=open(RCP_in_CO2file,'r').readlines()
RCP_ch4n2olines=open(RCP_in_CH4N2Ofile,'r').readlines()
RCP_aerosollines=open(RCP_in_Aerosolfile,'r').readlines()
RCP_nonco2lines=open(RCP_in_NonCO2file,'r').readlines()
RCP_qco2ch4n2olines=open(RCP_in_Qco2ch4n2ofile,'r').readlines()
#RCP_qOzoneWVlines=open(RCP_in_QstratoOzonefile,'r').readlines()
RCP_co2=np.zeros(rcp_nyr)
RCP_ch4=np.zeros(rcp_nyr)
RCP_n2o=np.zeros(rcp_nyr)
RCP_aerosol=np.zeros(rcp_nyr)
RCP_nonco2=np.zeros(rcp_nyr)
RCP_qco2=np.zeros(rcp_nyr)
RCP_qch4=np.zeros(rcp_nyr)
RCP_qn2o=np.zeros(rcp_nyr)
#RCP_qozone=np.zeros(rcp_nyr)
#RCP_qstratoWV=np.zeros(rcp_nyr)
RCP_years=np.arange(rcp_start,rcp_end+1)
for i_yr in range(rcp_nyr):
   RCP_co2[i_yr]=float(RCP_co2lines[i_yr].split()[1]) 
   RCP_ch4[i_yr]=float(RCP_ch4n2olines[i_yr].split()[1]) 
   RCP_n2o[i_yr]=float(RCP_ch4n2olines[i_yr].split()[2]) 
   RCP_aerosol[i_yr]=float(RCP_aerosollines[i_yr].split()[1]) 
   RCP_nonco2[i_yr]=float(RCP_nonco2lines[i_yr].split()[1]) 
   RCP_qco2[i_yr]=float(RCP_qco2ch4n2olines[i_yr].split()[1])
   RCP_qch4[i_yr]=float(RCP_qco2ch4n2olines[i_yr].split()[2])
   RCP_qn2o[i_yr]=float(RCP_qco2ch4n2olines[i_yr].split()[3])
   #RCP_qozone[i_yr]=float(RCP_qOzoneWVlines[i_yr].split()[1])
   #RCP_qstratoWV[i_yr]=float(RCP_qOzoneWVlines[i_yr].split()[2])


# Substitute the etminan ch4 and n2o
RCP_qn2o_etminan=delQ.etminan_N2O(RCP_n2o,RCP_co2,RCP_ch4) #,
#                        ch4_ppb_0=RCP_ch4[0],n2o_ppb_0=RCP_n2o[0],co2_ppm_0=RCP_co2[0])
RCP_qch4_etminan=delQ.etminan_CH4(RCP_ch4,RCP_n2o)#,ch4_ppb_0=RCP_ch4[0],n2o_ppb_0=RCP_n2o[0])
#RCP_qSozone_collins=delQ.collins_CH4StratoOzone(RCP_ch4)

RCP_nonco2_etminan = RCP_nonco2-RCP_qch4        -RCP_qn2o        \
                                +RCP_qch4_etminan+RCP_qn2o_etminan
#RCP_qnonco2_etminan = RCP_nonco2-RCP_qch4        -RCP_qn2o        -RCP_qozone-RCP_qstratoWV   \
#                                +RCP_qch4_etminan+RCP_qn2o_etminan+RCP_qSozone_collins


#Read SSP data in:
SSP_lines=open(SSP_in_file,'r').readlines()
for line in SSP_lines[1:]:
    split=line.split(',')
    if split[0]=='year':
        SSP_years_in=np.array([ int(yr) for yr in split[1:]])
    elif split[0]=='T':
        SSP_T_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='co2':
        SSP_co2_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='ch4':
        SSP_ch4_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='n2o':
        SSP_n2o_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='qnonco2':
        SSP_nonco2_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='qCO2':
        SSP_qco2_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='qCH4':
        SSP_qch4_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='qN2O':
        SSP_qn2o_in=np.array([ float(dat) for dat in split[1:]])
    elif split[0]=='qAerosol':
        SSP_aerosol_in=np.array([ float(dat) for dat in split[1:]])

ssp_nyr=len(SSP_T_in)

# Create aymptote, scenario specific:
if SCENARIO=='1p5deg':
    SSP_T_in[-2:]=SSP_T_in[-3]
elif SCENARIO=='2deg':
    SSP_T_in[-1:]=SSP_T_in[-2]

# Interpolate SSP data to yearly
SSP_years_interp=np.arange(SSP_years_in[0],SSP_years_in[-1]+1)
SSP_T_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_T_in)
SSP_co2_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_co2_in)
SSP_ch4_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_ch4_in)
SSP_n2o_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_n2o_in)
SSP_aerosol_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_aerosol_in)
SSP_nonco2_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_nonco2_in)
SSP_qco2_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_qco2_in)
SSP_qch4_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_qch4_in)
SSP_qn2o_interp=np.interp(SSP_years_interp,SSP_years_in,SSP_qn2o_in)

# Smooth data 
SSP_T_smooth=STs.meanbox_iterate(SSP_T_interp,iterations=10)
SSP_co2_smooth=STs.meanbox_iterate(SSP_co2_interp,iterations=10)
SSP_ch4_smooth=STs.meanbox_iterate(SSP_ch4_interp,iterations=10)
SSP_n2o_smooth=STs.meanbox_iterate(SSP_n2o_interp,iterations=10)
SSP_aerosol_smooth=STs.meanbox_iterate(SSP_aerosol_interp,iterations=10)
SSP_nonco2_smooth=STs.meanbox_iterate(SSP_nonco2_interp,iterations=10)
SSP_qco2_smooth=STs.meanbox_iterate(SSP_qco2_interp,iterations=10)
SSP_qch4_smooth=STs.meanbox_iterate(SSP_qch4_interp,iterations=10)
SSP_qn2o_smooth=STs.meanbox_iterate(SSP_qn2o_interp,iterations=10)

# Substitute the etminan ch4 and n2o and Collin WV and Ozone
SSP_qn2o_etminan_smooth=delQ.etminan_N2O(SSP_n2o_smooth,SSP_co2_smooth,SSP_ch4_smooth)
SSP_qch4_etminan_smooth=delQ.etminan_CH4(SSP_ch4_smooth,SSP_n2o_smooth)

#SSP_ozoneWV = delQ.SSP_CH4StratoOzone(SSP_ch4_interp)
#SSP_ozoneWV_collins = delQ.collins_CH4StratoOzone(SSP_ch4_interp)

SSP_nonco2_etminan_smooth = SSP_nonco2_smooth-SSP_qch4_smooth-SSP_qn2o_smooth \
                                +SSP_qch4_etminan_smooth+SSP_qn2o_etminan_smooth

splice_year=2005
SSP_start_index = np.where(SSP_years_interp==splice_year)[0][0]
RCP_start_index = np.where(RCP_years==splice_year)[0][0]

HadCRUT_T_trans=HadCRUT_T[RCP_start_index]
delT = maxT-HadCRUT_T_trans

SSP_T_a=SSP_T_smooth-SSP_T_smooth[SSP_start_index]
ratio = delT/SSP_T_a[-1]

SSP_T=np.zeros_like(RCP_co2)
SSP_T[:RCP_start_index]=HadCRUT_T[:RCP_start_index]
SSP_T[RCP_start_index:]=HadCRUT_T_trans+(SSP_T_a[SSP_start_index:]*ratio)
SSP_T[RCP_start_index-10:]=STs.meanbox_iterate(SSP_T[RCP_start_index-10:],iterations=10)
outf=open(SSP_out_Tfile+'_'+SCENARIO+'.dat','w')
for i_yr in range(len(SSP_T)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_T[i_yr]))
outf.close()

quit()


SSP_aero=np.zeros_like(RCP_co2)
SSP_aero[:RCP_start_index]=RCP_aerosol[:RCP_start_index]
SSP_aero[RCP_start_index:]=RCP_aerosol[RCP_start_index]+\
                                (SSP_aerosol_a[SSP_start_index:]*ratio)
SSP_aero[RCP_start_index-5:]=STs.meanbox_iterate(SSP_aero[RCP_start_index-5:],
                                                       iterations=10)
outf=open(SSP_out_Aerosolfile+'_'+SCENARIO+'.dat','w')
for i_yr in range(len(SSP_aero)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_aero[i_yr]))
outf.close()

