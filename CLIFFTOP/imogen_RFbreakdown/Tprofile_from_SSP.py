#!/bin/env python2.7

#import netCDF4 as nc
import numpy as np
#import matplotlib.pyplot as plt
#import sys,os
import ipdb

from imogen import delQ
from maths_tools import SmoothTools as STs

start_year,end_year=1850,2100
n_yrs=end_year-start_year

SCENARIO='1p5deg'
#SCENARIO_DIR='/users/eow/edwcom/CLIFFTOP/IMOGEN/scenarios/'
SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/'

HadCRUT_fname=SCENARIO_DIR+'HadCRUT.4.5.0.0.annual_ns_avg_smooth.rtf'

RCP_in_CO2file=SCENARIO_DIR+'rcp2.6_concs_co2_vn1p1.txt'
RCP_in_CH4N2Ofile=SCENARIO_DIR+'rcp2.6_concs_ch4_n2o_vn1p1.txt'
RCP_in_Aerosolfile=SCENARIO_DIR+'rcp2.6_dq_aerosol.txt'
RCP_in_NonCO2file=SCENARIO_DIR+'rcp2.6_qnonco2_vn1p1.txt'
RCP_in_Qco2ch4n2ofile=SCENARIO_DIR+'rcp2.6_q_co2_ch4_n2o_vn1p0.txt'
RCP_in_QstratoOzonefile=SCENARIO_DIR+'rcp2.6_qOzone_qstratoWV.txt'

SSP_in_Tfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_temp_anomaly_raw.dat'
SSP_in_CO2file=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2_raw.txt'
SSP_in_CH4N2Ofile=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o_raw.txt'
SSP_in_Aerosolfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_qaerosol_raw.txt'
SSP_in_NonCO2file=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2_raw.txt'
SSP_in_Qco2ch4n2ofile=SCENARIO_DIR+'SSP2-2.6_IMAGE_q_co2_ch4_n2o_vn1p0.txt'

SSP_out_Tfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_temp_anomaly' #.dat'
SSP_out_CO2file=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2'# .txt'
SSP_out_CH4N2Ofile=SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o'#.txt'
SSP_out_Aerosolfile=SCENARIO_DIR+'SSP2-2.6_IMAGE_qaerosol'#.txt'
SSP_out_NonCO2file=SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2'#.txt'

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
#ipdb.set_trace()
for i_yr in range(rcp_nyr):
   print(i_yr)
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
SSP_Tlines=open(SSP_in_Tfile,'r').readlines()
SSP_co2lines=open(SSP_in_CO2file,'r').readlines()
SSP_ch4n2olines=open(SSP_in_CH4N2Ofile,'r').readlines()
SSP_aerosollines=open(SSP_in_Aerosolfile,'r').readlines()
SSP_nonco2lines=open(SSP_in_NonCO2file,'r').readlines()
SSP_qco2ch4n2olines=open(SSP_in_Qco2ch4n2ofile,'r').readlines()
ssp_nyr=len(SSP_co2lines)
SSP_years_in=np.zeros(ssp_nyr)
SSP_T_in=np.zeros(ssp_nyr)
SSP_co2_in=np.zeros(ssp_nyr)
SSP_ch4_in=np.zeros(ssp_nyr)
SSP_n2o_in=np.zeros(ssp_nyr)
SSP_aerosol_in=np.zeros(ssp_nyr)
SSP_nonco2_in=np.zeros(ssp_nyr)
SSP_qco2_in=np.zeros(ssp_nyr)
SSP_qch4_in=np.zeros(ssp_nyr)
SSP_qn2o_in=np.zeros(ssp_nyr)
for i_yr in range(ssp_nyr):
   SSP_years_in[i_yr]=float(SSP_co2lines[i_yr].split()[0]) 
   SSP_T_in[i_yr]=float(SSP_Tlines[i_yr].split()[1]) 
   SSP_co2_in[i_yr]=float(SSP_co2lines[i_yr].split()[1]) 
   SSP_ch4_in[i_yr]=float(SSP_ch4n2olines[i_yr].split()[1]) 
   SSP_n2o_in[i_yr]=float(SSP_ch4n2olines[i_yr].split()[2]) 
   SSP_aerosol_in[i_yr]=float(SSP_aerosollines[i_yr].split()[1]) 
   SSP_nonco2_in[i_yr]=float(SSP_nonco2lines[i_yr].split()[1]) 
   SSP_qco2_in[i_yr]=float(SSP_qco2ch4n2olines[i_yr].split()[1])
   SSP_qch4_in[i_yr]=float(SSP_qco2ch4n2olines[i_yr].split()[2])
   SSP_qn2o_in[i_yr]=float(SSP_qco2ch4n2olines[i_yr].split()[3])

# Substitute the etminan ch4 and n2o and Collin WV and Ozone
SSP_qn2o_etminan=delQ.etminan_N2O(SSP_n2o_in,SSP_co2_in,SSP_ch4_in)
SSP_qch4_etminan=delQ.etminan_CH4(SSP_ch4_in,SSP_n2o_in) #,ch4_ppb_0=722.,n2o_ppb_0=270.)

SSP_T_in[-2:]=SSP_T_in[-3]

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
delT_1p5 = 1.5-HadCRUT_T_trans
delT_2   = 2.-HadCRUT_T_trans

SSP_T_a=SSP_T_smooth-SSP_T_smooth[SSP_start_index]
SSP_co2_a=SSP_co2_smooth-SSP_co2_smooth[SSP_start_index]
SSP_aerosol_a=SSP_aerosol_smooth-SSP_aerosol_smooth[SSP_start_index]
SSP_nonco2_a=SSP_nonco2_etminan_smooth-SSP_nonco2_etminan_smooth[SSP_start_index]

ratio_1p5 = delT_1p5/SSP_T_a[-1]
ratio_2   = delT_2/SSP_T_a[-1]

SSP_T_raw=np.zeros_like(RCP_co2)
SSP_T_raw[:RCP_start_index]=HadCRUT_T[:RCP_start_index]
SSP_T_raw[RCP_start_index:]=HadCRUT_T_trans+(SSP_T_a[SSP_start_index:])
SSP_T_raw[RCP_start_index-10:]=STs.meanbox_iterate(SSP_T_raw[RCP_start_index-10:],iterations=10)
outf=open(SSP_out_Tfile+'.dat','w')
for i_yr in range(len(SSP_T_raw)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_T_raw[i_yr]))
outf.close()
quit()

SSP_T_1p5=np.zeros_like(RCP_co2)
SSP_T_1p5[:RCP_start_index]=HadCRUT_T[:RCP_start_index]
SSP_T_1p5[RCP_start_index:]=HadCRUT_T_trans+(SSP_T_a[SSP_start_index:]*ratio_1p5)
SSP_T_1p5[RCP_start_index-10:]=STs.meanbox_iterate(SSP_T_1p5[RCP_start_index-10:],iterations=10)
outf=open(SSP_out_Tfile+'_1p5deg.dat','w')
for i_yr in range(len(SSP_T_1p5)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_T_1p5[i_yr]))
outf.close()


SSP_aero_1p5=np.zeros_like(RCP_co2)
SSP_aero_1p5[:RCP_start_index]=RCP_aerosol[:RCP_start_index]
SSP_aero_1p5[RCP_start_index:]=RCP_aerosol[RCP_start_index]+\
                                (SSP_aerosol_a[SSP_start_index:]*ratio_1p5)
SSP_aero_1p5[RCP_start_index-5:]=STs.meanbox_iterate(SSP_aero_1p5[RCP_start_index-5:],
                                                       iterations=10)
outf=open(SSP_out_Aerosolfile+'_1p5deg.dat','w')
for i_yr in range(len(SSP_aero_1p5)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_aero_1p5[i_yr]))
outf.close()


SSP_nonco2_1p5=np.zeros_like(RCP_co2)
SSP_nonco2_1p5[:RCP_start_index]=RCP_nonco2_etminan[:RCP_start_index]
SSP_nonco2_1p5[RCP_start_index:]=RCP_nonco2_etminan[RCP_start_index]+\
                                (SSP_nonco2_a[SSP_start_index:]*ratio_1p5)
SSP_nonco2_1p5[RCP_start_index-5:]=STs.meanbox_iterate(SSP_nonco2_1p5[RCP_start_index-5:],
                                                       iterations=10)
outf=open(SSP_out_NonCO2file+'_1p5deg.dat','w')
for i_yr in range(len(SSP_nonco2_1p5)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_nonco2_1p5[i_yr]))
outf.close()

quit()

SSP_T_2=np.zeros_like(RCP_co2)
SSP_T_2[:RCP_start_index]=HadCRUT_T[:RCP_start_index]
SSP_T_2[RCP_start_index:]=HadCRUT_T_trans+(SSP_T_a[SSP_start_index:]*ratio_2)
SSP_T_2[RCP_start_index-10:]=STs.meanbox_iterate(SSP_T_2[RCP_start_index-10:],iterations=10)

outf=open(SSP_out_Tfile+'_2deg.dat','w')
for i_yr in range(len(SSP_T_2)):
    outf.write('%4i %8.4f\n'%(RCP_years[i_yr],SSP_T_2[i_yr]))
outf.close()

quit()


SSP_co2_1p5=np.zeros_like(RCP_co2)
SSP_co2_1p5[:SSP_start_index]=RCP_co2[:SSP_start_index]
SSP_co2_1p5[SSP_start_index:]=RCP_co2[SSP_start_index]+(SSP_co2_a*ratio_1p5)





SSP_ch4_1p5=np.zeros_like(RCP_co2)
SSP_n2o_1p5=np.zeros_like(RCP_co2)
SSP_nonco2_1p5=np.zeros_like(RCP_co2)
SSP_aerosol_1p5=np.zeros_like(RCP_co2)





SSP_T_2=np.zeros_like(RCP_co2)
SSP_co2_2=np.zeros_like(RCP_co2)
SSP_ch4_2=np.zeros_like(RCP_co2)
SSP_n2o_2=np.zeros_like(RCP_co2)
SSP_nonco2_2=np.zeros_like(RCP_co2)
SSP_aerosol_2=np.zeros_like(RCP_co2)






