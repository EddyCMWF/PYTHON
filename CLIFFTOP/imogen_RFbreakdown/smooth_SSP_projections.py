#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
import numpy as np
from copy import deepcopy
from imogen import delQ

# Smooth options:
smooth_window = 2
smooth_Npass = 5

version_tag = 'vn3p0'
SCENARIO_DIR='/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/ssp_projections/'
out_dir = SCENARIO_DIR+'SYNTHESIS_projections/'
raw_SSP_dir = SCENARIO_DIR+'raw_SSP_data/'

SCENARIOS = [ 'SSP2-1.9_IMAGE', 'SSP2-baseline_IMAGE' ]

projections = { 'concs_ch4_n2o': { 'template_file':SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_ch4_n2o.txt',
                                   'store_vars':['ch4','n2o'] },
                'qnonco2': { 'template_file':SCENARIO_DIR+'SSP2-2.6_IMAGE_qnonco2.txt',
                              'store_vars':['qnonco2'] },
                'concs_co2': { 'template_file':SCENARIO_DIR+'SSP2-2.6_IMAGE_concs_co2.txt',
                              'store_vars':['co2'] },
                }


DFs = {}
for scenario in SCENARIOS:
    in_rawfile = raw_SSP_dir+scenario+'.csv'
    inf_raw = open(in_rawfile,'r')
    inlines_raw = inf_raw.readlines()
    inf_raw.close()
    print('Reading raw data for: ', inlines_raw.pop(0))
    
    years = inlines_raw.pop(0).split(',')[1:]
    years = [int(year) for year in years]
    years_obj = [ dt.date(year,1,1) for year in years]
    in_raw_data={}
    for line in inlines_raw:
        split = line.split(',')
        in_raw_data[split[0]] = [float(val) for val in split[1:]] 
    in_raw_df = pd.DataFrame(in_raw_data,index=years_obj)
    SP_obj = in_raw_df.index[0]
    EP_obj = in_raw_df.index[-1]
    scen_DF_list = []    
    for project in projections.keys():
        proj_dict = projections[project]
        template_inf = open(proj_dict['template_file'], 'r')
        template_lines = template_inf.readlines()
        template_inf.close()
        temp_data_dict = { var:[] for var in proj_dict['store_vars'] }
        temp_years_obj = []
        for line in template_lines:
            split = line.split()
            temp_years_obj.append(dt.date(int(split[0]),1,1))
            for val,var in zip(split[1:],proj_dict['store_vars']):
                temp_data_dict[var].append(float(val))

        output_df0=  pd.DataFrame(temp_data_dict,index=temp_years_obj)
        output_df =  pd.DataFrame(temp_data_dict,index=temp_years_obj)
        output_df2=  pd.DataFrame(temp_data_dict,index=temp_years_obj)
        reindex = output_df[var][SP_obj:EP_obj].index
        SP = np.where(output_df.index == SP_obj)[0][0]
        EP = np.where(output_df.index == EP_obj)[0][0]
        
        for var in proj_dict['store_vars']:
            start_diff = output_df[var][SP_obj] - in_raw_df[var][SP_obj]
            in_raw_df[var] += start_diff
            output_df[var][SP_obj:EP_obj] = in_raw_df[var].reindex(reindex).interpolate()
            output_df[var][EP_obj:] =output_df[var][EP_obj]
            output_df2[var][SP_obj:EP_obj] = in_raw_df[var].reindex(reindex).interpolate()
            output_df2[var][EP_obj:] =output_df2[var][EP_obj]
          
            for counter in range(smooth_Npass):
                for ipt in range(SP, EP+1):
                    m_index=output_df.index[ipt]
                    s_index=output_df.index[ipt-smooth_window]
                    e_index=output_df.index[ipt+smooth_window]
                    output_df[var][m_index] = np.mean(output_df2[var][s_index:e_index])
                output_df2[var][SP:EP+1] = output_df[var][SP:EP+1]
            # Now rescale to maintain any peak in time-series:
            #output_df[var][SP:EP+1] = output_df[var][SP:EP+1] * (in_raw_df[var].max()/output_df[var][SP:EP+1].max())
        scen_DF_list.append(output_df)

        # output smooth data to file:
        noutvars = len(proj_dict['store_vars'])
        outf = open(out_dir+scenario+'_'+project+'_'+version_tag+'.txt', 'w')
        for ipt in output_df2.index:
            outf.write("%4i"%(ipt.year) + 
                 noutvars*" %8.3f"%tuple([output_df[var][ipt] for var in proj_dict['store_vars']])+'\n' )
        outf.close()
        fig,axes = plt.subplots(nrows=noutvars,figsize=(8,4*noutvars))
        if not hasattr(axes,'__iter__'): axes = [axes]
        for i,var in enumerate(proj_dict['store_vars']):
            output_df0[var].plot(ax=axes[i],label='Template')
            in_raw_df[var].plot(ax=axes[i],label='SSP-raw')
            output_df[var].plot(ax=axes[i],label='Smoothed')
        handles,labels = axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc=8,ncol=3)
        fig.savefig(out_dir+scenario+'_'+project+'_'+version_tag+'.png',bbox_inches='tight')
        plt.close()
    
    DFs[scenario] = pd.concat(scen_DF_list,axis=1)

BigDF = pd.concat(DFs,axis=1)

# Now create qnonCO2 for RCP1.9 with baseline CH4:

rcp19_N2O = BigDF[SCENARIOS[0]]['n2o']
rcp19_CH4 = BigDF[SCENARIOS[0]]['ch4']

rcpBL_CH4 = BigDF[SCENARIOS[1]]['ch4']

# Write out rcpBL_CH4 and rcp19_N2O to file
outf = open(out_dir+SCENARIOS[0]+'-baselineCH4_conc_ch4_n2o_'+version_tag+'.txt','w')
#outf = open(out_dir+'SSP2-rcpBL_CH4-rcp19_N2O_'+version_tag+'.txt','w')
for iyear in BigDF.index:
    outf.write('%4i %8.3f %8.3f\n'%(iyear.year,rcpBL_CH4[iyear],rcp19_N2O[iyear]))
outf.close()

#CH4_diff = BigDF[SCENARIOS[1]]['ch4'] - BigDF[SCENARIOS[0]]['ch4']

qCH4_diff = delQ.etminan_CH4(rcpBL_CH4, rcp19_N2O, ch4_ppb_0=rcp19_CH4, n2o_ppb_0=rcp19_N2O)
#qCH4_diff = delQ.etminan_CH4(BigDF[SCENARIOS[1]]['ch4'], BigDF[SCENARIOS[0]]['n2o'],
#                        ch4_ppb_0=BigDF[SCENARIOS[0]]['ch4'], n2o_ppb_0=BigDF[SCENARIOS[0]]['n2o'])
qCH4ozone_diff = delQ.collins_CH4StratoOzone(rcpBL_CH4, ch4_ppb_0=rcp19_CH4)
#qCH4ozone_diff = delQ.collins_CH4StratoOzone(BigDF[SCENARIOS[1]]['ch4'], ch4_ppb_0=BigDF[SCENARIOS[0]]['ch4'])
qnonCO2_merged = BigDF[SCENARIOS[0]]['qnonco2'] + qCH4_diff + qCH4ozone_diff
qnonCO2_merged2 = deepcopy(qnonCO2_merged) 
qnonCO2_merged0 = deepcopy(qnonCO2_merged) 

# smooth the new qnonco2 profile:
for counter in range(smooth_Npass**2):
    for ipt in range(SP, EP+1):
        m_index=output_df.index[ipt]
        s_index=output_df.index[ipt-smooth_window]
        e_index=output_df.index[ipt+smooth_window]
        qnonCO2_merged[m_index] = np.mean(qnonCO2_merged2[s_index:e_index])
    qnonCO2_merged2 = deepcopy(qnonCO2_merged)

#pd.concat([qnonCO2_merged0,qnonCO2_merged,qnonCO2_merged2],axis=1).plot()
fig,ax = plt.subplots(figsize=(8,4))
qnonCO2_merged0.plot(ax=ax,label='original')
qnonCO2_merged.plot(ax=ax,label='smoothed')
fig.legend(loc=8,ncol=2)
fig.savefig(out_dir+SCENARIOS[0]+'-baselineCH4_qnonco2_'+version_tag+'.png',bbox_inches='tight')
plt.close()

# Write out SSP2-1.9 qnonco2 with additional SSP2-BL CH4 contributions
outf = open(out_dir+SCENARIOS[0]+'-baselineCH4_qnonco2_'+version_tag+'.txt','w')
for iyear in BigDF.index:
    outf.write('%4i %8.3f\n'%(iyear.year, qnonCO2_merged[iyear]))
outf.close()



