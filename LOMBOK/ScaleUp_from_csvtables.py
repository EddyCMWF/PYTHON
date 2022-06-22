#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:09:41 2019

@author: edwcom
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.patches as patches

def box_and_whisker(datastats,xpos,ecolor,ax,boxwidth=0.1,
                    med_color=None,fill=False,fcolor=False,hatch=None):
    if med_color==None: med_color=ecolor
    if fcolor==None: fcolor=ecolor
    patch = patches.Rectangle( (xpos-(boxwidth/2.),datastats['Q1']),
                                boxwidth, datastats['Q3']-datastats['Q1'],
                                fill=fill, edgecolor=ecolor, facecolor=fcolor, 
                                linewidth=1.5, hatch=hatch)
    line = ax.add_patch(patch)
    #line = ax.plot( [xpos-(boxwidth/2.),xpos+(boxwidth/2.)],
    ax.plot( [patch.get_x(), patch.get_x()+patch.get_width()],
             [datastats['Median'],datastats['Median']], c=med_color, lw=2)
    #ax.plot( [xpos,xpos],[datastats['Q3'],datastats['Maximum']], c=color, lw=1.5)
    #ax.plot( [xpos,xpos],[datastats['Q1'],datastats['Minimum']], c=color, lw=1.5)
    return line



data_dir = '/users/eow/edwcom/LOMBOK/ScaleUpObs/'

land_cover_file = data_dir+'land_cover.csv'
emission_rates_file = data_dir+'emission_rates.csv'

# Read in emissions rates inmto panda dataframe
ER_inf = open(emission_rates_file)
ER_lines = ER_inf.readlines()
ER_inf.close()
ER_headers = ER_lines.pop(0).replace('\n','').split(',')
ER_data = {} # hdr:[] for hdr in headers } 
for line in ER_lines:
    split=line.replace('\n','').replace('\r','').split(',')
    ER_data[split[0]+'-'+split[2]] = [ float(val) for val in split[3:]]
    ER_data[split[0]+'-'+split[2]].append(split[1])
ER_DF = pd.DataFrame(ER_data,index=ER_headers[3:]+[ER_headers[1]])


# Read land covers into panda Dataframe
LC_inf = open(land_cover_file)
LC_lines = LC_inf.readlines()
LC_inf.close()
LC_headers = LC_lines.pop(0).replace('\n','').replace('\r','').split(',')
LC_data = { hdr:[] for hdr in LC_headers[1:]}
years = []
for line in LC_lines:
    split=line.replace('\n','').replace('\r','').split(',')
    for hdr,val in zip(LC_headers,split):
        if hdr!='Year':
            LC_data[hdr].append(float(val))
        else:
            years.append(int(val))

LC_DF = pd.DataFrame(LC_data,index=years)

# conversion factors
m2_to_ha = 1e-4
ug_to_kg = 1e-9
hr_to_yr = 1./(365.*24.)
kg_to_Mt = 1e-6
ug_to_Mt = ug_to_kg * kg_to_Mt
scale_factor = ug_to_Mt/(m2_to_ha*hr_to_yr)

# Construct N2O Dataframe, 
#  Index = Years
#  Columns = N2O-Oil Palm, N2O-All Forest / [Minimum,Q1,Median,Q3,Maximum,Mean]
#stats = ['Q1', 'Median', 'Q3']
stats = ['Minimum','Q1', 'Median', 'Q3','Maximum']
LC_type = 'Oil Palm'
ER_var  = 'N2O-N-OP'
#N2O_OP_DF = pd.concat( [ LC_DF[LC_type]*ER_DF[ER_var][stat]*scale_factor
                            #for stat in stats], axis=0, keys=stats )
N2O_OP_DF = pd.concat( [ LC_DF[LC_type][year]*ER_DF[ER_var][stats]*scale_factor
                            for year in years], axis=0, keys=years )
N2O_OP_DF.columns=stats
LC_type = 'All Forest (adjusted)'
ER_var  = 'N2O-N-Forest'
#N2O_AllF_DF = pd.concat( [ LC_DF[LC_type]*ER_DF[ER_var][stat]*scale_factor
#                              for stat in stats],axis=1 )
N2O_AllF_DF = pd.concat( [ LC_DF[LC_type][year]*ER_DF[ER_var][stats]*scale_factor
                            for year in years], axis=0, keys=years )
N2O_AllF_DF.columns=stats
N2O_DF = pd.concat({'Forest':N2O_AllF_DF, 'Oil Palm':N2O_OP_DF,
                    'Total': N2O_AllF_DF+N2O_OP_DF},axis=1)


# Construct CH4 Dataframe, 
#  Index = Years
#  Columns = N2O-Oil Palm, N2O-All Forest / [Minimum,Q1,Median,Q3,Maximum,Mean]
stats = ['Minimum','Q1', 'Median', 'Q3','Maximum']
LC_type = 'Oil Palm'
ER_var  = 'CH4-OP'
#CH4_OP_DF = pd.concat( [ LC_DF[LC_type]*ER_DF[ER_var][stat]*scale_factor
#                            for stat in stats],axis=1 )
CH4_OP_DF = pd.concat( [ LC_DF[LC_type][year]*ER_DF[ER_var][stats]*scale_factor
                            for year in years], axis=0, keys=years )
CH4_OP_DF.columns=stats
LC_type = 'All Forest (adjusted)'
ER_var  = 'CH4-Forest'
CH4_AllF_DF = pd.concat( [ LC_DF[LC_type]*ER_DF[ER_var][stat]*scale_factor
                              for stat in stats],axis=1 )
CH4_AllF_DF = pd.concat( [ LC_DF[LC_type][year]*ER_DF[ER_var][stats]*scale_factor
                            for year in years], axis=0, keys=years )
CH4_AllF_DF.columns=stats
CH4_DF = pd.concat({'Forest':CH4_AllF_DF, 'Oil Palm':CH4_OP_DF, 
                    'Total':CH4_AllF_DF+CH4_OP_DF},axis=1)




# Plot options:
LC_names = ['Forest', 'Oil Palm', 'Total' ]
LCarea_names = ['All Forest (adjusted)', 'Oil Palm' ]
Legend_labels = [LC+' emissions' for LC in LC_names]+['Forest Area', 'Oil Palm Area' ]
LC_colours  = ['darkgreen', 'gold', 'darkred']  
LC_fcolours = ['lightgreen', 'lightyellow', 'pink']  
LC_xpos_ref = [-0.2, -0.05, 0.15]  
LC_boxwidth = [ 0.1, 0.1, 0.2]
LCarea_xpos_ref = 0.35
LCarea_width = 0.03
nLCs = len(LC_names)
nyears = len(years)

fig,axes = plt.subplots(nrows=2,ncols=1,figsize=[12,8])
Species_Names = ['N$_2$O', 'CH$_4$']
Species_DFs = [N2O_DF,CH4_DF]
nSpecies = len(Species_Names)
for iSp in range(nSpecies):
    ax = axes[iSp]
    # second axis Land Cover area 
    ax2=ax.twinx()
    DF=Species_DFs[iSp]
    for i in range(nyears):
        year=years[i]
        xpos_centre = i+0.5
        lines=[]
        for j in range(nLCs):
            LCname=LC_names[j]
            lines.append(box_and_whisker(DF[LCname][year], xpos_centre+LC_xpos_ref[j], 
                                         LC_colours[j], ax, boxwidth=LC_boxwidth[j], 
                                         fill=True, fcolor=LC_fcolours[j]))
        basept = 0.0
        for j in range(nLCs-1):
            LCname=LC_names[j]
            LCareaname=LCarea_names[j]
            lines.append(ax2.bar(xpos_centre+LCarea_xpos_ref, LC_DF[LCareaname][year]*m2_to_ha,
                                 color=LC_colours[j], width=LCarea_width, bottom=basept, 
                                 label=LCname+' area'  ) )
            basept = basept+LC_DF[LCareaname][year]*m2_to_ha
            
    ax.set_xticks(np.arange(0.5,nyears,1.))
    ax.set_xticklabels(years)
    xlims = ax.get_xlim()
    ax.plot(xlims,[0,0],c='k',lw=0.5,ls='--')
    ax.set_xlim(xlims)
    ax.set_ylabel(Species_Names[iSp]+' (Mt per year)',fontsize=14)
    ax2.set_ylabel('Land Cover Areal (Ha)',fontsize=14)
fig.legend(lines,Legend_labels,fontsize=13,ncol=nLCs+2,loc=8)
fig.savefig(data_dir+'ScaledUp_Sabah_Soil_Emissions.png',bbox_inches='tight')


LC_names = ['Forest', 'Oil Palm', 'Total' ]
LCarea_names = ['All Forest (adjusted)', 'Oil Palm' ]
Legend_labels = [LC+' emissions' for LC in LC_names]+['Forest Area', 'Oil Palm Area' ]
LC_hatches  = ['//', "\\\\", None ]  
LC_fcolours = ['lightgrey', 'darkgrey', 'dimgrey' ]
LC_xpos_ref = [-0.2, -0.05, 0.15]  
LC_boxwidth = [ 0.1, 0.1, 0.2]
LCarea_xpos_ref = 0.35
LCarea_width = 0.05
nLCs = len(LC_names)
nyears = len(years)

fig,axes = plt.subplots(nrows=2,ncols=1,figsize=[12,8])
Species_Names = ['N$_2$O', 'CH$_4$']
Species_DFs = [N2O_DF,CH4_DF]
nSpecies = len(Species_Names)
for iSp in range(nSpecies):
    ax = axes[iSp]
    # second axis Land Cover area 
    ax2=ax.twinx()
    DF=Species_DFs[iSp]
    for i in range(nyears):
        year=years[i]
        xpos_centre = i+0.5
        lines=[]
        for j in range(nLCs):
            LCname=LC_names[j]
            lines.append(box_and_whisker(DF[LCname][year], xpos_centre+LC_xpos_ref[j], 
                                         'k', ax, boxwidth=LC_boxwidth[j], 
                                         fill=True, fcolor=LC_fcolours[j])) #,hatch=LC_hatches[j]))
        basept = 0.0
        for j in range(nLCs-1):
            LCname=LC_names[j]
            LCareaname=LCarea_names[j]
            lines.append(ax2.bar(xpos_centre+LCarea_xpos_ref, LC_DF[LCareaname][year]*m2_to_ha,
                                 color=LC_fcolours[j], hatch=LC_hatches[j], edgecolor='k',
                                 width=LCarea_width, bottom=basept ) ) 
            basept = basept+LC_DF[LCareaname][year]*m2_to_ha
            
    ax.set_xticks(np.arange(0.5,nyears,1.))
    ax.set_xticklabels(years)
    xlims = ax.get_xlim()
    ax.plot(xlims,[0,0],c='k',lw=0.5,ls='--')
    ax.set_xlim(xlims)
    ax.set_ylabel(Species_Names[iSp]+' (Mt per year)',fontsize=14)
    ax2.set_ylabel('Land Cover Areal (Ha)',fontsize=14)

fig.legend(lines,Legend_labels,fontsize=13,ncol=nLCs+2,loc=8)
fig.savefig(data_dir+'ScaledUp_Sabah_Soil_Emissions_Greyscale.png',bbox_inches='tight')




