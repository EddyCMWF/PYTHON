#!/usr/bin/python
#
# Script to compare JULES-ecosse output with standar jules
#
#

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import netcdftime as nctime


albmar_file='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/albmar_crich_operative.day.nc'
#ecosse_file='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/test_namelists/output/JULES_ECOSSE_TEST_RUNS.day.nc'
#ecosse_file2='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/test_namelists/output/JULES_ECOSSE_TEST_RUNS.day_ecosse.nc'

ecosse_file='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/with_soilN/output/JULES_ECOSSE_with_soilN.day.nc'
ecosse_file2='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/with_soilN/output/JULES_ECOSSE_with_soilN.day_ecosse.nc'

output_dir='/users/eow/edwcom/GREENHOUSE/jules_v4.3.1_ecosse_site_runs/plots/'

alb_inf = nc.Dataset(albmar_file,'r')
TIME=alb_inf.variables['time'][:]
TIME_obj = nctime.num2date(TIME, \
                           units=alb_inf.variables['time'].units, \
                           calendar=alb_inf.variables['time'].calendar)

eco_inf = nc.Dataset(ecosse_file,'r')
eco_inf2 = nc.Dataset(ecosse_file2,'r')

var_list = [ str(var) for var in alb_inf.variables ]
for i in [-1,12,11,10,9,8,7,6,3,2,1,0]:
     var_list.pop(i)

FIG=plt.figure(figsize=[19,11])

for i in range(len(var_list)):
    var=var_list[i]
    print var, np.mean(alb_inf.variables[var][:]-eco_inf.variables[var][:])
    AX=FIG.add_subplot(6,4,i+1)
    if (var=='c_veg') | (var=='lai') | (var=='canht'):
        albline=AX.plot(TIME_obj,alb_inf.variables[var][:,2,0].squeeze(),label='albmar-run')
        ecoline=AX.plot(TIME_obj,eco_inf.variables[var][:,2,0].squeeze(),label='ecosse-run')
    elif var=='resp_s':
        albmar_resp_s = np.sum(alb_inf.variables[var][:].squeeze(),axis=1)
        albline=AX.plot(TIME_obj,albmar_resp_s,label='albmar-run')
        ecoline=AX.plot(TIME_obj,eco_inf2.variables['co2_soil_gb'][:].squeeze(),label='ecosse-run')
    elif var=='cs':
        albline=AX.plot(TIME_obj,alb_inf.variables[var][:,0].squeeze(),label='albmar-run')
        ecoline=AX.plot(TIME_obj,eco_inf2.variables['c_dpm_total_gb'][:].squeeze(),label='ecosse-run')
    else:
        albline=AX.plot(TIME_obj,alb_inf.variables[var][:,0].squeeze(),label='albmar-run')
        ecoline=AX.plot(TIME_obj,eco_inf.variables[var][:,0].squeeze(),label='ecosse-run')
    
    AX.set_title(var)
    ticklocs = AX.xaxis.get_majorticklocs()
    AX.xaxis.set_ticks( ticklocs[::2] )

AX.legend( bbox_to_anchor=(1.4,0.8),loc=2,borderaxespad=0.)
FIG.tight_layout(rect=[0,0,1,0.94])
FIG.suptitle('Comp. of stndrd JULES and JULES-ECOSSE at Crichton fieldsite', \
             fontsize=30)
FIG.savefig(output_dir+'Time-Series-Comparison-TESTRUN.png', \
            bbox_inches='tight')

plt.close()


#Plot the CS pools
FIG=plt.figure(figsize=[12,10])

#DPM
AX=FIG.add_subplot(2,2,1)
albline=AX.plot(TIME_obj,alb_inf.variables['cs'][:,0].squeeze(),label='albmar-run')
ecoline=AX.plot(TIME_obj,eco_inf2.variables['c_dpm_total_gb'][:].squeeze(),label='ecosse-run')
AX.set_title('DPM')
AX.set_ybound(lower=0)
ticklocs = AX.xaxis.get_majorticklocs()
AX.xaxis.set_ticks( ticklocs[::2] )

#RPM
AX=FIG.add_subplot(2,2,2)
albline=AX.plot(TIME_obj,alb_inf.variables['cs'][:,1].squeeze(),label='albmar-run')
ecoline=AX.plot(TIME_obj,eco_inf2.variables['c_rpm_total_gb'][:].squeeze(),label='ecosse-run')
AX.set_title('RPM')
AX.set_ybound(lower=0)
ticklocs = AX.xaxis.get_majorticklocs()
AX.xaxis.set_ticks( ticklocs[::2] )

#BIO
AX=FIG.add_subplot(2,2,3)
albline=AX.plot(TIME_obj,alb_inf.variables['cs'][:,2].squeeze(),label='albmar-run')
ecoline=AX.plot(TIME_obj,eco_inf2.variables['c_bio_total_gb'][:].squeeze(),label='ecosse-run')
AX.set_title('Bio')
AX.set_ybound(lower=0)
ticklocs = AX.xaxis.get_majorticklocs()
AX.xaxis.set_ticks( ticklocs[::2] )

#HUM
AX=FIG.add_subplot(2,2,4)
albline=AX.plot(TIME_obj,alb_inf.variables['cs'][:,3].squeeze(),label='albmar-run')
ecoline=AX.plot(TIME_obj,eco_inf2.variables['c_hum_total_gb'][:].squeeze(),label='ecosse-run')
AX.set_title('Hum')
AX.set_ybound(lower=0)
ticklocs = AX.xaxis.get_majorticklocs()
AX.xaxis.set_ticks( ticklocs[::2] )

AX.legend( bbox_to_anchor=(-0.1,-0.2),loc=10,borderaxespad=0.,ncol=2)
FIG.suptitle('Comp. of soil CS from JULES and JULES-ECOSSE at Crichton',\
              fontsize=30 )
FIG.tight_layout(rect=[0,0.06,1,0.94])

FIG.savefig(output_dir+'CS-Time-Series-Comparison-TESTRUN.png', \
            bbox_inches='tight')

plt.close()





