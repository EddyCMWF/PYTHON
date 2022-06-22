#! /bin/python
#
# Script to run JULES for the various PALS sites.
#
##################################################

import os, sys, glob
import csv
import numpy as np
import datetime as dt
import netCDF4 as nc

# Search for optional parsed arguments and pop them from list of args
if '-prespun' in sys.argv:
    prespun='Y'
    argloc = sys.argv.index('-prespun')
    temp=sys.argv.pop(argloc)
    del temp
    del argloc
else:
    prespun='N'

# Search for optional parsed arguments and pop them from list of args
if '-respin' in sys.argv:
    respin='Y'
    argloc = sys.argv.index('-respin')
    temp=sys.argv.pop(argloc)
    del temp
    del argloc
else:
    respin='N'

Default_soil_file='/users/eow/edwcom/GREENHOUSE/PALS_comparison/default_soil_params.csv'
if '-VG' in sys.argv:
    soil_name='VanGenuchten'
    soil_file='/users/eow/edwcom/WFD_EI/SoilProperties_wfdei_VanGenuchten_2D.nc'
    argloc= sys.argv.index('-VG')
    temp=sys.argv.pop(argloc)
    del temp
    del argloc
elif '-BC' in sys.argv:
    soil_name='BC'
    soil_file='/users/eow/edwcom/WFD_EI/SoilProperties_wfdei_BrookesCorey_2D.nc'
    argloc= sys.argv.index('-BC')
    temp=sys.argv.pop(argloc)
    del temp
    del argloc
else:
    soil_name='Default'
    soil_file=Default_soil_file



# Read parsed run option
run = sys.argv[1]
print prespun, respin, run
# Accepted options: 'Jvn4.3.1' standard jules at vn4.3.1
#                   'Jvn4.3.1-E-F' standard jules-ecosse-fun @
PALS_DIR  = '/prj/GREENHOUSE/PALS_comparison/'

WORK_DIR  = PALS_DIR+run+'/'
if not os.path.isdir(WORK_DIR):
    os.mkdir(WORK_DIR)

DATA_DIR  = PALS_DIR+'met_data/'
if not os.path.isdir(DATA_DIR):
    os.mkdir(DATA_DIR)

output_dir = WORK_DIR+'output/'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

# Select Jules executable based on run selected
if run=='Jvn4.3.1':
    JULES_EXECUTABLE='/users/eow/edwcom/JULES/jules_vn4.3.1_ECP/exes/jules_std.exe'
elif 'Jvn4.3.1-E' in run:
    JULES_EXECUTABLE='/users/eow/edwcom/JULES/r2806_vn4.3.1_ecosse_fun_ECP/exes/jules_ecosse_fun_std.exe'
elif 'Jvn4.5' in run:
    JULES_EXECUTABLE='/users/eow/edwcom/JULES/jules_vn4.5_trunk/exes/jules.exe'
else: 
    JULES_EXECUTABLE='/users/eow/edwcom/JULES/jules_vn4.3.1_ECP/exes/jules_std.exe'
    print 'run not recognised, using jules executable located: '
    print JULES_EXECUTABLE

N_Deposition_file = PALS_DIR+'ancil_data/deposition_N_annual_1g_year.dat'

site_metadata_file='sites_meta.csv'
site_md_list = list(csv.reader(open(PALS_DIR+site_metadata_file,'r'),delimiter=','))
site_md_headers=site_md_list.pop(0)

site_frac_file = 'sites_julesfracs.dat'
site_frac_dat = list(csv.reader(open(PALS_DIR+site_frac_file,'r'),delimiter='-'))
frac_list = [line[1] for line in site_frac_dat]

os.chdir(WORK_DIR+'work_namelists/')

# Open Soil file and create 2D array of soil properties
lines=open(Default_soil_file,'r').readlines()
Soil_vars=lines[0].split(',')[:-1]
Soil_vals=[float(val) for val in lines[1].split(',')[:-1]]
Soil_Values_line=lines[1][:-1]
CS='5'

if soil_name!='Default':
    Soil_Dict={}
    inf=nc.Dataset(soil_file,'r')
    for var in Soil_vars:
        Soil_Dict[var]=inf.variables[var][:]
    CS_grid=inf.variables['cs'][:]
    inf.close()


#for site_md,frac in zip(site_md_list,frac_list):
for site_md,frac in zip(site_md_list[19:],frac_list[19:]):
    # strip relevant meta data of blank spaces
    site       = site_md[0].strip()
    VEG_type   = site_md[2].strip()
    canht      = site_md[3].strip()
    Start_date = site_md[5].strip()
    End_date   = site_md[6].strip()
    latitude   = site_md[7].strip()
    longitude  = site_md[8].strip()
    jules_longitude = '0.0'   # JULES fudge to account for local time met data
    
    if soil_name!='Default':
        lat_index=int((float(latitude)+90.)*2)
        lon_index=int((float(longitude)+180.)*2)
        Soil_vals=[str(Soil_Dict[var][lat_index,lon_index]) for var in Soil_vars]
        Soil_Values_line=''
        for val in Soil_vals:
            Soil_Values_line+=val+','
        CS=str(CS_grid[lat_index,lon_index])

    print 'SoilVals: '+Soil_Values_line
    # Create datetime object of end_date and subtract a day
    End_date_obj = dt.datetime.strptime(End_date, '%Y-%m-%d %H:%M:%S') \
                    -dt.timedelta(days=1)
    # Update string value
    End_run_date = str(End_date_obj)
    
    # Create spin end day by adding on a year
    Spin_end_date= str( dt.datetime.strptime(Start_date, '%Y-%m-%d %H:%M:%S') \
                            + dt.timedelta(days=365) )
    
    # Create Drive Filename from metadata:
    drive_file = DATA_DIR+site+'Fluxnet.1.4_met.nc'
    
    # Create Output settings
    run_id     = site
    print run_id
    if 'variableN' in run:
        N_Deposition_file=PALS_DIR+'ancil_data/site_depostion_files/' \
                            + site+'_Ndeposition_mean.dat'

    # Copy template namelists into working directory
    os.system('cp '+WORK_DIR+'template_namelists/*.nml '\
                   +WORK_DIR+'work_namelists/'          )
    
    # If prespun then copy the modified initial conditions and timesteps nmls
    if (prespun=='Y'):
        os.system('cp '+WORK_DIR+'template_namelists/prespun/*.nml '\
                  +WORK_DIR+'work_namelists/'          )
        prespun_ic_file=glob.glob(output_dir+site+'.dump.'+\
                             Start_date[:4]+Start_date[5:7]+Start_date[8:10] +\
                             '*.nc')[0]
        os.system('sed -i s?START_DUMP_DUMMY?"'+prespun_ic_file+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    if (respin=='Y'):
        os.system('cp '+WORK_DIR+'template_namelists/prespun/initial_conditions.nml '\
                  +WORK_DIR+'work_namelists/'          )
        respin_ic_file=glob.glob(output_dir+site+'.dump.'+\
                             Start_date[:4]+Start_date[5:7]+Start_date[8:10] +\
                             '*.nc')[0]
        os.system('sed -i s?START_DUMP_DUMMY?"'+respin_ic_file+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    # sed the metadata into the namelists
    os.system('sed -i s?START_DATE_DUMMY?"'+Start_date+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?START_SPIN_DATE_DUMMY?"'+Start_date+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    os.system('sed -i s?END_DATE_DUMMY?"'+End_date+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?END_SPIN_DATE_DUMMY?"'+Spin_end_date+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?END_RUN_DATE_DUMMY?"'+End_run_date+'"?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    
    os.system('sed -i s?LATITUDE_DUMMY?'+latitude+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?LONGITUDE_DUMMY?'+jules_longitude+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    os.system('sed -i s?DRIVE_FILE_DUMMY?'+drive_file+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    os.system('sed -i s?OUTPUT_DIR_DUMMY?'+output_dir+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?RUN_ID_DUMMY?'+run_id+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i s?WORK_DIR_DUMMY?'+WORK_DIR+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    
    os.system('sed -i s?N_DEPSOTION_FILE_DUMMY?'+N_Deposition_file+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    print 'sed -i s?SOIL_PROPERTIES_DUMMY?"'+Soil_Values_line+'"?g '\
                           + WORK_DIR+'work_namelists/*nml'
    os.system('sed -i s?SOIL_PROPERTIES_DUMMY?'+Soil_Values_line+'?g '\
               + WORK_DIR+'work_namelists/*nml')
    os.system('sed -i  s?SOIL_CARBON_DUMMY?"'+CS+'"?g '\
                + WORK_DIR+'work_namelists/*nml' )
    
    # Copy template initial conditions
    os.system('cp '+WORK_DIR+'initial_conditions.dat.template '\
                   +WORK_DIR+'initial_conditions.dat'           )
    
    # sed replace dummy string with with frac
    os.system('sed -i  s?FRACTIONAL_COVER_DUMMY?"'+frac+'"?g '\
                + WORK_DIR+'initial_conditions.dat' )
    
    if (prespun!='Y'):
        os.system('tar -czf archive_nmls/'+site+'_nmls.tgz *.nml')
    #quit()
    print 'Executing: '
    print JULES_EXECUTABLE+' > z_logs/'+site+'_run.log 2> z_logs/'+site+'_err.log'
    print 'From: '
    print os.getcwd()
    os.system(JULES_EXECUTABLE+' > z_logs/'+site+'_run.log 2> z_logs/'+site+'_err.log')
    
    #quit()
    os.system('rm *.nml')
    os.system('rm '+output_dir+site+'.dump.spin*.nc')
    

