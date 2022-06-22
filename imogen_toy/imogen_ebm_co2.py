# This is a climate Energy Balance Model (EBM) based around Huntingford and Cox   
# "An analogue model to derive additional climate change scenarios from existing GCM simulations" 
# Climate Dynamics, 2000. Vol 16, pp575-586.
# C. Huntingford (16th July 2015)

# READS IN CO2 concentrations instead.

# This needs to run on CEH linux box wllf017 or higher (to get newest version of matrix inversion). 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys,os
from scipy import sparse
from scipy.sparse.linalg import spsolve 
from scipy.interpolate import interp1d

import ocean_co2

# Set oceanic effective diffusion, clim sens over land, ocean and land differential warming rate
kappa = 474.0                                  # (W/m/K)
lambda_l = 1.486                               # (W/m2/K)
lambda_o = 1.227                               # (W/m2/K)
mu = 1.513                                     # Unitless

# If want second run with modified forcings, n_change=2, otherwise n_change=1
n_change=1

# Set volumetric heat capacity for salt water, and fraction of Earth with ocean
cp = 4.04E6                                    # (J/K/m3)
f = 0.7066                                     # Unitless
ocean_area = 3.627E14                          # (m2)

# Set initial CO2 concentration, oceanic temperature and "memory" of CO2-ocean interactions
co2_atmos_ppmv_init = 284.725                  # (ppm)
conv=0.471                                     # Converts (net) emissions GtC to change in ppm (ppm/GtC)
t_ocean_init = 289.28                          # (K) 

# Set up arrays to collect for plotting. Assumming that the emissions are affecting year after (in timestep), then yrs are 1861-2100 inclusive. 
nyr = 241                                      # Number of years of simulation 
nyr_plot=nyr-1
yr_plot=np.zeros(nyr_plot)
yr_plot_drive=np.zeros(nyr_plot)
emiss_plot=np.zeros(nyr_plot)
non_co2_rf_plot=np.zeros(nyr_plot)
co2_plot=np.zeros(nyr_plot)
temp_global_plot=np.zeros(nyr_plot)            # Holds the global mean temperature change (K)
if (n_change==2):
   yr_plot_change=np.zeros(nyr_plot)
   yr_plot_drive_change=np.zeros(nyr_plot)
   emiss_plot_change=np.zeros(nyr_plot)
   non_co2_rf_plot_change=np.zeros(nyr_plot)
   co2_plot_change=np.zeros(nyr_plot)
   temp_global_plot_change=np.zeros(nyr_plot)            # Holds the global mean temperature change (K)

for i_change in range(0,n_change):
   fa_ocean=np.zeros(10000)

# Now read in both "business-as-usual" emissions CO2 emissions.
# The timestepping is that temp is calculated 1860,.....,2100 
# Effects of emissions on CO2 cycle (and so Q) are updated after 1860.....after 2099
   emiss_in = np.zeros(nyr-1)                     # Array holds emissions in (GtC/yr)
   f_in = open('emiss_co2.dat', 'r')
   for iyr in range(0,nyr-1): 
      nums_in_string=f_in.readline()
      emiss_in[iyr]= np.float(nums_in_string.split()[1]) 
   f_in.close()

# Now add in non-CO2 radiative forcing
   q_non_co2 = np.zeros(nyr-1)                    # Array holds non-CO2 radiative forcing (W/m2)
   f_in = open('/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/SSP2-2.6_IMAGE_qnonco2_vn2p0.txt', 'r')
   for iyr in range(0,nyr-1): 
      nums_in_string=f_in.readline()
      q_non_co2[iyr]= np.float(nums_in_string.split()[1]) 
   f_in.close()
#   q_non_co2 = q_non_co2 -q_non_co2[0]

# Now add in CO2 driver instead 
   co2_input = np.zeros(nyr-1)                      # Array holds non-CO2 radiative forcing (W/m2)
   q_co2_input = np.zeros(nyr-1)                    # Array holds non-CO2 radiative forcing (W/m2)
   #f_in = open('/users/global/chg/imogen_test/trans_cen_BCC_mod_bcc-csm1-1/CO2.dat', 'r')
   f_in = open('/prj/CLIFFTOP/COMMON_DATA/SCENARIOS/SSP2-2.6_IMAGE_concs_co2_vn2p0.txt','r')
   for iyr in range(0,nyr-1): 
      nums_in_string=f_in.readline()
      co2_input[iyr]= np.float(nums_in_string.split()[1]) 
      q_co2_input[iyr] = 5.397*np.log(co2_input[iyr]/co2_input[0])
   #print('co2',co2_input)
   #print('q_co2',q_co2_input)
   f_in.close()

# User change here emissions and non-CO2 RF (if i_change=1)
# Now do spline from 2020 onwards. Force it to start at 2020, go through decade mid-points
# Original value below - option to change underneath 

# Emiss 2020-2030 = 13.5 Gt/yr  RF non-co2 2020-2030 = 1.06 W/m2 
# Emiss 2030-2040 = 15.4 Gt/yr  RF non-co2 2030-2040 = 1.14 W/m2 
# Emiss 2040-2050 = 16.7 Gt/yr  RF non-co2 2040-2050 = 1.24 W/m2 
# Emiss 2050-2060 = 18.3 Gt/yr  RF non-co2 2050-2060 = 1.34 W/m2 
# Emiss 2060-2070 = 20.0 Gt/yr  RF non-co2 2060-2070 = 1.43 W/m2 
# Emiss 2070-2080 = 22.0 Gt/yr  RF non-co2 2070-2080 = 1.50 W/m2 
# Emiss 2080-2090 = 24.7 Gt/yr  RF non-co2 2080-2090 = 1.56 W/m2 
# Emiss 2090-2100 = 27.6 Gt/yr  RF non-co2 2090-2100 = 1.62 W/m2 

# User supplied numbers here:
   if (i_change==1):
      x_vals=[2020, 2025, 2035, 2045, 2055, 2065, 2075, 2085, 2095, 2099]
      new_2025_co2 = 12.0
      new_2035_co2 = 10.0
      new_2045_co2 = 9.0 
      new_2055_co2 = 8.0 
      new_2065_co2 = 4.0 
      new_2075_co2 = 3.0 
      new_2085_co2 = 2.5 
      new_2095_co2 = 1.0 
      y_vals_co2=[emiss_in[160], new_2025_co2, new_2035_co2, new_2045_co2, new_2055_co2, new_2065_co2, 
                                        new_2075_co2, new_2085_co2, new_2095_co2, new_2095_co2]
      f_inter_function_co2 = interp1d(x_vals, y_vals_co2, kind='cubic')
      emiss_new_from_2020=f_inter_function_co2(np.linspace(2020,2099,80))
      emiss_in[160:nyr-1]=emiss_new_from_2020[:]

      new_2025_rf = 1.0
      new_2035_rf = 0.8
      new_2045_rf = 0.5 
      new_2055_rf = 0.5 
      new_2065_rf = 0.3 
      new_2075_rf = 0.25 
      new_2085_rf = 0.1 
      new_2095_rf = 0.1 
      y_vals_rf=[q_non_co2[160], new_2025_rf, new_2035_rf, new_2045_rf, new_2055_rf, new_2065_rf, 
                                        new_2075_rf, new_2085_rf, new_2095_rf, new_2095_rf]
      f_inter_function_rf = interp1d(x_vals, y_vals_rf, kind='cubic')
      rf_new_from_2020=f_inter_function_rf(np.linspace(2020,2099,80))
      q_non_co2[160:nyr-1]=rf_new_from_2020[:]

# Initialise atmospheric CO2 concentration, and increment
   co2_atmos_ppmv = co2_atmos_ppmv_init
   co2_atmos_change_ppmv = 0.0

# Now set up the conditions for solution to the parabolic diffusion equation. 
# Set number of calls of ocean diffusion PDE solver per year.
   n_pde = 20
   dt = (1.0/float(n_pde))*60.0*60.0*24.0*365.0   # Model timestep (seconds)
   temp_ocean_top = np.zeros(nyr)                 # Array holds oceanic warming (K)

# Set number of vertical layers and ocean depth 
   n_vert = 254 
   depth = 5000.0                                 # (metres)
   dz = depth/float(n_vert)                       # (metres)

   s=(kappa/cp)*(dt/(dz*dz))
# Notation. Jankowska uses (i) h for dz, (ii) k for dt, (iii) alpha^2 for (kappa/cp)
# Jankowska uses (iv) lambda for s. Their lambda=alpha^2*(k/h^2), so gives our "s"

# Now do the pde solving. First set up the arrays (and zero'd)
   t_ocean_old=np.zeros(n_vert+1)
   t_ocean_new=np.zeros(n_vert+1)

#######################################################
# The numerical scheme is below (follow Jankowska "An interval finite difference method of Crank-Nicolson
# type for solving the one-dimensional heat conduction equation with mixed boundary conditions
# 
# -(s/2)u_{i-1,j+1} + (1+s)u_{i,j+1} - (s/2)u_{u+1,j+1} = (s/2)u_{i-1,j} + (1-s)u_{i,j} + (s/2)u_{u+1,j}
# for timesteps j=1,2,..... and depth i=0,1,.......,n
# 
# Then boundary conditions are for the imaginary point at bottom of ocean u_{n+1,j}=u_{n-1,j} (and similarly for j+1) 
# And at the top of the ocean is u_{-1,j}=u_{1,j} = 2*dz*(N/(f*kappa)) [**CHECK - don't think did this in the end**]
#
# Note - assume u idential to zero at t=0. (No need to worry about Eqn(31))
# Set up Eqn to solve CU_{j+1}=DU_{j}+E where C, D are [n+1 x n+1] and U is [n+1, 1] and E is [n+1, 1] 
#######################################################

   C=np.zeros([n_vert+1, n_vert+1])
   D=np.zeros([n_vert+1, n_vert+1])
   E=np.zeros(n_vert+1)

# Conservation check for end of run - add up energy in. 
   q_energy=0.0

# Now start loop over the different years. 
   for j in range(0,nyr):

# Update forcing before iterating over sub-years in the PDE solver.
      q2co2=3.74
      #q = (q2co2/np.log(2.0)) * np.log(co2_atmos_ppmv/co2_atmos_ppmv_init)
      q = q_co2_input[j] + q_non_co2[j] 
      print('q, q_co2_input[j], q_non_co2[j]')
      print(q, q_co2_input[j], q_non_co2[j])
# Add in previous year's non-co2 Q
      #if (j >= 1):
#     #    q = q + q_non_co2[j-1]
      #   q = q_co2_input[j-1] + q_non_co2[j-1]

# Two factors on rhs of Eqn (10), Huntingford and Cox (2000). These are mixed boundary conditions.
# Writing with division by kappa and sign change, then our Eqn (10) becomes:
# dT/dz = factor1 + factor2 * T(z=0)
# Hence in Jankowska notation, their psi_1 is our factor 1, their A is our factor 2.
      factor1=-q/(kappa*f)
      if (j==0):
         temp_ocean_top_local = 0.0
      else:
         temp_ocean_top_local = temp_ocean_top[j-1]
      factor2=((1.0-f)*lambda_l*mu)/(f*kappa)  + (lambda_o/kappa)
      print('factor1,factor2=',factor1,factor2)

      for k in range(0,n_pde):
         ival=j*n_pde+k
# First set t_old as t_new ready for next iteration
         t_ocean_old=t_ocean_new
# Sort out the top points, bottom points and then all for C  
#      C[0, 0] = (1.0+s) 
         C[0, 0] = s*(1.0+dz*factor2)+1 
         C[0, 1] = -s
         C[n_vert, n_vert-1] = -s
         C[n_vert, n_vert] = (1.0+s)
         for m in range(1, n_vert):
            C[m, m-1] = -s/2.0
            C[m, m+1] = -s/2.0
            C[m, m] = 1.0 + s 
# Sort out the top points, bottom points and then all for D  
#      D[0, 0] = (1.0-s) 
         D[0, 0] = -s*(1.0+dz*factor2)+1 
         D[0, 1] = s
         D[n_vert, n_vert] = (1.0-s)
         D[n_vert, n_vert-1] = s
         for m in range(1, n_vert):
            D[m, m-1] = s/2.0
            D[m, m+1] = s/2.0
            D[m, m] = 1.0 - s 
#      E[0]=((s*dz)/(f*kappa))*(toa_smoothed_norm_sub[ival]+toa_smoothed_norm_sub[ival+1])
#      E[0]=-(2.0*s*dz)*factor1     # Assume slow variation in Q
         E[0]=-(kappa/cp)*(dt/dz)*(factor1+factor1)     # Assume slow variation in Q

# Now solve for U_{j+1}. First put bits in to tri-diagonal things for "sparse.spdiags".
         C_sub=np.zeros(n_vert+1)
         C_main=np.zeros(n_vert+1)
         C_super=np.zeros(n_vert+1)
         for m in range(0,n_vert+1): 
             C_main[m]=C[m,m]
         for m in range(0,n_vert): 
             C_sub[m]=C[m+1,m]
             C_super[m]=C[m,m+1]

# Now calculate the right-hand side (called b_rhs)
         d_mat=np.mat(D)
         e_mat=np.mat(E).transpose()
         t_ocean_old_mat=np.mat(t_ocean_old).transpose()
         b_rhs=np.dot(D,t_ocean_old_mat)+e_mat
         b_rhs=np.ravel(b_rhs)

# Now perform the calculation to update the oceanic temperatures 
         t_ocean_new=spsolve(C,b_rhs)
         print(t_ocean_new[:4])
# Update the cumulative heat entering ocean, and also save to top ocean temperature
#      q_total=q_total+dt*(toa_smoothed_norm_sub[ival]/f)
         q_energy=q_energy+dt*kappa*(-factor1-(t_ocean_new[0]*factor2))

      temp_ocean_top[j]=t_ocean_new[0]
# Now at end of year, update the CO2 concentration. Don't update at end of last year temperature is calculated.
# For vegetation, temporarily assume draws down at all times 25% of emissions.
      if (j <= nyr-2):
         co2_previous = co2_atmos_ppmv
         co2_atmos_ppmv = co2_atmos_ppmv + conv*emiss_in[j]
         d_land_atmos = -conv*0.25*emiss_in[j]
         co2_atmos_ppmv = co2_atmos_ppmv + d_land_atmos

# Now calculate the oceanic draw-down
# At end of year, now update atmospheric CO2 concentration and radiative forcing
         nfarray=10000                                   # Number of timesteps of history of CO2 draw-down for Greens fn
#         d_ocean_atmos = ocean_co2.ocean_co2(j, co2_atmos_ppmv, co2_atmos_ppmv_init, t_ocean_new[0],
#            fa_ocean,ocean_area,co2_atmos_change_ppmv,nyr,t_ocean_init,nfarray)
         d_ocean_atmos = ocean_co2.ocean_co2(j, co2_input[j], co2_atmos_ppmv_init, t_ocean_new[0],
            fa_ocean,ocean_area,co2_input[j]-co2_input[j-1],nyr,t_ocean_init,nfarray)

         co2_atmos_ppmv = co2_atmos_ppmv + d_ocean_atmos
         co2_atmos_change_ppmv = co2_atmos_ppmv-co2_previous

#         print 'year, CO2 concs, temp_0, total Q = ', j+1860, co2_atmos_ppmv, t_ocean_new[0], q
# Collect data for plotting 1861-2100 
      if (j > 0):
         if (i_change == 0):
            yr_plot[j-1]=j+1860
            co2_plot[j-1]=co2_input[j-1]
            temp_global_plot[j-1]=f*t_ocean_new[0] + (1-f)*mu*t_ocean_new[0]
            yr_plot_drive[j-1]=j+1860-1
            emiss_plot[j-1]=emiss_in[j-1]
            non_co2_rf_plot[j-1]=q_non_co2[j-1]
         else: 
            yr_plot_change[j-1]=j+1860
            co2_plot_change[j-1]=co2_input[j-1]
            temp_global_plot_change[j-1]=f*t_ocean_new[0] + (1-f)*mu*t_ocean_new[0]
            yr_plot_drive_change[j-1]=j+1860-1
            emiss_plot_change[j-1]=emiss_in[j-1]
            non_co2_rf_plot_change[j-1]=q_non_co2[j-1]

      if j==3: quit()
      print 'year, CO2 concs, emiss, d_land_atmos, d_ocean_atmos, temp'
      print j+1850, co2_input[j], conv*emiss_in[j], d_land_atmos, d_ocean_atmos, t_ocean_new[0:5]
      print 'q,q_co2_input[j], q_non_co2[j]:'
      print q,q_co2_input[j], q_non_co2[j] 

# Check conservation of energy at end of run. Compare against total heat in the final profile.
   q_energy_derived=0.0
   for l in range(1,n_vert+1):
      q_energy_derived=q_energy_derived+cp*0.5*(t_ocean_new[l]+t_ocean_new[l-1])*dz 

   print 'Heat conservation check (%) = ',100.0*(q_energy_derived/q_energy)

print 'final temperature = ',temp_global_plot[-1]
# Now create the plot
fig = plt.figure(figsize=(6.0,8.0))
mpl.rc('xtick',labelsize=9) ; mpl.rc('ytick',labelsize=9)
# Plot driver of emissions
ax1=plt.subplot(4,1,1)
ax1.plot(yr_plot_drive, emiss_plot,linewidth=2.0,color='blue')
if(n_change == 2):
   ax1.plot(yr_plot_drive_change, emiss_plot_change,linewidth=2.0,color='blue',linestyle='--')
ax1.set_title('Carbon Emissions',fontsize=11)
ax1.set_ylabel('(Gt / Year)', fontsize=10)
ax1.set_ylim(0.0,30.0)
ax1.text(1855,25,'Continuous IS92a, business-as-usual',fontsize=10)
ax1.text(1855,20,'Dashed, user supplied post 2020',fontsize=10)

# Plot driver of non-CO2 RF 
ax2=plt.subplot(4,1,2)
ax2.plot(yr_plot_drive, non_co2_rf_plot,linewidth=2.0,color='blue')
if(n_change == 2):
   ax2.plot(yr_plot_drive_change, non_co2_rf_plot_change,linewidth=2.0,color='blue',linestyle='--')
ax2.set_title('Non-CO2 Radiative Forcing (Methane etc)',fontsize=11)
ax2.set_ylabel('(W / m2)', fontsize=10)
ax2.set_ylim(0.0,1.8)

# Plot derived atmospheric CO2 concentration 
ax3=plt.subplot(4,1,3)
ax3.plot(yr_plot, co2_plot,linewidth=2.0,color='red')
if(n_change == 2):
   ax3.plot(yr_plot_change, co2_plot_change,linewidth=2.0,color='red',linestyle='--')
ax3.set_title('Calculated CO2 concentration',fontsize=11)
ax3.set_ylabel('(ppm)',fontsize=10)
ax3.set_ylim(200.0,900.0)

# Plot driver of non-CO2 RF 
ax4=plt.subplot(4,1,4)
ax4.plot(yr_plot, temp_global_plot,linewidth=2.0,color='red')
if(n_change == 2):
   ax4.plot(yr_plot_change, temp_global_plot_change,linewidth=2.0,color='red',linestyle='--')
ax4.set_title('Calculated global warming',fontsize=11)
ax4.set_ylabel('(K)',fontsize=10)
ax4.set_xlabel('Year', fontsize=11)
ax4.plot([1850.0,2100.0],[2.0,2.0],color='black',linestyle='--',linewidth=1.0)
ax4.set_ylim(0.0,5.0)

plt.subplots_adjust(bottom=0.08, top=0.95, left=0.12, right=0.94, hspace=0.40)
plt.savefig('imogen_ebm.pdf')
print 'Got to end'
sys.exit()
