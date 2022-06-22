#!/bin/env python2.7

import sys
import netCDF4 as nc
import numpy as np

def Dif(mix, z, depthuni, botdepth):
    if (z<=depthuni):
        dif=mix
    elif (z<botdepth):
        dif = mix*(1-(z-depthuni)/(botdepth-depthuni))
    else:
        dif = 1e-11
    return dif

def Rin(tauin, z, litter):
    return tauin*litter*np.exp(-tauin*z)

if len(sys.argv) < 3:
    print('ERROR: Three input filenames required: ./equilibrium_soilcarbon.py [INFILE] [DUMPFILE] [OUTFILE]')
    print('        optional inputs: -l_update_dumpfile, to update the soil carbon in the input dumpfile')
    quit()

if '-l_update_dumpfile' in sys.argv:
    L_update_dumpfile=True
else:
    L_update_dumpfile=False

infile=sys.argv[1]
Dinfile=sys.argv[2]
outfile=sys.argv[3]

outvars = ['cs','resp_s','lit_c','frac','t_soil','clay']

nland = 1631

dpm_rpm_ratio=(0.25,0.25,0.25,0.25,0.25,0.67,1.44,1.44,0.67,1.44,1.44,0.33,0.33)
npft=len(dpm_rpm_ratio)

dz = np.array([0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293,0.21517585,
                   0.23784142,0.25980762,0.28117066,0.30200527,0.32237098,0.34231625,0.36188121])
nz = len(dz)
z = np.zeros_like(dz)
for iz in range(1,nz):
    z[iz] = z[iz-1]+dz[iz-1]

kaps       = [10.,0.3,0.66,0.02]  # units of per year, not per second as in JULES:
kaps_rothC = [ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10]

tauout   = 2.0
tauin    = 5.0
botdepth = 3.0
depthuni = 1.0

inf=nc.Dataset(infile,'r')
CS     = inf.variables['cs'][:,:,:,0,:].astype('float64')
RESP_S = inf.variables['resp_s'][:,:,:,0,:].astype('float64')
LIT_C  = inf.variables['lit_c'][:,:,0,:].astype('float64')
FRAC   = inf.variables['frac'][:,:13,0,:].astype('float64')

ntime=CS.shape[0]
npools=CS.shape[1]
nland=CS.shape[3]

if L_update_dumpfile:
    Dinf = nc.Dataset(Dinfile,'a')
else:
    Dinf = nc.Dataset(Dinfile,'r')
clay = np.array([Dinf.variables['clay'][:].astype('float64') for i in range(ntime) ])   #0.23
beta = 1.0 / (4.0895 + 2.672*np.exp(-0.0786 * 100.0 * clay))
T_soil=Dinf.variables['t_soil'][:].astype('float64')

# Array to store equilibrium soil carbon
CS_equi = np.zeros_like(CS)

# Working array's for storing Cs-eq by pool
CS_dpm = np.zeros_like(CS[:,0,:])
CS_rpm = np.zeros_like(CS[:,1,:])
CS_bio = np.zeros_like(CS[:,2,:])
CS_hum = np.zeros_like(CS[:,3,:])

Carb   = np.zeros_like(CS[:,3,:])  # Working variable for BIO and HUM calcualtion

CS_dpm_old = np.copy(CS[:,0,:])

resps_dpm  = RESP_S[:,0,:,:]/kaps_rothC[0]
tmfact_dpm = resps_dpm/CS_dpm_old

T_soil_bot = T_soil[-1,:]
mix = np.zeros([ntime,nland])+0.0001
# Extra line to set points where bottom layer raises above 273.15 in the year:
mix[:,T_soil_bot<273.15]=0.0005

LIT_C *= FRAC
#LIT_C[LIT_C<1e-11]=1e-11
LIT_TOT = np.sum(LIT_C,axis=1)
dpm = np.zeros([ntime,nland])
for ipft in range(npft):
    dpm += LIT_C[:,ipft,:]*(dpm_rpm_ratio[ipft]/(1+dpm_rpm_ratio[ipft]))
fracdpm = dpm/np.sum(LIT_C,axis=1)
fracdpm[np.isnan(fracdpm)]=0.5

# arrays for DPM and RPM terms
a=np.zeros([ntime,nz,nland])
b=np.zeros([ntime,nz,nland])

# Arrays for BIO and HUM terms
aa=np.zeros([ntime,nz,nland])
bb=np.zeros([ntime,nz,nland])
d=np.zeros([ntime,nz,nland])
dd=np.zeros([ntime,nz,nland])

#######################################################################################
#  DPM 
############################
# Set bottom layer values of a and b:
a[:,nz-2,:]=1/(1+(  0.5*dz[-1]*(dz[-1]+dz[-2])*kaps[0]*tmfact_dpm[:,-1,:]
              / Dif(mix,z[-1],depthuni,botdepth) )   
            )

b[:,nz-2,:]= a[:,nz-2,:]*0.5*dz[-1]*(dz[-1]+dz[-2])*fracdpm      \
            * Rin(tauin,z[nz-1]+0.5*dz[nz-1],LIT_TOT)    \
            / Dif(mix,z[nz-1],depthuni,botdepth)

# * Upwards sweep for decomposable plant material (DPM) */
for i in range(nz-2,0,-1):
    
    a[:,i-1,:]= 1/(1 + (  dz[i]+dz[i-1])  
                        * (   Dif(mix,z[i+1],depthuni,botdepth)*(1-a[:,i,:])/(dz[i]+dz[i+1])    
                           + 0.5*dz[i]*kaps[0]*tmfact_dpm[:,i,:] ) 
                        / Dif(mix,z[i],depthuni,botdepth) )
    
    b[:,i-1,:] =   (dz[i]+dz[i-1]) * a[:,i-1,:]  \
                 * (  0.5 * dz[i] * fracdpm * Rin(tauin,z[i]+0.5*dz[i],LIT_TOT)
                     + Dif(mix,z[i+1],depthuni,botdepth)*b[:,i,:]/(dz[i]+dz[i+1]) ) \
                 / Dif(mix,z[i],depthuni,botdepth)
            
CS_dpm[:,0,:] =  (   b[:,0,:] + 0.5 * dz[0] * (dz[0]+dz[1]) * fracdpm * Rin(tauin,0.5*dz[0],LIT_TOT)
                   / Dif(mix,z[1],depthuni,botdepth)   ) \
                /(  1 - a[:,0,:] + 0.5*dz[0]*(dz[0]+dz[1])*kaps[0]*tmfact_dpm[:,0,:]
                / Dif(mix,z[1],depthuni,botdepth))

# * Downwards sweep for DPM *
for i in range(1,nz):
    CS_dpm[:,i,:] = CS_dpm[:,i-1,:]*a[:,i-1,:] + b[:,i-1,:]
CS_dpm[np.isnan(CS_dpm)]=1e-11
CS_dpm[np.isinf(CS_dpm)]=1e-11

#######################################################################################
#  RPM 
############################
# Set bottom layer values of a and b:
a[:,nz-2,:]=1/(1+(  0.5*dz[-1]*(dz[-1]+dz[-2])*kaps[1] * tmfact_dpm[:,nz-1,:] 
              / Dif(mix,z[-1],depthuni,botdepth) )   
            )

b[:,nz-2,:]= a[:,nz-2,:]*0.5*dz[-1]*(dz[-1]+dz[-2])*(1-fracdpm)  \
            * Rin(tauin,z[nz-1]+0.5*dz[nz-1],LIT_TOT)    \
            / Dif(mix,z[nz-1],depthuni,botdepth)

# * Upwards sweep for decomposable plant material (RPM) */
for i in range(nz-2,0,-1):
    a[:,i-1,:]= 1/(1 + (  dz[i]+dz[i-1])  
                        * (   Dif(mix,z[i+1],depthuni,botdepth)*(1-a[:,i,:])/(dz[i]+dz[i+1])    
                           + 0.5*dz[i]*kaps[1]*tmfact_dpm[:,i,:] ) 
                        / Dif(mix,z[i],depthuni,botdepth) )
    
    b[:,i-1,:] =   (dz[i]+dz[i-1]) * a[:,i-1,:]   \
                 * (  0.5 * dz[i] * (1-fracdpm) * Rin(tauin,z[i]+0.5*dz[i],LIT_TOT)
                     + Dif(mix,z[i+1],depthuni,botdepth)*b[:,i,:]/(dz[i]+dz[i+1]) )  \
                 / Dif(mix,z[i],depthuni,botdepth)
            
CS_rpm[:,0,:] =  (   b[:,0,:] + 0.5 * dz[0] * (dz[0]+dz[1]) * (1-fracdpm) * Rin(tauin,0.5*dz[0],LIT_TOT)
                   / Dif(mix,z[1],depthuni,botdepth)   )  \
                /(  1 - a[:,0,:] + 0.5*dz[0]*(dz[0]+dz[1])*kaps[1]*tmfact_dpm[:,0,:]
               / Dif(mix,z[1],depthuni,botdepth))

# * Downwards sweep for RPM *
for i in range(1,nz):
    CS_rpm[:,i,:] = CS_rpm[:,i-1,:]*a[:,i-1,:] + b[:,i-1,:]
CS_rpm[np.isnan(CS_rpm)]=1e-11
CS_rpm[np.isinf(CS_rpm)]=1e-11


#######################################################################################
#  BIO  and HUM
############################
# Now we define d, a', b' and d' as well as a and b So that we can solve for BIO and HUM (which are coupled)
# separately define kC + k'C' 
for i in range(nz):
    Carb[:,i,:] =  CS_dpm[:,i,:]*kaps[0]*tmfact_dpm[:,i,:] \
                 + CS_rpm[:,i,:]*kaps[1]*tmfact_dpm[:,i,:] 

A = 1 + ( kaps[2]*tmfact_dpm[:,-1,:] \
        * 0.5*dz[-1]*(dz[-1]+dz[-2])*(1-(0.46*beta))/Dif(mix,z[-1],depthuni,botdepth) )
B = 0.5*dz[-1]*(dz[-1]+dz[-2])*0.46*beta*Carb[:,-1,:]/Dif(mix,z[-1],depthuni,botdepth)
C = 0.5*dz[-1]*(dz[-1]+dz[-2])*0.46*beta*kaps[3]*tmfact_dpm[:,-1,:] \
        / Dif(mix,z[-1],depthuni,botdepth)

D = 1 + ( kaps[3]*tmfact_dpm[:,-1,:] \
        * 0.5*dz[-1]*(dz[-1]+dz[-2])*(1-(0.54*beta))/Dif(mix,z[-1],depthuni,botdepth) )
E = 0.5*dz[-1]*(dz[-1]+dz[-2])*0.54*beta*Carb[:,-1,:]/Dif(mix,z[-1],depthuni,botdepth)
F = 0.5*dz[-1]*(dz[-1]+dz[-2])*0.54*beta*kaps[2]*tmfact_dpm[:,-1,:]\
        / Dif(mix,z[-1],depthuni,botdepth)

# /* Upwards sweep for BIO and HUM */
for i in range(nz-2,-1,-1):
    
    a[:,i,:] = 1/(A - (C*F/D))
    b[:,i,:] = a[:,i,:]*(B + (C*E/D))
    d[:,i,:] = a[:,i,:]*C/D
    aa[:,i,:] = 1/(D - (C*F/A))
    bb[:,i,:] = aa[:,i,:]*(E + (F*B/A))
    dd[:,i,:] = aa[:,i,:]*F/A
    if (i>0):
        A = 1 + (  0.5*(dz[i]+dz[i-1]) 
                 * ( (kaps[2]*tmfact_dpm[:,i,:]*dz[i]*(1-0.46*beta))
                    +(Dif(mix,z[i+1],depthuni,botdepth)*(1-a[:,i,:])*2/(dz[i]+dz[i+1])) )
                 / Dif(mix,z[i],depthuni,botdepth)  )
        B = (    0.5 * (dz[i]+dz[i-1])
              *  ( (dz[i]*0.46*beta*Carb[:,i,:]) + (b[:,i,:]*Dif(mix,z[i+1],depthuni,botdepth)*2/(dz[i]+dz[i+1])) ) 
              /  Dif(mix,z[i],depthuni,botdepth)  )
        C = (    0.5*(dz[i]+dz[i-1])
              *  (  (dz[i]*0.46*beta*kaps[3]*tmfact_dpm[:,i,:]) 
                  + (d[:,i,:]*Dif(mix,z[i+1],depthuni,botdepth)*2/(dz[i]+dz[i+1])) )
              / Dif(mix,z[i],depthuni,botdepth) )
        D = 1 + (  0.5*(dz[i]+dz[i-1]) 
                 * ( (kaps[3]*tmfact_dpm[:,i,:] * dz[i] * (1-0.54*beta))
                    +(Dif(mix,z[i+1],depthuni,botdepth)*(1-aa[:,i,:])*2/(dz[i]+dz[i+1])) )
                 /Dif(mix,z[i],depthuni,botdepth) )
        E = (    0.5*(dz[i]+dz[i-1])
              *  ( (dz[i]*0.54*beta*Carb[:,i,:]) + (bb[:,i,:]*Dif(mix,z[i+1],depthuni,botdepth)*2/(dz[i]+dz[i+1])) )
              / Dif(mix,z[i],depthuni,botdepth) )
        F = (    0.5*(dz[i]+dz[i-1])
              *  ( (dz[i]*0.54*beta*kaps[2]*tmfact_dpm[:,i,:])
                  +(dd[:,i,:]*Dif(mix,z[i+1],depthuni,botdepth)*2/(dz[i]+dz[i+1])) )
              / Dif(mix,z[i],depthuni,botdepth)  )

a[np.isnan(a)] = 1e-11
b[np.isnan(b)] = 1e-11
d[np.isnan(d)] = 1e-11
aa[np.isnan(aa)] = 1e-11
bb[np.isnan(bb)] = 1e-11
dd[np.isnan(dd)] = 1e-11

#* Now we have equations of the form: Cbio[i+1] = Cbio[i]*a[i] + b[i] + d[i]Chum[i] and Chum[i+1] = Chum[i]*aa[i] + bb[i] + dd[i]Cbio[i]
#           Now we just need to solve for Cbio[0] and Chum[0] and we can get all the rest from the above
A = a[:,0,:] - 1 - 0.5*dz[0]*(dz[1]+dz[0])*kaps[2]*tmfact_dpm[:,0,:]*(1-0.46*beta)/Dif(mix,z[1],depthuni,botdepth)
B = b[:,0,:] + 0.5*dz[0]*(dz[1]+dz[0])*0.46*beta*Carb[:,0,:]/Dif(mix,z[1],depthuni,botdepth)
C = d[:,0,:] + 0.5*dz[0]*(dz[1]+dz[0])*0.46*beta*kaps[3]*tmfact_dpm[:,0,:]/Dif(mix,z[1],depthuni,botdepth)

D = aa[:,0,:] - 1 - 0.5*dz[0]*(dz[1]+dz[0])*kaps[3]*tmfact_dpm[:,0,:]*(1-0.54*beta)/Dif(mix,z[1],depthuni,botdepth)
E = bb[:,0,:] + 0.5*dz[0]*(dz[1]+dz[0])*0.54*beta*Carb[:,0,:]/Dif(mix,z[1],depthuni,botdepth)
F = dd[:,0,:] + 0.5*dz[0]*(dz[1]+dz[0])*0.54*beta*kaps[2]*tmfact_dpm[:,0,:]/Dif(mix,z[1],depthuni,botdepth)


CS_hum[:,0,:] = ((A*E) - (B*F))/((C*F) - (A*D))
CS_bio[:,0,:] = ((B*D) - (C*E))/((C*F) - (A*D))

for i in range(1,nz):

    CS_bio[:,i,:] = CS_bio[:,i-1,:]*a[:,i-1,:] + b[:,i-1,:] + d[:,i-1,:]*CS_hum[:,i-1,:]
    CS_hum[:,i,:] = CS_hum[:,i-1,:]*aa[:,i-1,:] + bb[:,i-1,:] + dd[:,i-1,:]*CS_bio[:,i-1,:]

# Debugging: To print output in same format as Sarah Chadburn for comparisons
#print('Writing ascii file to: '+outfile.replace('.nc','.dat'))
#outf=open(outfile.replace('.nc','.dat'),'w')
#for iland in range(nland):
#    outf.write(nz*'%10.5f ' % tuple([ CS_dpm[0,iz,iland] for iz in range(nz)])+'\n')
#    outf.write(nz*'%10.5f ' % tuple([ CS_rpm[0,iz,iland] for iz in range(nz)])+'\n')
#    outf.write(nz*'%10.5f ' % tuple([ CS_bio[0,iz,iland] for iz in range(nz)])+'\n')
#    outf.write(nz*'%10.5f ' % tuple([ CS_hum[0,iz,iland] for iz in range(nz)])+'\n')
#outf.close()

# Convert to kg m^-2
for iz in range(nz):
    CS_dpm[:,iz,:] *= dz[iz]
    CS_rpm[:,iz,:] *= dz[iz]
    CS_bio[:,iz,:] *= dz[iz]
    CS_hum[:,iz,:] *= dz[iz]



CS_equi[:,0,:]=CS_dpm
CS_equi[:,1,:]=CS_rpm
CS_equi[:,2,:]=CS_bio
CS_equi[:,3,:]=CS_hum
CS_equi[np.isnan(CS_equi)]=1e-6

print('Producing output file: '+outfile)
outf=nc.Dataset(outfile,'w')

for var in outvars:
    if var in inf.variables.keys():
        invar=inf.variables[var]
        tempfile=inf
    elif var in Dinf.variables.keys():
        invar=Dinf.variables[var]
        tempfile=Dinf
    else:
        print('No template variable in input file for: '+var)
        break
    
    for dim in invar.dimensions:
        if dim not in outf.dimensions:
            outf.createDimension(dim,len(tempfile.dimensions[dim]))

    outvar=outf.createVariable(var,'float32',invar.dimensions)
    for att in invar.ncattrs():
        outvar.setncattr( str(att),invar.getncattr(str(att)) )

    if var == 'cs':
        outvar[:]=CS_equi 
    else:
        outvar[:]=invar[:]

outf.setncattr("Title","Layered Soil Carbon Equilibrium calculated using Sarah Chadburn's analytic soultion")
outf.setncattr("Note","Created using Eddy Comyn-Platt's (edwcom@ceh.ac.uk) python script")

outf.close()

if L_update_dumpfile:
    print('Updating Cs in dumpfile: '+Dinfile)
    Dinf.variables['cs'][:]=CS_equi
Dinf.close()

inf.close()

print('equilibrium_soilcarbon.py finished normally')

