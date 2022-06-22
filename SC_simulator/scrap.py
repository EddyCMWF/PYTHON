
plots_DIR = '/users/eow/edwcom/SC_simulator/output/plots/single_timestep/'


n96_landindexfile='/users/eow/edwcom/CRUNCEP/n96/ancil/jules_land_index.nc'

# Open and read the land index file
inf=nc.Dataset(n96_landindexfile,'r')
index=inf.variables['index_2D'][:]
lats_2D=inf.variables['lats_2D'][:]
lons_2D=inf.variables['lons_2D'][:]
inf.close()


Leaf_lit_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Leaf_lit_origgrid[:,soil_points]=Leaf_LF
Leaf_lit_origgrid=np.ma.masked_equal(Leaf_lit_origgrid,fill_value)

Root_lit_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Root_lit_origgrid[:,soil_points]=Root_LF
Root_lit_origgrid=np.ma.masked_equal(Root_lit_origgrid,fill_value)

Stem_lit_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Stem_lit_origgrid[:,soil_points]=Stem_LF
Stem_lit_origgrid=np.ma.masked_equal(Stem_lit_origgrid,fill_value)

Leaf_lit_2D=Leaf_lit_origgrid[:,index]
Root_lit_2D=Root_lit_origgrid[:,index]
Stem_lit_2D=Stem_lit_origgrid[:,index]

for i in range(5):
    Leaf_lit_2D[i,:,:][np.where(index.mask==True)]=fill_value
    Root_lit_2D[i,:,:][np.where(index.mask==True)]=fill_value
    Stem_lit_2D[i,:,:][np.where(index.mask==True)]=fill_value


Leaf_lit_2D=np.ma.masked_equal(Leaf_lit_2D,-999.)
Root_lit_2D=np.ma.masked_equal(Root_lit_2D,-999.)
Stem_lit_2D=np.ma.masked_equal(Stem_lit_2D,-999.)



plt.subplot(1,3,1)
plt.imshow(Leaf_lit_2D[0,:,:],origin='bottom')
plt.colorbar()
plt.title('Leaf')
plt.subplot(1,3,2)
plt.imshow(Root_lit_2D[0,:,:],origin='bottom')
plt.colorbar()
plt.title('Root')
plt.subplot(1,3,3)
plt.imshow(Stem_lit_2D[0,:,:],origin='bottom')
plt.colorbar()
plt.title('Stem')
plt.show()


Lit_C_ECP_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Lit_C_ECP_origgrid[:,soil_points]=Lit_C_ECP
Lit_C_ECP_origgrid=np.ma.masked_equal(Lit_C_ECP_origgrid,fill_value)

Lit_C_ECP_2D=Lit_C_ECP_origgrid[:,index]
for i in range(5):
    Lit_C_ECP_2D[i,:,:][np.where(index.mask==True)]=fill_value

Lit_C_ECP_2D=np.ma.masked_equal(Lit_C_ECP_2D,-999.)


LocLit_C_ECP_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
LocLit_C_ECP_origgrid[:,soil_points]=LocLit_C_ECP
LocLit_C_ECP_origgrid=np.ma.masked_equal(LocLit_C_ECP_origgrid,fill_value)

LocLit_C_ECP_2D=LocLit_C_ECP_origgrid[:,index]
for i in range(5):
    LocLit_C_ECP_2D[i,:,:][np.where(index.mask==True)]=fill_value

LocLit_C_ECP_2D=np.ma.masked_equal(LocLit_C_ECP_2D,-999.)


Tleaf_J_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Tleaf_J_origgrid[:,soil_points]=Tleaf_JULES
Tleaf_J_origgrid=np.ma.masked_equal(Tleaf_J_origgrid,fill_value)

Tleaf_J_2D=Tleaf_J_origgrid[:,index]
for i in range(5):
    Tleaf_J_2D[i,:,:][np.where(index.mask==True)]=fill_value

Tleaf_J_2D=np.ma.masked_equal(Tleaf_J_2D,-999.)




Cveg_origgrid=np.zeros([5,lats_orig.shape[0]])+fill_value
Cveg_origgrid[:,soil_points]=Cveg_ECP
Cveg_origgrid=np.ma.masked_equal(Cveg_origgrid,fill_value)

Cveg_2D=Cveg_origgrid[:,index]
for i in range(5):
    Cveg_2D[i,:,:][np.where(index.mask==True)]=fill_value

Cveg_2D=np.ma.masked_equal(Cveg_2D,-999.)


data=Cveg_2D[0,:,:]
plt.subplot(1,2,1)
plt.imshow(data,origin='bottom')
plt.colorbar()

data=LocLit_C_ECP_2D[0,:,:]
plt.subplot(1,2,2)
plt.imshow(data,origin='bottom')
plt.colorbar()
plt.show()



C_4pools_Q10t_2D= C_4pools_Q10t[:,index]
for i in range(4):
    C_4pools_Q10t_2D[i,:,:][np.where(index.mask==True)]=fill_value

C_4pools_Q10t_2D=np.ma.masked_equal(C_4pools_Q10t_2D,-999.)


Lit_SCpools_2D= Lit_SCpools_origgrid[:,index]
for i in range(4):
    Lit_SCpools_2D[i,:,:][np.where(index.mask==True)]=fill_value

Lit_SCpools_2D=np.ma.masked_equal(Lit_SCpools_2D,-999.)


Soil_Resp_Fact_2D= Soil_Resp_Fact_Q10t_origgrid[:,index]
for i in range(4):
    Soil_Resp_Fact_2D[i,:,:][np.where(index.mask==True)]=fill_value

Soil_Resp_Fact_2D=np.ma.masked_equal(Soil_Resp_Fact_2D,-999.)



plt.subplot(2,3,1)
plt.imshow(C4pools_Q10t_2D[0,:,:],origin='bottom')
plt.title='Soil Carbon'
plt.colorbar(fraction=0.05)

plt.subplot(2,3,2)
plt.imshow(Lit_SCpools_2D[0,:,:],origin='bottom')
plt.title='Litter'
plt.colorbar(fraction=0.05)

plt.subplot(2,3,3)
plt.imshow(Soil_Resp_Fact_2D[0,:,:]*(3600.*24.*360.),origin='bottom')
plt.title='Soil Resp. Factor'
plt.colorbar(fraction=0.05)


plt.subplot(2,3,4)
plt.imshow(C4pools_Q10t_2D[1,:,:],origin='bottom')
plt.colorbar(fraction=0.05)


plt.subplot(2,3,5)
plt.imshow(Lit_SCpools_2D[1,:,:],origin='bottom')
plt.colorbar(fraction=0.05)


plt.subplot(2,3,6)
plt.imshow(Soil_Resp_Fact_2D[1,:,:]*(3600.*24.*360.),origin='bottom')
plt.colorbar(fraction=0.05)

plt.show()





colors=['y','g','orange','saddlebrown']
labels=['DPM','RPM','BIO','HUM']

fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothCt[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothCt[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()



fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_ECP[i,:],lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10t[i,:]*1e6,lats_JULES,ls='',marker='.',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('ukgC y-1')
#ax=fig.add_subplot(1,3,3)
#ax.plot(C_4pools_RothCt[1,:],lats_JULES,ls='',marker='.')
#ax.set_title('Carbon Pool')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10t[i,:],lats_JULES,ls='',marker='.',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()



lat_binsize=5.
minlat=-60.
maxlat=90.
lat_bins=np.arange(minlat,maxlat,lat_binsize)

Lit_SCpools_Zonal=np.zeros([2,len(lat_bins)])
Soil_Resp_Fact_Q10_Zonal=np.zeros([4,len(lat_bins)])
Soil_Resp_Fact_RothC_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_Q10_Zonal=np.zeros([4,len(lat_bins)])
C_4pools_RothC_Zonal=np.zeros([4,len(lat_bins)])


for bin in range(len(lat_bins)):
    tempdex=np.where( (lats_JULES>=lat_bins[bin]        ) & \
                    (lats_JULES< lat_bins[bin]+lat_binsize) )
    #
    Lit_SCpools_Zonal[:,bin]=np.mean(Lit_SCpools_ECP[:,tempdex],axis=1)
    Soil_Resp_Fact_Q10_Zonal[:,bin]=np.mean(Soil_Resp_Fact_Q10t[:,tempdex],axis=1)
    Soil_Resp_Fact_RothC_Zonal[:,bin]=np.mean(Soil_Resp_Fact_RothCt[:,tempdex],axis=1)
    C_4pools_Q10_Zonal[:,bin]=np.mean(C_4pools_Q10t[:,tempdex],axis=1)
    C_4pools_RothC_Zonal[:,bin]=np.mean(C_4pools_RothCt[:,tempdex],axis=1)




fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_Q10_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_Q10_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
#ax.set_xlim(0,1e10)

fig.suptitle('Q10',fontsize=18)
plt.show()



fig=plt.figure()
ax=fig.add_subplot(1,3,1)
for i in [1,0]:
    ax.plot(Lit_SCpools_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Litter')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')

ax=fig.add_subplot(1,3,2)
for i in  [0,3,2,1]:
    ax.plot(Soil_Resp_Fact_RothC_Zonal[i,:],lat_bins,ls='-',color=colors[i])

ax.set_title('Soil Resp.')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC y-1')
ax=fig.add_subplot(1,3,3)
for i in [3,1,2,0]:
    ax.plot(C_4pools_RothC_Zonal[i,:],lat_bins,ls='-',label=labels[i],color=colors[i])

ax.set_title('Carbon Pool')
ax.set_ylabel('latitude')
ax.set_xlabel('kgC')
ax.legend(loc='lower right')
ax.set_xlim(0,1e10)

fig.suptitle('RothC',fontsize=18)
plt.show()

