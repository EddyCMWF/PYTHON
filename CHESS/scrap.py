
SMmin=SMwilt_sthu*1.7
SMopt=(SMwilt_sthu+1)*0.5

plt.subplot(2,4,1)
plt.imshow(SMwilt_sthu[grindex]*grimask,origin='bottom',vmax=1,vmin=0)
plt.colorbar()

plt.subplot(2,4,2)
plt.imshow(SMmin[grindex]*grimask,origin='bottom',vmax=1,vmin=0)
plt.colorbar()

plt.subplot(2,4,3)
plt.imshow(SMopt[grindex]*grimask,origin='bottom',vmax=1,vmin=0)
plt.colorbar()

plt.subplot(2,4,4)
plt.imshow((SMopt-SMmin)[grindex]*grimask,origin='bottom',vmax=1,vmin=-1)
plt.colorbar()

plt.subplot(2,4,5)
plt.imshow(SMsat[grindex]*grimask,origin='bottom',vmax=1,vmin=0)
plt.colorbar()

plt.subplot(2,4,6)
plt.imshow(np.mean(SM_sthu_climat,axis=0)[grindex]*grimask,origin='bottom',vmax=1,vmin=0)
plt.colorbar()

plt.subplot(2,4,7)
plt.imshow(np.mean(Rs_SM_f,axis=0)[grindex]*grimask,origin='bottom')
plt.colorbar()

plt.subplot(2,4,8)
plt.imshow(np.mean(Soil_Resp_Fact_Q10t[0],axis=0)[grindex]*grimask,origin='bottom')
plt.colorbar()

plt.show()

