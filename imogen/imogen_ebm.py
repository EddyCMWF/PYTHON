import numpy as np
from scipy.sparse.linalg import spsolve

def forward_ebm(dq_all,kappa,f,lambda_l,lambda_o,nu,
                   cp=4.04e6,n_pde=20,n_vert=254,depth=5000.):
    
    # Number of years is number of dq
    nyr = len(dq_all)
    
    dt = (1.0/float(n_pde))*60.0*60.0*24.0*365.0   # Model timestep (seconds)
    
    # thickness of ocean layers
    dz = depth/float(n_vert)
    
    temp_ocean_top = np.zeros(nyr)                 # Array holds oceanic warming (K)
    t_ocean_old = np.zeros(n_vert+1)
    t_ocean_new = np.zeros(n_vert+1)
    
    # Terms for solving the tri-diagonal PDE matrix :
    s=(kappa/cp)*(dt/(dz*dz))
    C=np.zeros([n_vert+1, n_vert+1])
    D=np.zeros([n_vert+1, n_vert+1])
    E=np.zeros(n_vert+1)

    q_energy=0.0
    for j in range(nyr):
        dq = dq_all[j]

        # Following equation 10 of Huntingford and Cox (2000)
        # dividing through by kappa gives a solvable for:
        # dT/dz = factor1 + factor2*T(z=0)
        
        # factor1 = -dQ / (kappa*f)
        factor1 = -dq /(kappa*f)
        
        factor2 = (((1.0-f)*lambda_l*nu)/(kappa*f))  + (lambda_o/kappa)
        
        
        for k in range(n_pde):
            ival= j*n_pde + k
            
            # First Set new to old
            t_ocean_old = np.copy(t_ocean_new)

            # Sort out the top points, bottom points and then all for C  
            # C[0, 0] = (1.0+s) 
            C[0, 0] = s*(1.0+dz*factor2)+1 
            C[0, 1] = -s
            C[n_vert, n_vert-1] = -s
            C[n_vert, n_vert] = (1.0+s)
            for m in range(1, n_vert):
                C[m, m-1] = -s/2.0
                C[m, m+1] = -s/2.0
                C[m, m] = 1.0 + s 

            # Sort out the top points, bottom points and then all for D  
            # D[0, 0] = (1.0-s) 
            D[0, 0] = -s*(1.0+dz*factor2)+1 
            D[0, 1] = s
            D[n_vert, n_vert] = (1.0-s)
            D[n_vert, n_vert-1] = s
            for m in range(1, n_vert):
                D[m, m-1] = s/2.0
                D[m, m+1] = s/2.0
                D[m, m] = 1.0 - s 
            
            # E[0]=((s*dz)/(f*kappa))*(toa_smoothed_norm_sub[ival]+toa_smoothed_norm_sub[ival+1])
            # E[0]=-(2.0*s*dz)*factor1     # Assume slow variation in Q
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
            
            # Update the cumulative heat entering ocean, and also save to top ocean temperature
            #  q_total=q_total+dt*(toa_smoothed_norm_sub[ival]/f)
            q_energy=q_energy+dt*kappa*(-factor1-(t_ocean_new[0]*factor2))

        # Store top layer in output array
        temp_ocean_top[j]=t_ocean_new[0]


    # Check conservation of energy at end of run. Compare against total heat in the final profile.
    q_energy_derived=0.0
    for l in range(1,n_vert+1):
        q_energy_derived=q_energy_derived+cp*0.5*(t_ocean_new[l]+t_ocean_new[l-1])*dz 

    print('Heat conservation check (%) = ',100.0*(q_energy_derived/q_energy))
    
    return temp_ocean_top

def forward_ebm_year(dq, t_ocean_old, kappa,f,lambda_l,lambda_o,nu,
                     cp=4.04e6,n_pde=20,n_vert=254,depth=5000.):
    
    # Number of gcm is years is number of dq vals
    ngcm = len(dq)
    # If ngcm=1, convert values into array and reshape t_ocean_arrays 
    # to give dummy gcm dimension:
    if ngcm==1:
        t_ocean_old = t_ocean_old.reshape([ngcm]+list(t_ocean_old.shape))
        dq=np.array(dq)
        kappa=np.array(kappa)
        f=np.array(f)
        lambda_l=np.array(lambda_l)
        lambda_o=np.array(lambda_o)
        nu=np.array(nu)

    dt = (1.0/float(n_pde))*60.0*60.0*24.0*365.0   # Model timestep (seconds)
     
    # thickness of ocean layers
    dz = depth/float(n_vert)
    
    #temp_ocean_top = np.zeros(nyr)                 # Array holds oceanic warming (K)
    # Now given t_ocean_old = np.zeros(n_vert+1)
    t_ocean_new = np.copy(t_ocean_old)
    print(t_ocean_old.shape) 
    
    # Terms for solving the tri-diagonal PDE matrix :
    s=(kappa/cp)*(dt/(dz*dz))
    C=np.zeros([ngcm, n_vert+1, n_vert+1])
    D=np.zeros([ngcm, n_vert+1, n_vert+1])
    E=np.zeros([ngcm, n_vert+1])

    q_energy=0.0

    # Following equation 10 of Huntingford and Cox (2000)
    # dividing through by kappa gives a solvable for:
    # dT/dz = factor1 + factor2*T(z=0)
        
    # factor1 = -dQ / (kappa*f)
    factor1 = -dq /(kappa*f)
    factor2 = (((1.0-f)*lambda_l*nu)/(kappa*f))  + (lambda_o/kappa)
    print('factor1,factor2=',factor1,factor2)

    for k in range(n_pde):
        ival= k
        
        # First Set new to old
        t_ocean_old = np.copy(t_ocean_new)
        #print(t_ocean_old[0,0,:2])
        # Sort out the top points, bottom points and then all for C  
        # C[0, 0] = (1.0+s) 
        C[:, 0, 0] = s*(1.0+dz*factor2)+1 
        C[:, 0, 1] = -s
        C[:, n_vert, n_vert-1] = -s
        C[:, n_vert, n_vert] = (1.0+s)
        for m in range(1, n_vert):
            C[:, m, m-1] = -s/2.0
            C[:, m, m+1] = -s/2.0
            C[:, m, m] = 1.0 + s 

        # Sort out the top points, bottom points and then all for D  
        # D[0, 0] = (1.0-s) 
        D[:, 0, 0] = -s*(1.0+dz*factor2)+1 
        D[:, 0, 1] = s
        D[:, n_vert, n_vert] = (1.0-s)
        D[:, n_vert, n_vert-1] = s
        for m in range(1, n_vert):
            D[:, m, m-1] = s/2.0
            D[:, m, m+1] = s/2.0
            D[:, m, m] = 1.0 - s 
            
        # E[:,0]=((s*dz)/(f*kappa))*(toa_smoothed_norm_sub[ival]+toa_smoothed_norm_sub[ival+1])
        # E[:,0]=-(2.0*s*dz)*factor1     # Assume slow variation in Q
        E[:,0]=-(kappa/cp)*(dt/dz)*(factor1+factor1)     # Assume slow variation in Q
        
        # Now solve for U_{j+1}. First put bits in to tri-diagonal things for "sparse.spdiags".
        C_sub=np.zeros([ngcm,n_vert+1])
        C_main=np.zeros([ngcm,n_vert+1])
        C_super=np.zeros([ngcm,n_vert+1])
        for m in range(0,n_vert+1): 
            C_main[:,m]=C[:,m,m]
        for m in range(0,n_vert): 
            C_sub[:,m]=C[:,m+1,m]
        C_super[:,m]=C[:,m,m+1]
        
        for i_gcm in range(ngcm):
            # Now calculate the right-hand side (called b_rhs)
            d_mat=np.mat(D[i_gcm,:,:])
            e_mat=np.mat(E[i_gcm,:]).transpose()
            t_ocean_old_mat=np.mat(t_ocean_old[i_gcm,:]).transpose()
            #b_rhs=np.dot(D[i_gcm,:,:],t_ocean_old_mat)+e_mat
            b_rhs=np.dot(d_mat,t_ocean_old_mat)+e_mat
            b_rhs=np.ravel(b_rhs)
            
            # Now perform the calculation to update the oceanic temperatures 
            t_ocean_new[i_gcm,:]=spsolve(C[i_gcm,:],b_rhs)
        print(t_ocean_new[0,0,:4])
        
        ## Update the cumulative heat entering ocean, and also save to top ocean temperature
        ##  q_total=q_total+dt*(toa_smoothed_norm_sub[ival]/f)
        #q_energy=q_energy+dt*kappa*(-factor1-(t_ocean_new[0]*factor2))
        

    # Check conservation of energy at end of run. Compare against total heat in the final profile.
    #q_energy_derived=0.0
    #for l in range(1,n_vert+1):
    #    q_energy_derived=q_energy_derived+cp*0.5*(t_ocean_new[l]+t_ocean_new[l-1])*dz 
    #
    #print('Heat conservation check (%) = ',100.0*(q_energy_derived/q_energy))
    
    return t_ocean_new

