def parabolic(kappa, delta_temp_ocean_yearly):
    # This is an implicit solver for parabolic ocean PDE, with Dirichlet condition satisfied at the surface 
    # C. Huntingford (11th March 2017)

    # Variables in are kappa=ocean diffusion (W/m/K), yearly top level ocean temperature delta_temp_ocean (K) 
    # When run in CEH, on linux box wllf017 or higher. To get newest version of matrix inversion. 
    import numpy as np
    import sys,os
    from scipy import sparse
    from scipy.sparse.linalg import spsolve 
    from scipy.interpolate import interp1d

    # Set volumetric heat capacity for salt water, and fraction of Earth with ocean
    cp = 4.04E6                                    # (J/K/m3)

    # Now set up the conditions for solution to the parabolic diffusion equation. 
    # Set number of calls of ocean diffusion PDE solver per year.
    n_pde = 20
#    n_pde = 100
    dt = (1.0/float(n_pde))*60.0*60.0*24.0*365.0   # Model timestep (seconds)

    # Set number of vertical layers and ocean depth 
    n_vert = 254                                   # So n_vert+1 grid points
#    n_vert = 1000                                  # So n_vert+1 grid points
    depth = 5000.0                                 # (metres)
    dz = depth/float(n_vert)                       # (metres)

    # Set the output variable
    nyr=delta_temp_ocean_yearly.shape[0]
    print(nyr)
    temp_gradient_out=np.zeros(nyr)                # dT/dz at the ocean surface (K/m) 

    # Prepare for solving PDE
    s=(kappa/cp)*(dt/(dz*dz))
    print(s)
    t_ocean_old=np.zeros(n_vert-1)                 # Does not include top level (prescribed) or bottom level
    t_ocean_new=np.zeros(n_vert-1)                 # Does not include top level (prescribed) or bottom level

    #######################################################
    # The numerical scheme with Direchlet conditions is as below: 
    # Solve for t_ocean_old for 1 <= i < n_vert-1 
    # Solve for t_ocean_new for 1 <= i < n_vert-1 
    # 
    # Then A x t_ocean_new + C = B x t_ocean_old + D 
    # A is (n_vert-1 x n_vert-1) 
    # B is (n_vert-1 x n_vert-1) 

    # A = [2(1+s)  -s  ..     ..    ..   0     0]   C = [-s x delta_temp_ocean(t+dt) ] 
    #     [ -s    2(1+s)  -s  ..    ..   0     0]       [0]
    #     [ ..   ..       ..  ..    ..   ..   ..]       [0]
    #     [ ..   ..       ..  ..   -s 2(1+s)  -s]       [0] 
    #     [ ..   ..       ..  ..    ..-s  2(1+s)]       [0] 
    # 
    # 
    # B = [2(1-s)   s  ..     ..    ..   0     0]   D = [s x delta_temp_ocean(t) ] 
    #     [  s    2(1-s)   s  ..    ..   0     0]       [0]
    #     [ ..   ..       ..  ..    ..   ..   ..]       [0]
    #     [ ..   ..       ..  ..    s 2(1-s)   s]       [0] 
    #     [ ..   ..       ..  ..    .. s  2(1-s)]       [0] 
    # 
    #######################################################

    A=np.zeros([n_vert-1, n_vert-1])
    C=np.zeros(n_vert-1)
    B=np.zeros([n_vert-1, n_vert-1])
    D=np.zeros(n_vert-1)

    # A and B are invariant in time, so can be set now
    for i_pos in range(0,n_vert-1):
        A[i_pos, i_pos]=2.0*(1.0+s)
    for i_pos in range(0,n_vert-1-1):
        A[i_pos+1, i_pos] = -s
        A[i_pos, i_pos+1] = -s

    for i_pos in range(0,n_vert-1):
        B[i_pos, i_pos]=2.0*(1.0-s)
    for i_pos in range(0,n_vert-1-1):
        B[i_pos+1, i_pos] = s
        B[i_pos, i_pos+1] = s

    A_mat=np.mat(A)
    B_mat=np.mat(B)

    print(A[:4,:4])
    print(B[:4,:4])
    # Conservation check for end of run - add up energy in. 
    q_energy_derived=0.0 ; q_total=0.0

     # Now start loop over the different years. 
    for j in range(1,nyr):
        for k in range(0,n_pde):
            #print(k)
            delta_temp_ocean_yearly_sub_yearly_t =  delta_temp_ocean_yearly[j-1] + \
                (np.float(k)/np.float(n_pde))*(delta_temp_ocean_yearly[j] - delta_temp_ocean_yearly[j-1])
            delta_temp_ocean_yearly_sub_yearly_t_dt =  delta_temp_ocean_yearly[j-1] + \
                (np.float(k+1)/np.float(n_pde))*(delta_temp_ocean_yearly[j] - delta_temp_ocean_yearly[j-1])

            C[0] = -s*delta_temp_ocean_yearly_sub_yearly_t_dt
            D[0] = s*delta_temp_ocean_yearly_sub_yearly_t
            #print(k, C[0],D[0])

            C_mat=np.mat(C).transpose()
            D_mat=np.mat(D).transpose()

            # Now combine all the right-hand side terms
            t_ocean_old=t_ocean_new
            t_ocean_old_mat=np.mat(t_ocean_old).transpose()
            #print('t_o:',t_ocean_old[:4])
            b_rhs=np.dot(B,t_ocean_old_mat)
            #print('b.t_o:',(np.ravel(b_rhs[:4])))
            b_rhs = b_rhs + D_mat - C_mat
            b_rhs=np.ravel(b_rhs)                             # Flattens to 1-D array
            #print('b_rhs:', b_rhs[:4])

# Now perform the calculation to update the oceanic temperatures 
            t_ocean_new=spsolve(A,b_rhs)

# Calculate the heat flux at each timestep.
            q_energy = -kappa*((t_ocean_new[0] - delta_temp_ocean_yearly_sub_yearly_t_dt)/dz)
            q_total=q_total+dt*q_energy

            # Calculate yearly mean flux in to the ocean
            temp_gradient_out[j] = temp_gradient_out[j] + q_energy*(1.0/np.float(n_pde))
        print(delta_temp_ocean_yearly[j],t_ocean_new[0:4])
        print(temp_gradient_out,q_total)
# Check conservation of energy at end of each year. For a column of size one metre^2. 
        q_energy_derived=0.0
        for l in range(1,n_vert-1):                # Points where heat calculated
            q_energy_derived = q_energy_derived + cp*0.5*(t_ocean_new[l]+t_ocean_new[l-1])*dz 
        q_energy_derived = q_energy_derived + cp*0.5*(t_ocean_new[1]+delta_temp_ocean_yearly_sub_yearly_t_dt)*dz    # Top layer
        print( q_energy_derived)
        q_energy_derived = q_energy_derived + cp*0.5*(t_ocean_new[n_vert-2])*dz                                     # Bottom layer

#        print 'Heat conservation check (%) = ', j+1850, 100.0*(q_energy_derived/q_total)
    print( 'Heat conservation check final year (%) = ', j+1850, 100.0*(q_energy_derived/q_total) )

    return temp_gradient_out
