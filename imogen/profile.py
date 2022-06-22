import numpy as np
import os, sys

def profile(beta, dt_now, dt_limit, mu_zero, mu_one, yr_now=2015, yr_beta=1995, smooth_years=0, smooth_iters=1):
    # Builds temperature profile based on parameters given. 1850-2015 + 500 years 
    # C. Huntingford (13th March 2017)

    delta_temp_global=np.zeros(2015-1850+1+500)

    # Read in file of HadCRUT. First get number of lines. These are used up to year 1995.
    file_nam='/users/global/chg/CLUES/profile_paper/HadCRUT/HadCRUT.4.5.0.0.annual_ns_avg_smooth.rtf'
    f=open(file_nam, 'r')
    n_lines=0
    for line in f:
        n_lines = n_lines+1 
    f.close()

    yr=np.zeros(n_lines, dtype='int')
    dtemp_median=np.zeros(n_lines)

    # Read in the 1st column (years) and the 2nd column (median).
    f=open(file_nam, 'r')
    for i_yr in range(0,n_lines):
        line_in = f.readline()
        vals_in = line_in.split()
        yr[i_yr]=vals_in[0]
        dtemp_median[i_yr]=vals_in[1]
    
    if smooth_years!=0:
        smooth_range=int(smooth_years/2.)
        for i_smooth in range(smooth_iters):
            for i_yr in range(n_lines):
                if i_yr < smooth_range:
                    dtemp_median[i_yr]=np.mean(dtemp_median[:i_yr+smooth_range])
                elif i_yr > n_lines-smooth_range:
                    dtemp_median[i_yr]=np.mean(dtemp_median[i_yr-smooth_range:])
                else:
                    dtemp_median[i_yr]=np.mean(dtemp_median[i_yr-smooth_range:i_yr+smooth_range])

    # Normalise with the first 50 years
    start_51 = np.mean(dtemp_median[0:51])
    dtemp_median = dtemp_median - start_51

    # Extra the 1995 value from the median numbers. To ensure smooth transition, normalise so line up OK year 1995.
    nbetayrs=yr_now-yr_beta
    beta_index = np.where(yr==yr_beta)[0]
    now_index = beta_index+nbetayrs

    median_beta = dtemp_median[beta_index]
    curve_beta = dt_now - float(nbetayrs)*beta

    dtemp_median = dtemp_median*(curve_beta/median_beta)
    #print(beta_index, yr[beta_index])
    # Now fill in up to year 1995
    for i in range(0,beta_index,1):    # 10):
        delta_temp_global[i]=dtemp_median[i]
        #diff=dtemp_median[i+10]-dtemp_median[i]
        #for j in range(10):
        #    delta_temp_global[i+j] =  delta_temp_global[i]  \
        #                            + (diff*float(j)/10.) 

    # Now for year 1995-yr_now (changed to linear)
    for i in range(0,nbetayrs+1):
        delta_temp_global[beta_index+i] = dt_now - float(nbetayrs)*beta + np.float(i)*beta
        #delta_temp_global[145+i] = delta_temp_global[144+i] \
        #                        + ( (dt_now-delta_temp_global[145])/20.)
    
    #print(now_index, yr[now_index])
    # Now for year 2016-2515 inc 
    gamma = beta - mu_zero*(dt_limit-dt_now)
    for i in range(0,500):
        delta_temp_global[i+now_index] = \
             dt_now + gamma*np.float(i) - (1.0 - np.exp(-(mu_zero+mu_one*np.float(i))*np.float(i))) \
           * (gamma*np.float(i) - (dt_limit - dt_now))
          #   dt_now + gamma*np.float(i+1) - (1.0 - np.exp(-(mu_zero+mu_one*np.float(i+1))*np.float(i+1))) \
          # * (gamma*np.float(i+1) - (dt_limit - dt_now))


    return delta_temp_global


def prescribed_histo_ECP(dt_limit, mu_zero, mu_one,
                         delta_temp_global_histo,
                         dt_now=0.89,
                         start_year=1850,end_year=2100,
                         trans_year=1995,now_year=2015):
    # Builds temperature profile based on parameters given. 1850-2015 + 500 years 
    # C. Huntingford (13th March 2017)
    print(start_year,end_year)
    years=np.arange(start_year,end_year+1)
    delta_temp_global=np.zeros(end_year-start_year+1)

    trans_index = np.where(years==trans_year)[0][0]
    now_index = np.where(years==now_year)[0][0]
    trans_range=(now_index-trans_index)
    end_index=len(delta_temp_global)-1
    curve_range=(end_index-now_index)


    # Now fill in up to trans_year
    delta_temp_global[:trans_index]=delta_temp_global_histo[:trans_index]

    # To ensure smooth transition, linearise between 1995 and 2015
    #dt_now   = delta_temp_global[now_index]
    dt_trans = delta_temp_global[trans_index]
    temp_range = dt_now-dt_trans
    beta   = temp_range/trans_range
    
    # Now for year 1995-2015 (changed to linear)
    for i in range(trans_range+1):
        index=i+trans_index
        t_increment = (float(i)/float(trans_range))
        delta_temp_global[index] = dt_now - float(trans_range)*beta + np.float(i)*beta
        
    # Now for year 2016-2515 inc 
    gamma = beta - mu_zero*(dt_limit-dt_now)
    for i in range(curve_range):
        index=i+now_index+1
        delta_temp_global[index] = \
            dt_now + gamma*np.float(i+1) \
          -   (1.0 - np.exp(-(mu_zero+mu_one*np.float(i+1))*np.float(i+1))) \
            * (gamma*np.float(i+1) - (dt_limit - dt_now))

    return beta, delta_temp_global


def prescribed_histo_ECP_2(dt_limit, mu_zero, mu_one,
                           delta_temp_global_histo,
                           now_year=2015,
                           start_year=1850,end_year=2100,
                           beta_range=5):

    # Builds temperature profile based on parameters given. 1850-2015 + 500 years 
    # C. Huntingford (13th March 2017)
    print(start_year,end_year)
    years=np.arange(start_year,end_year+1)
    delta_temp_global=np.zeros(end_year-start_year+1)

    now_index = np.where(years==now_year)[0][0]
    end_index=len(delta_temp_global)-1
    curve_range=(end_index-now_index)

    # Now fill in up to trans_year
    delta_temp_global[:now_index+1]=delta_temp_global_histo[:now_index+1]

    # Now calculate beta over the beta range
    dt_now=delta_temp_global[now_index]
    temp_range = dt_now-delta_temp_global[now_index-beta_range]
    beta   = temp_range/beta_range
    
    # Now for year 2016-2515 inc 
    gamma = beta - mu_zero*(dt_limit-dt_now)
    for i in range(curve_range):
        index=i+now_index+1
        delta_temp_global[index] = \
            dt_now + gamma*np.float(i+1) \
          -   (1.0 - np.exp(-(mu_zero+mu_one*np.float(i+1))*np.float(i+1))) \
            * (gamma*np.float(i+1) - (dt_limit - dt_now))

    return beta, delta_temp_global


def prescribed_histo_BC(dt_limit, mu_zero, mu_one,
                       delta_temp_global_histo,
                       dt_now=0.89,now_year=2015,
                       start_year=1850,end_year=2100,
                       beta_range=5):
    # Builds temperature profile based on parameters given. 1850-2015 + 500 years 
    # C. Huntingford (13th March 2017)
    years=np.arange(start_year,end_year+1)
    delta_temp_global=np.zeros(end_year-start_year+1)

    now_index = np.where(years==now_year)[0][0]
    end_index=len(delta_temp_global)-1
    curve_range=(end_index-now_index)

    #calculate the difference between histo and required dt_now
    delta_dt_now = dt_now-delta_temp_global_histo[now_index]

    # Now fill in up to trans_year
    delta_temp_global[:now_index+1]=delta_temp_global_histo[:now_index+1]+delta_dt_now

    # Now calculate beta over the beta range
    temp_range = dt_now-delta_temp_global[now_index-beta_range]
    beta   = temp_range/beta_range
    
    # Now for year 2016-2515 inc 
    gamma = beta - mu_zero*(dt_limit-dt_now)
    for i in range(curve_range):
        index=i+now_index+1
        delta_temp_global[index] = \
            dt_now + gamma*np.float(i+1) \
          -   (1.0 - np.exp(-(mu_zero+mu_one*np.float(i+1))*np.float(i+1))) \
            * (gamma*np.float(i+1) - (dt_limit - dt_now))


    return beta, delta_dt_now, delta_temp_global


def prescribed_histo_BC_2(dt_limit_in, mu_zero, mu_one,
                          delta_temp_global_histo,
                          dt_now_in=0.89,now_year=2015,
                          start_year=1850,end_year=2100,
                          beta_range=5):
    # Builds temperature profile based on parameters given. 1850-2015 + 500 years 
    # C. Huntingford (13th March 2017)
    years=np.arange(start_year,end_year+1)
    delta_temp_global=np.zeros(end_year-start_year+1)

    now_index = np.where(years==now_year)[0][0]
    end_index=len(delta_temp_global)-1
    curve_range=(end_index-now_index)

    # Now fill in up to trans_year
    delta_temp_global[:now_index+1]=delta_temp_global_histo[:now_index+1]
    
    # Calculate dt_limit as the desired limit plus the discrepancy between now temperatures
    delta_dt_now = delta_temp_global[now_index]-dt_now_in
    dt_limit=dt_limit_in+delta_dt_now
    dt_now=dt_now_in+delta_dt_now
    
    # Now calculate beta over the beta range
    temp_range = delta_temp_global[now_index]-delta_temp_global[now_index-beta_range]
    beta   = temp_range/beta_range

    # Now for year 2016-2515 inc 
    gamma = beta - mu_zero*(dt_limit-dt_now)
    for i in range(curve_range):
        index=i+now_index+1
        delta_temp_global[index] = \
            dt_now + gamma*np.float(i+1) \
          -   (1.0 - np.exp(-(mu_zero+mu_one*np.float(i+1))*np.float(i+1))) \
            * (gamma*np.float(i+1) - (dt_limit - dt_now))


    return beta, delta_dt_now, delta_temp_global


