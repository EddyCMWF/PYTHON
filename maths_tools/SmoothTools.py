# ECP - various maths tools
    
import numpy as np

def FFT_smooth(y, x_res,         \
               smooth_factor=5.  ):

    import scipy.fftpack
    N = y.size
    w = scipy.fftpack.rfft(y)
    f = scipy.fftpack.rfftfreq(N, x_res)
    spectrum = w**2
    
    actual_factor= smooth_factor  # 1./(smooth_factor*0.04)
    
    cutoff_idx = spectrum < (spectrum.max()/actual_factor)
    w2 = w.copy()
    w2[cutoff_idx] = 0
    
    y2 = scipy.fftpack.irfft(w2)
    
    return y2


###############################################################
# STDDEV_mask - returns a mask where the standard deviation 
#                of a window of data of length window exceeds 
#                a n amount defined std_limit
#  Only deals with 1D data at present
###############################################################
def STDDEV_mask(data, window=20, std_limit=None):
    
    # number of data to loop over
    nDATA = len(data)

    # Half window for indexing
    half_window = int(window/2.)
    
    # if standard deviation limit not set, set to std of the data array
    if std_limit==None:
        std_limit=np.std(data)

    # define mask as empty boolean array of lenght nDATA
    mask = np.zeros(nDATA,dtype='bool')
        
    # loop over elements from half_window to nDATA-half_window,
    #   set mask to True if std of window is greater than std_limit
    for i in range(half_window,nDATA-half_window+1):
        mask[i] = np.std(data[i-half_window:i+half_window]) > std_limit


    # fill start and end of array with first and last calculated values.
    #   Could add other options for dealing with ends if desired
    mask[:half_window]       = mask[half_window]
    mask[nDATA-half_window:] = mask[nDATA-half_window]
    
    return mask


def meanbox_iterate(data,window_radius=1,iterations=1):
    
    npts=len(data)
    for i_iter in range(iterations):
        new_data=np.copy(data)
        #print(new_data[:5])
        for i_pt in range(1,npts-1):
            if i_pt<window_radius:
                #print(i_pt)
                #print(data[:i_pt+window_radius])
                new_data[i_pt]=np.mean(data[:i_pt+window_radius+1])# 
            elif i_pt>(npts-window_radius):
                new_data[i_pt]=np.mean(data[i_pt-window_radius:])
            else:
                new_data[i_pt]=np.mean(data[i_pt-window_radius:i_pt+window_radius+1])
        data=np.copy(new_data)
    
    return new_data

