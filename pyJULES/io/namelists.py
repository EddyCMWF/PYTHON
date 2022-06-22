###############################################################################
# Module: namelist
# Owner: Edward Comyn-Platt
# Purpose: Module to handle a JULES namelists
###############################################################################

###############################################################################
# Function: read_pft_params
# Author: Edward Comyn-Platt
# Purpose: to read a pft_params namelist
# Inputs: None
# Optional Inputs: nml_dir='./' - namelist directory, 
#                  filename= 'pft_params.nml' - pft_params filename,
# Output: PFT_PARAMS - dictioonary containing the pft_params
###############################################################################
def read_pft_params(nml_dir='.',filename='pft_params.nml'):
    # Read lines of pft file
    lines=[line.rstrip() for line in open(nml_dir+'/'+filename)] 
    #Create dictionary to store params
    PFT_PARAMS={}
    #loop araound lines
    for line in lines:
        # ignore lines that do not contain an =
        if '=' in line:
            # Remove ',' from the end of the line
            if line[-1]==',':
                line=line[:-1]
            
            # split line into name and values string
            name,vals_string=line.split('=')
            # split vals_string
            vals_string=vals_string.split(',')
            
            # loop arounds vals_string, appending values as floats to list
            vals=[]
            for val_s in vals_string:
                # deal with multiple assignment
                if '*' in val_s:
                    N,val=val_s.split('*')
                    vals=vals+[float(val) for i in range(int(N))]
                else:
                    vals=vals+[float(val_s)]
            # apend vals to PFT_PARAMS dictionary
            PFT_PARAMS[name]=list(vals)

    return PFT_PARAMS

###############################################################################
# Function: read_default_pft_params
# Author: Edward Comyn-Platt
# Purpose: to read a pft_params namelist
# Inputs: None
# Output: PFT_PARAMS - dictioonary containing the pft_params
###############################################################################
def read_default_pft_params():
    import os
    nml_dir=os.path.dirname(__file__)+'/default_namelists'
    print nml_dir 
    return read_pft_params(nml_dir=nml_dir)




