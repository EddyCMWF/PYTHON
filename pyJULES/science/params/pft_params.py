#############################################################################
#
#   JULES pft_params.py
#   
#   Owner: Edward Comyn-Platt
#   
#   Pupose: Various pft_params configurations
#
#   Description: Based on various JULES configurations
#
#############################################################################


#########################################################################
# Function: read_pft_params_nml
# Author: Edward Comyn-Platt
# Purpose: Read a jules pft_params.nml
# Inputs: 
# Optional Inputs - file = 'pft_params.nml' (in current working dir.
#                   PFT_names = Names for pfts, if not given they are indexed 
#                               1 to Npfts in file
# Outputs: PFT_params - Dictionary containing all the PFT params.
#########################################################################
def read_pft_params_nml(nml_file='pft_params.nml'):
    # Open file and read lines
    inlines=open(nml_file,'r').readlines()
    inlines = filter(None,[ line.replace('\n','') for line in inlines ])
    # get rid of header lines which start with a '#' or '&', or
    #   have zero length when spaces are removed
    headers=[]
    while inlines[0]!='&JULES_PFTPARM':
        headers.append(inlines.pop(0))
    # pop first line:
    headers.append(inlines.pop(0))
        
    pft_param_dict={}
    # '/' signifies end of namelist
    while inlines[0]!='/':
        line=inlines.pop(0)
        if line[0]!='!':
            split1=line.replace(' ','').split('=')
            param=split1[0]
            vals=filter(None,split1[1].split(','))
            print(param,vals)
            pft_param_dict[param]=[float(val) for val in vals]
    
    return pft_param_dict





#def JULES_C_BL():
#    pft_params = 

#def JULES_C_5PFTS():
    
