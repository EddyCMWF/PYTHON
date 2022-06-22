#!/usr/local/sci/bin/python2.7
# -*- coding: iso-8859-1 -*-
#import iris
import numpy as np
import pandas as pd
import xarray as xr
#import statsmodels.api as sm
#import matplotlib.pyplot as plt
#import pickle as pkl
#from netCDF4 import Dataset
#import sys
#from scipy import stats
from copy import deepcopy

class hwsd():
    def __init__(self, latrangeOut=(-90,90), lonrangeOut=(-180,180), 
                       latresOut=0.5, lonresOut=0.5,
                       output_vars = ['T_SAND','T_SILT','T_CLAY','T_OC']):

#        self.fill = True # fill in holes
        self.fill = False

#        self.hwsdFile = "hwsd.bil"
        #self.hwsdFile = "/project/jchmr/hadng/hwsd/hwsd_rjel.bil"
        self.hwsdFile = "/local/localscratch/edwcom/data/HWSD/hwsd.bil"
        self.hwsdLUTfile = "/local/localscratch/edwcom/data/HWSD/HWSD_DATA.txt"
        
        # In grid parameters
        self.nyIn = 21600
        self.latminIn = -90.
        self.latmaxIn = 90.
        self.latrangeIn = (self.latminIn, self.latmaxIn)
        self.latresIn = np.abs(self.latmaxIn-self.latminIn)/self.nyIn
        
        self.nxIn = 43200
        self.lonminIn = -180.
        self.lonmaxIn = 180.
        self.lonrangeIn = (self.lonminIn, self.lonmaxIn)
        self.lonresIn = np.abs(self.lonmaxIn-self.lonminIn)/self.nxIn
        
        # Out grid parameters:
        if latresOut=='native':
            self.latresOut = self.latresIn
        else:
            self.latresOut = latresOut
        self.set_latrangeOut(latrangeOut)
        
        if lonresOut=='native':
            self.lonresOut = self.lonresIn
        else:
            self.lonresOut = lonresOut
        self.set_lonrangeOut(lonrangeOut)     
        
        self.xstartIn = int((self.lonminOut - self.lonminIn)/self.lonresIn)
        self.ystartIn = int((self.latminOut - self.latminIn)/self.latresIn)
        self.xendIn = int((self.lonmaxOut - self.lonminIn)/self.lonresIn)
        self.yendIn = int((self.latmaxOut - self.latminIn)/self.latresIn)
        
        self.dy = int(float(self.latresOut)/float(self.latresIn))
        self.dx = int(float(self.lonresOut)/float(self.lonresIn))
        
        self.dxy = 5 # set value for infilling. 5= normal
        
        self.fill_value=-9999
        
        self.output_vars = output_vars
        self.n_output_vars = len(output_vars)

    def set_latrangeOut(self,latrange):
        self.latrangeOut = latrange
        self.latminOut = self.latrangeOut[0]
        self.latmaxOut = self.latrangeOut[1]
        self.nyOut = int( (self.latmaxOut-self.latminOut)/self.latresOut )
      
    def set_lonrangeOut(self,lonrange):
        self.lonrangeOut = lonrange
        self.lonminOut = self.lonrangeOut[0]
        self.lonmaxOut = self.lonrangeOut[1]
        self.nxOut = int( (self.lonmaxOut-self.lonminOut)/self.lonresOut )
      
    def set_latresOut(self,latres):
        self.latresOut = latres
        self.nyOut = int( (self.latmaxOut-self.latminOut)/self.latresOut )
      
    def set_lonresOut(self,lonres):
        self.lonresOut = lonres
        self.nxOut = int( (self.lonmaxOut-self.lonminOut)/self.lonresOut )
        
    def geoGrid(self):
        self.lonsOut = np.arange(self.lonminOut,self.lonmaxOut,self.lonresOut)+self.lonresOut/2.
        self.latsOut = np.arange(self.latminOut,self.latmaxOut,self.latresOut)+self.lonresOut/2.
    
    def get_geoGrid(self):
        if not hasattr(hwsd,'lonsOut'): self.geoGrid()
        return self.lonsOut, self.latsOut
    
    def original_output_vars(self):
        # Output vars from the script I orignally received from Nic.
        self.output_vars = ['T_SAND','T_SILT','T_CLAY','T_OC','S_OC','T_GRAVEL','S_GRAVEL',
                    'T_BULK_DENSITY','S_BULK_DENSITY','T_REF_BULK_DENSITY','S_REF_BULK_DENSITY',
                    'T_MOD_BULK_DENSITY','S_MOD_BULK_DENSITY',
                    'TcarbMass','TcarbMassRef','TcarbMassMod',
                    'ScarbMass','ScarbMassRef','ScarbMassMod']
##
    def readBIL(self):

        print("starting to read bil file: "+self.hwsdFile)
        fid = open( self.hwsdFile, "rb")
        ## file seems to be little endian 16 bit integers
        self.hwsdGrid = np.fromfile( fid, 
                                     count = self.nxIn * self.nyIn,
                                     dtype = np.int16)
        ## now reshape to the grid size
        self.hwsdGrid = self.hwsdGrid.reshape(( self.nyIn, self.nxIn))
        ## Flip north and south:
        self.hwsdGrid = self.hwsdGrid[::-1,:]
        # Extract the section to be regridded:
        self.hwsdGrid = self.hwsdGrid[self.ystartIn:self.yendIn,self.xstartIn:self.xendIn]
        
        fid.close()
##

#    def create_output_arrays(self):
#        self.fao90code = np.zeros( ( self.nyOut, self.nxOut), np.str)
#        
#        self.outputData = { var: np.ma.masked_equal( 
#                np.zeros((self.nyOut,self.nxOut),np.float)+self.fill_value, self.fill_value )
#                        for var in self.output_vars }
    
    def create_output_Xarrays(self):
        self.fao90code = np.zeros( ( self.nyOut, self.nxOut), np.str)
        if not hasattr(hwsd,'lonsOut'): self.geoGrid()
        #empty_array = np.ma.masked_equal(
        #        np.zeros((self.nyOut,self.nxOut),np.float)+self.fill_value,self.fill_value)
        empty_array = np.zeros((self.nyOut,self.nxOut),np.float)+np.nan
        self.outputXRData = xr.Dataset( 
                { var: (['lat','lon'], deepcopy(empty_array) )
                        for var in self.output_vars },
                coords={'lat':(['lat'],self.latsOut),
                        'lon':(['lon'],self.lonsOut)}
                )
        
        self.outputXRData['lat'].attrs['long_name']='Latitude'
        self.outputXRData['lat'].attrs['units']='Degrees North'
        self.outputXRData['lon'].attrs['long_name']='Longitude'
        self.outputXRData['lon'].attrs['units']='Degrees East'
##
    def nativeHWSD(self):
        #self.create_output_arrays()
        self.create_output_Xarrays()
        
        min_var = self.hwsdGrid[self.hwsdGrid!=0].min()
        max_var = self.hwsdGrid[self.hwsdGrid!=0].max()
        for isoil in xrange(min_var,max_var+1):
            mask = self.hwsdGrid==isoil
            
            if mask.sum()>0:
                self.fao90code[mask] = self.LUT['fao90'][isoil]
                for var in self.output_vars:
                    #self.outputData[var][mask] = self.LUT[var][isoil]
                    self.outputXRData[var].data[mask] = self.LUT[var][isoil]
        
        for var in self.output_vars:
            self.outputXRData[var].data[self.outputXRData[var].data==self.fill_value]=np.nan
        
    def regridHSWD(self):
        
        #self.create_output_arrays()
        self.create_output_Xarrays()
        
        dy = self.dy
        dx = self.dx
        
        for y in xrange( self.nyOut):
## how many pixels to miss out before starting in the vertical
            #y0 = self.ystartIn + (y * dy)   # Now extract a smaller original grid so start from edge
            y0 = y * dy
            y1 = y0 + dy
            for x in xrange( self.nxOut):
                # x0 = self.xstartIn + (x * dx)
                x0 = x * dx
                x1 = x0 + dx
                # Subsection of hwsd to be averaged:
                var2Mean = self.hwsdGrid[ y0:y1, x0:x1]
                ii = var2Mean > 0
                if ii.sum() > 0:
                    # lowest and highest indexes soil classes:
                    min_var=min(var2Mean[ii])
                    max_var=max(var2Mean[ii])

                    binb=min_var
                    bine=max_var+1
                    # edges of bins for np.hist (set to 0.1 less than values to avoid confusion):
                    bins=np.arange(binb,bine+1)-0.1  
                    # actual values for bins/var2mean
                    bins_act=np.arange(binb,bine)
                    # use histogram for speed
                    varfreq,varval = np.histogram(var2Mean[ii],bins=bins) 
                        
                    # remove those which have zero occurance:
                    iii = varfreq > 0
                    text_sbins=bins_act[iii]
                    #text_svarval=varval[iii]    # Numpy must have changed as this no longer works
                    #text_svarval=varval[:-1][iii]   # Use this instead of the above
                    text_svarfreq=varfreq[iii]
                    
                    #remove those which have no texture values:                       
                    soilText = np.zeros(len(text_sbins),dtype=int)
                    for ic,sbin in enumerate(text_sbins):
                        soilText[ic] = ( self.LUT['T_SAND'][sbin]+self.LUT['T_SILT'][sbin]
                                        +self.LUT['T_CLAY'][sbin] )
                    iv = soilText > 0.0
                    if iv.sum() > 0:
                        sbins=text_sbins[iv]
                        #svarval=text_svarval[iv]
                        svarfreq=text_svarfreq[iv]
    
                        # Find dominant soil group, now stored as string:
                        soilGroups = self.LUT['fao90'][sbins]
                        soilGroups_unique = self.LUT['fao90'][sbins].unique()
                        countSoilGroup = np.zeros_like(soilGroups_unique)
                        for ic,SG in enumerate(soilGroups_unique):
                            ind = np.where(soilGroups == soilGroups_unique[ic])[0]
                            # now count how many times this soil group occurs:
                            countSoilGroup[ic] = svarfreq[ind].sum()
                        
                        # In case of multiple occurrences of the maximum values, the indices
                        #  corresponding to the first occurrence are returned.
                        i_soilMain = np.argmax(countSoilGroup)
                        soilMain = soilGroups_unique[i_soilMain]
                        self.fao90code[ y, x] = soilMain[0]

                        #######################
                        ## Calculate avg of all soils in grid boxes:
                        sLUT = self.LUT.loc[sbins][self.output_vars]
                        sLUT = sLUT.mask(sLUT<0.)
                        
                        #avg_dictionary = {}
                        for var in self.output_vars:
                            svarfreq_avg = deepcopy(svarfreq)
                            # Remove any non available soil properties
                            svarfreq_avg[sLUT[var].mask==True]=0
                            scount = svarfreq_avg.sum()
                            svar_weight = svarfreq_avg.astype(float)/float(scount)
                            #self.outputData[var][ y, x] = (sLUT[var]*svar_weight).sum()
                            self.outputXRData[var].data[ y, x] = (sLUT[var]*svar_weight).sum()
                            
        for var in self.output_vars:
            self.outputXRData[var].data[self.outputXRData[var].data==self.fill_value]=np.nan


##
    def readLookUpGroups(self):
# Routine to find and label the major soil groups

        fid = open( self.hwsdLUTfile, "r") 
        #fid = open( "/local/localscratch/edwcom/data/HWSD/HWSD_DATA.txt","r") ## HWSD lookup table
        LUT_lines = fid.readlines()
        fid.close()
        
        #header = fid.readline()
        header = LUT_lines.pop(0)
        header_cols =  header.replace('"',"").replace("\r\n","").split(",")
        LUT_dictionary = { hdr_col:[] for hdr_col in header_cols }
        
        for line in LUT_lines:
            split=line.replace('"','').replace("\r\n","").split(",")
            if split[1] not in LUT_dictionary['MU_GLOBAL']:     # Only count the first record for each MU_GLOBAL 
                for i,hdr_col in enumerate(header_cols):
                    LUT_dictionary[hdr_col].append(split[i])
        
        for hdr_col in header_cols: 
            LUT_dictionary[hdr_col] = np.array(LUT_dictionary[hdr_col])
        
        int_vars = ["MU_GLOBAL","IL","ID","ROOTS","DRAINAGE","REF_DEPTH",
                    "PHASE1","PHASE2","SEQ","SWR","ADD_PROP","AWC_CLASS",
                    "T_SAND","T_SILT","T_CLAY","T_GRAVEL","T_TEXTURE","T_USDA_TEX_CLASS",
                    "S_SAND","S_SILT","S_CLAY","S_GRAVEL","S_USDA_TEX_CLASS",
                    ]
                    # "S_TEXTURE",
        for var in int_vars: LUT_dictionary[var] = LUT_dictionary[var].astype(int)
        
        flt_vars = ["SHARE",
                    "T_OC","T_BULK_DENSITY","T_REF_BULK_DENSITY","T_BS",
                    "T_ESP","T_ECE","T_TEB",
                    "T_CEC_SOIL","T_CEC_CLAY",
                    "T_PH_H2O","T_CACO3","T_CASO4",
                    "S_OC","S_BULK_DENSITY","S_REF_BULK_DENSITY","S_BS",
                    "S_ESP","S_ECE","S_TEB",
                    "S_CEC_SOIL","S_CEC_CLAY",
                    "S_PH_H2O","S_CACO3","S_CASO4",
                    ]
        for var in flt_vars: LUT_dictionary[var] = LUT_dictionary[var].astype(float)
        
        bool_vars = ["ISSOIL",
                     ]
        for var in bool_vars: LUT_dictionary[var] = LUT_dictionary[var].astype(bool)
        
        str_vars = [ "SU_SYM74",  "SU_SYM85",  "SU_SYM90",
                     "SU_CODE74", "SU_CODE85", "SU_CODE90",
                    ]
        # Not required as already strings:
        #for var in str_vars: LUT_dictionary[hdr_col] = LUT_dictionary[hdr_col].astype(float)
        
        # Checksum for the sand silt clay fractions:
        CSUM = LUT_dictionary['T_SAND']+LUT_dictionary['T_SILT']+LUT_dictionary['T_CLAY']
        bad_CSUM = np.where((CSUM != 0)&(CSUM != 100)&(CSUM != -29997))[0]
        for i in bad_CSUM: print('Check 100% soil textures: ',i,
                              LUT_dictionary['T_SAND'][i],LUT_dictionary['T_SILT'][i],
                              LUT_dictionary['T_CLAY'][i],CSUM[i])
        
        
        # Get the FAO code and fill in where necessary:        
        fao90=deepcopy(LUT_dictionary["SU_SYM90"])
        # Correct error in infile:
        fao90[fao90 == "1.45"] = "RGe"

        ## where there is no code, take appropriate code in FAO74 and do conversion
        # Conversion arrays:
        sym74=[ "G", "R", "I", "Q", "E", "U", "T", "V", "Z", "S", "Y" 
                   , "X", "K", "C", "H", "M", "B", "L", "D", "P", "W" 
                   , "A", "N", "F", "O", "J"]
        sym90=["GL","RG","LP","AR","-1","-2", "AN","VR","SC","SN","-3" 
                   ,"-4","KS","CH","PH","GR","CB","LV","PD","PZ", "PL"
                   ,"AC","NT","FR","HS","FL"]
        
        no_fao90_code = np.where(fao90 == "-9999")[0]
        
        for icode in no_fao90_code:
            fao74_code = LUT_dictionary["SU_SYM74"][icode][0]
            ## just take first letter as this is all that's used in FAO74 soil groups:
            ind=sym74.index(fao74_code)
            fao90[icode]=sym90[ind]

        # ensure only have first two characters:
        for icode in range(len(fao90)):
            if len(fao90[icode])>2:
                fao90[icode] = fao90[icode][:2]
        #fao90_list = fao90.tolist()
        LUT_dictionary['fao90'] = fao90
        
        
        # Calulate the Carbon contents by mass:
        ##############################################
        TcarbPerc=LUT_dictionary['T_OC']
        ScarbPerc=LUT_dictionary['S_OC']
        
        # Get gravel content, and fill gaps with zeros.
        TGr = deepcopy(LUT_dictionary['T_GRAVEL'])
        TGr[TGr<0.] = 0.
        SGr = deepcopy(LUT_dictionary['S_GRAVEL'])
        SGr[SGr<0.] = 0.
        
        # Bulk Density:
        TBD=LUT_dictionary['T_BULK_DENSITY']
        Tmask = (TcarbPerc < 0.)|(TBD<0.)
        LUT_dictionary['TcarbMass'] = TcarbPerc/100.*TBD*1000.0*(1.0-0.01*TGr)
        LUT_dictionary['TcarbMass'][Tmask] = self.fill_value
        SBD=LUT_dictionary['S_BULK_DENSITY']
        Smask = (ScarbPerc < 0.)|(SBD<0.)
        LUT_dictionary['ScarbMass'] = ScarbPerc/100.*SBD*1000.0*(1.0-0.01*SGr)
        LUT_dictionary['ScarbMass'][Smask] = self.fill_value
        
        # Refernce Bulk Density:
        TRefBD=LUT_dictionary['T_REF_BULK_DENSITY']
        TmaskRef = (TcarbPerc < 0.)|(TRefBD<0.)
        LUT_dictionary['TcarbMassRef'] = TcarbPerc/100.*TRefBD*1000.0*(1.0-0.01*TGr)
        LUT_dictionary['TcarbMassRef'][TmaskRef] = self.fill_value
        SRefBD=LUT_dictionary['S_REF_BULK_DENSITY']
        SmaskRef = (ScarbPerc < 0.)|(SRefBD<0.)
        LUT_dictionary['ScarbMassRef'] = ScarbPerc/100.*SRefBD*1000.0*(1.0-0.01*SGr)
        LUT_dictionary['ScarbMassRef'][SmaskRef] = self.fill_value
        
        # Modified carbon mass treats the Histo- and ferro-sols differently.
        ind = (fao90 == 'HS') | (fao90 == 'AN')
        TModBD = deepcopy(TRefBD)
        TModBD[ind] = TBD[ind]
        TmaskMod = (TcarbPerc < 0.)|(TModBD<0.)
        LUT_dictionary['T_MOD_BULK_DENSITY']=TModBD
        LUT_dictionary['TcarbMassMod'] = TcarbPerc/100.*TModBD*1000.0*(1.0-0.01*TGr)
        LUT_dictionary['TcarbMassMod'][TmaskMod] = self.fill_value
        SModBD = deepcopy(SRefBD)
        SModBD[ind] = SBD[ind]
        SmaskMod = (ScarbPerc < 0.)|(SModBD<0.)
        LUT_dictionary['S_MOD_BULK_DENSITY']=SModBD
        LUT_dictionary['ScarbMassMod'] = ScarbPerc/100.*SModBD*1000.0*(1.0-0.01*SGr)
        LUT_dictionary['ScarbMassMod'][SmaskMod] = self.fill_value
        
        # Store the LUT in the class as a panda dataframe:
        self.LUT = pd.DataFrame(LUT_dictionary,index=LUT_dictionary['MU_GLOBAL'])
        
    #I= Lithosols = very limited by depth over hard rock
    
    # Y & X yermosols  and xerosols "A major change to the legend has been
    # the removal of two first level classes that were defined by an aridic
    # soil moisture regime, Yermosols and Xerosols. This change was based on
    # one of FAO's general principles of their classification system, "that no
    # climatic criteria would be used to define the soil units." The two
    # classes were originally established because there were no better
    # separation criteria. Accumulation of calcium carbonate and gypsum
    # are now used as additional separation criteria to deal with the
    # aridic problem. Calcisols and Gypsisols classes were introduced for
    # this purpose. These soils occur predominately under arid and
    # semi-arid conditions (FAO, 1988, p. 5-6)."
    # Therefore in sym90 could just assume they are calcisols CL
    # LP=lithic leptosols = very limited by depth over hard rock


#    def fillHoles(self,fill_vars=None):
##        dxy = 5
#        dxy = self.dxy
#        # If not specified, fill in all variables:
#        if fill_vars == None:
#            fill_vars = self.output_vars
#        
#        for var in fill_vars:
#            temp_data = deepcopy(self.outputXRData[var])
#            
#            bad_mask = temp_data==self.fill_value
#            bad_points = np.where(bad_mask)
#            bad_points = np.where(~bad_mask)
#            
#            for bdpt in bad_points:
#                
#            
#        temp_fld = deepcopy(fld)
#        for it in xrange(1000):
#
#            ii = ( self.mask > 0) & ( fld < 0) 
#            if ii.sum() < 1:
#                print( "Filled all the holes for fld!",it)
#                break
#            print( it, dxy, ii.sum())
#            for y in xrange( self.nyOut):
#                for x in xrange( self.nxOut):
#                    if ii[y,x] > 0:
#                        ind =temp_fld[y-dxy:y+dxy,x-dxy:x+dxy] >= 0
#                        if ind.sum() > 0 :
#                            vfld = temp_fld[y-dxy:y+dxy,x-dxy:x+dxy][ind]
#                            fld[y,x] = np.mean(vfld)
#
##            dxy += 10
#            dxy += 2*self.dxy
#
#
#        return fld
##

    def doStuff(self):
        # Populate all the fields using default set up.
        self.readLookUpGroups()
        
        # read in the raw 16 bit integer file with the soil codes
        self.readBIL()
        
        if (self.latresIn==self.latresOut)&(self.lonresIn==self.lonresOut):
            # if new resolution same as HWSD resolution, just extract data.
            self.nativeHWSD()
        else:
            # elsi Regrid to new dimensions
            self.regridHSWD()
        
        # set up lat lon grids.
        self.geoGrid()
        
        
