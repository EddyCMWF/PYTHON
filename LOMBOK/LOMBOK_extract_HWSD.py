#!/usr/local/sci/bin/python2.7
# -*- coding: iso-8859-1 -*-
#import iris
import numpy as np
#import statsmodels.api as sm
import matplotlib.pyplot as plt
#import pickle as pkl
from netCDF4 import Dataset
#import sys
#from scipy import stats
from copy import deepcopy

class test():
    def __init__(self):

#        self.fill = True # fill in holes
        self.fill = False

#        self.hwsdFile = "hwsd.bil"
        self.hwsdFile = "/project/jchmr/hadng/hwsd/hwsd_rjel.bil"
        self.nyIn = 21600
        self.nxIn = 43200

        self.nyOut = 360
        self.nxOut = 720

        self.dxy = 5 # set value for infilling. 5= normal
#       self.dxy = 30 # set value for infilling. 5= normal
##
##
    def readBIL(self):

        print("starting to read bil file")
        fid = open( self.hwsdFile, "rb")
## file seems to be little endian 16 bit integers
        self.hwsdGrid = np.fromfile( fid, \
                                     count = self.nxIn * self.nyIn,\
                                     dtype = np.int16)
## now reshape to the grid size
        self.hwsdGrid = self.hwsdGrid.reshape(( self.nyIn, self.nxIn))
##
##
    def dataAvg(self,fldout,ind,fldin,svarfreq):

# Average the data within each coarse grid box:

        svarfreqUse = deepcopy(svarfreq)
        fldinUse = fldin
# if missing data (set by ind) then remove from averaging:
        if ind.sum() > 0:
            fldinUse[ind] = 0.
            svarfreqUse[ind] = 0 

        scountTot = svarfreqUse.sum()
        if scountTot > 0:
            fldout = np.array(fldinUse*svarfreqUse/scountTot).sum()

        return fldout
##
##
    def regridHSWD(self, dx, dy):
        self.newGrid = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.int16)
        self.sand = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.silt = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.clay = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TcarbPerc = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.ScarbPerc = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TcarbMass = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.ScarbMass = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TcarbMassRef = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.ScarbMassRef = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TcarbMassMod = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.ScarbMassMod = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.SBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TRefBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.SRefBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TModBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.SModBD = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.TGr = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)
        self.SGr = np.zeros( ( self.nyOut, self.nxOut),\
                                 np.float)

        self.sand[:] = -999
        self.silt[:] = -999
        self.clay[:] = -999
        self.TcarbPerc[:] = -999
        self.ScarbPerc[:] = -999
        self.TcarbMass[:] = -999
        self.ScarbMass[:] = -999
        self.TcarbMassRef[:] = -999
        self.ScarbMassRef[:] = -999
        self.TcarbMassMod[:] = -999
        self.ScarbMassMod[:] = -999
        self.TBD[:] = -999
        self.SBD[:] = -999
        self.TRefBD[:] = -999
        self.SRefBD[:] = -999
        self.TModBD[:] = -999
        self.SModBD[:] = -999
        self.TGr[:] = -999
        self.SGr[:] = -999

##
        ferrasol = self.fao90_index('FR')

        for y in xrange( self.nyOut):
## how many pixels to miss out before starting in the vertical
            y0 = y * dy
            y1 = ( y + 1) * dy
            for x in xrange( self.nxOut):
                x0 = x * dx
                x1 = ( x + 1) * dx
                var2Mean = self.hwsdGrid[ y0:y1, x0:x1]
                ii = var2Mean > 0
                if ii.sum() > 0:
                    min_var=min(var2Mean[ii])
                    max_var=max(var2Mean[ii])

                    binb=min_var
                    bine=max_var+1
# edges of bins for np.hist (set to 0.1 less than values to avoid confusion) :
                    bins=np.arange(binb,bine+1)-0.1  
# actual values for bins/var2mean
                    bins_act=np.arange(binb,bine)
# use histogram for speed
                    varfreq,varval = np.histogram(var2Mean[ii],bins=bins) 
                        
# remove those which have zero occurance:
                    iii = varfreq > 0
                    text_sbins=bins_act[iii]
                    text_svarval=varval[iii]
                    text_svarfreq=varfreq[iii]

                    soilText = np.zeros(len(text_sbins),dtype=int)
                    for ic in xrange(0,len(text_sbins)) :
                        soilText[ic] = int(self.Textures[str(text_sbins[ic])][0]) + \
                            int(self.Textures[str(text_sbins[ic])][1]) + int(self.Textures[str(text_sbins[ic])][2])

#remove those which have no texture values:                       
                    iv = soilText > 0.0
                    if iv.sum() > 0:
                        sbins=text_sbins[iv]
                        svarval=text_svarval[iv]
                        svarfreq=text_svarfreq[iv]
                        
# get corresponding soil groups for all mu numbers with a non-zero freq
                        soilGroup = np.zeros(len(sbins),dtype=int)
                        for ic in xrange(0,len(sbins)) :
                            soilGroup[ic] = self.soilInd[sbins[ic]]

                        fao90_list = []
                        countSoilGroup = [0]*len(sbins)
                        for ic in xrange(0,len(sbins)) :

                            old = soilGroup[ic] in fao90_list
                            if old == False :
                                ind = soilGroup - soilGroup[ic] == 0
# now count how many times this soil group occurs:
                                countSoilGroup[ic] = svarfreq[ind].sum() 
                            fao90_list.append(soilGroup[ic])

# Find dominant soil group:
                        max_val = countSoilGroup == max(countSoilGroup)
                        soilMain=soilGroup[max_val]
                        if len(soilMain) > 1 :
                            print("warning soilMain=",soilMain)
                            soilMain=soilMain[0] # just take first value for now

                        self.newGrid[ y, x] = soilMain

#######################
## Calculate avg of all soils in grid boxes:
                        svarfreq = np.array(svarfreq) # convert from list to array so can index multiple values in one
                        scountSoilGroup = svarfreq.sum()

                        sand = np.zeros(len(sbins),dtype=float)
                        silt = np.zeros(len(sbins),dtype=float)
                        clay = np.zeros(len(sbins),dtype=float)
                        TcarbPerc = np.zeros(len(sbins),dtype=float)
                        ScarbPerc = np.zeros(len(sbins),dtype=float)
                        TcarbMass = np.zeros(len(sbins),dtype=float)
                        ScarbMass = np.zeros(len(sbins),dtype=float)
                        TcarbMassRef = np.zeros(len(sbins),dtype=float)
                        ScarbMassRef = np.zeros(len(sbins),dtype=float)
                        TcarbMassMod = np.zeros(len(sbins),dtype=float)
                        ScarbMassMod = np.zeros(len(sbins),dtype=float)
                        TBD = np.zeros(len(sbins),dtype=float)
                        SBD = np.zeros(len(sbins),dtype=float)
                        TModBD = np.zeros(len(sbins),dtype=float)
                        SModBD = np.zeros(len(sbins),dtype=float)
                        TRefBD = np.zeros(len(sbins),dtype=float)
                        SRefBD = np.zeros(len(sbins),dtype=float)
                        TGr = np.zeros(len(sbins),dtype=float)
                        SGr = np.zeros(len(sbins),dtype=float)

                        for ic in xrange(0,len(sbins)) :
                            sand[ic] = self.Textures[str(sbins[ic])][0]
                            silt[ic] = self.Textures[str(sbins[ic])][1]
                            clay[ic] = self.Textures[str(sbins[ic])][2]
                            TcarbPerc[ic] =self.Textures[str(sbins[ic])][3]
                            ScarbPerc[ic] =self.Textures[str(sbins[ic])][7]
                            
                            TBD[ic] = self.Textures[str(sbins[ic])][6]
                            SBD[ic] = self.Textures[str(sbins[ic])][8]
                            TRefBD[ic] = self.Textures[str(sbins[ic])][9]
                            SRefBD[ic] = self.Textures[str(sbins[ic])][10]
                            TGr[ic] = self.Textures[str(sbins[ic])][11]
                            SGr[ic] = self.Textures[str(sbins[ic])][12]
 
                        ind = TGr < 0
                        if ind.sum() > 0:
                            TGr[ind] = 0
                        ind = SGr < 0
                        if ind.sum() > 0:
                            SGr[ind] = 0

                        TModBD = deepcopy(TRefBD)
                        SModBD = deepcopy(TRefBD)
                        ind = (soilGroup == self.histosol) | (soilGroup == self.andosol)
                        if ind.sum() > 0:
                            TModBD[ind] = TBD[ind]
                            SModBD[ind] = TBD[ind]

                        TcarbMass = TcarbPerc/100.*TBD*1000.0*(1.0-0.01*TGr)
                        ScarbMass = ScarbPerc/100.*SBD*1000.0*(1.0-0.01*SGr)
                        TcarbMassRef = TcarbPerc/100.*TRefBD*1000.0*(1.0-0.01*TGr)
                        ScarbMassRef = ScarbPerc/100.*SRefBD*1000.0*(1.0-0.01*SGr)
                        TcarbMassMod = TcarbPerc/100.*TModBD*1000.0*(1.0-0.01*TGr)
                        ScarbMassMod = ScarbPerc/100.*SModBD*1000.0*(1.0-0.01*SGr)
                        
# Calculate the field before averaging takes place:

                        ind = (TcarbPerc < 0) | (TBD < 0) | (TGr < 0)
                        self.TcarbMass[y,x]=self.dataAvg(self.TcarbMass[y,x],ind,TcarbMass,svarfreq)
                        ind = (ScarbPerc < 0) | (SBD < 0) | (SGr < 0)
                        self.ScarbMass[y,x]=self.dataAvg(self.ScarbMass[y,x],ind,ScarbMass,svarfreq)

                        ind = (TcarbPerc < 0) | (TRefBD < 0) | (TGr < 0)
                        self.TcarbMassRef[y,x]=self.dataAvg(self.TcarbMassRef[y,x],ind,TcarbMassRef,svarfreq)
                        ind = (ScarbPerc < 0) | (SRefBD < 0) | (SGr < 0)
                        self.ScarbMassRef[y,x]=self.dataAvg(self.ScarbMassRef[y,x],ind,ScarbMassRef,svarfreq)

                        ind = (TcarbPerc < 0) | (TModBD < 0) | (TGr < 0)
                        self.TcarbMassMod[y,x]=self.dataAvg(self.TcarbMassMod[y,x],ind,TcarbMassMod,svarfreq)
                        ind = (ScarbPerc < 0) | (SModBD < 0) | (SGr < 0)
                        self.ScarbMassMod[y,x]=self.dataAvg(self.ScarbMassMod[y,x],ind,ScarbMassMod,svarfreq)

                            
# This must be done after tcarbmass* and scarbmass calcs above:
                        ind = (TcarbPerc < 0) 
                        self.TcarbPerc[y,x]=self.dataAvg(self.TcarbPerc[y,x], TcarbPerc < 0 ,TcarbPerc,svarfreq)
                        ind = (ScarbPerc < 0) 
                        self.ScarbPerc[y,x]=self.dataAvg(self.ScarbPerc[y,x], ScarbPerc < 0 ,ScarbPerc,svarfreq)

                        self.TBD[y,x]=self.dataAvg(self.TBD[y,x], TBD < 0 ,TBD,svarfreq)
                        self.SBD[y,x]=self.dataAvg(self.SBD[y,x], SBD < 0 ,SBD,svarfreq)

                        self.TRefBD[y,x]=self.dataAvg(self.TRefBD[y,x], TRefBD < 0 ,TRefBD,svarfreq)
                        self.SRefBD[y,x]=self.dataAvg(self.SRefBD[y,x], SRefBD < 0 ,SRefBD,svarfreq)

                        self.TModBD[y,x]=self.dataAvg(self.TModBD[y,x], TModBD < 0 ,TModBD,svarfreq)
                        self.SModBD[y,x]=self.dataAvg(self.SModBD[y,x], SModBD < 0 ,SModBD,svarfreq)


                        self.sand[ y, x] = np.array(sand*svarfreq/scountSoilGroup).sum()
                        self.silt[ y, x] = np.array(silt*svarfreq/scountSoilGroup).sum()
                        self.clay[ y, x] = np.array(clay*svarfreq/scountSoilGroup).sum()

                        
# Check there are no nans:
                        ind = np.isnan(self.TcarbMass[ y, x]) | np.isnan(self.ScarbMass[ y, x])
                        if ind.sum() > 0:
                            print("nans Mass",self.TcarbMass[ y, x],self.ScarbMass[ y, x])
                            print(ScarbMass,svarfreqSCarbMass)
                            print(scountSoilGroupSCarbMass)
                        ind = np.isnan(self.TcarbMassRef[ y, x]) | np.isnan(self.ScarbMassRef[ y, x])
                        if ind.sum() > 0:
                            print("nans MassRef",self.TcarbMassRef[ y, x],self.ScarbMassRef[ y, x])
                        ind = np.isnan(self.TcarbMassMod[ y, x]) | np.isnan(self.ScarbMassMod[ y, x])
                        if ind.sum() > 0:
                            print("nans MassMod",self.TcarbMassMod[ y, x],self.ScarbMassMod[ y, x])


# check sand silt and clay total= 100
                        checkSum = self.sand[ y, x]+self.silt[ y, x]+self.clay[ y, x]
                        if checkSum > 0:
                            if checkSum < 99.9:
                                print("min max",checkSum.min(),checkSum.max())
                                print(" self.sand[ y, x]", self.sand[ y, x])
                                print(" self.silt[ y, x]", self.silt[ y, x])
                                print(" self.clay[ y, x]", self.clay[ y, x])
                                print("sand",sand[:])
                                print("silt",silt[:])
                                print("clay",clay[:])
                                print("svarfreq",svarfreq)
                                print("scountSoilGroup",scountSoilGroup)
                                quit()
                            if checkSum > 100.1:
                                print("min max",checkSum.min(),checkSum.max())
                                print(" self.sand[ y, x]", self.sand[ y, x])
                                print(" self.silt[ y, x]", self.silt[ y, x])
                                print(" self.clay[ y, x]", self.clay[ y, x])
                                quit()
                                
######################
## if the dominant soil is a! ferrasol calculate mean soil properties (assuming values for all ferrasols):
                        if int(self.newGrid[ y, x]) == ferrasol :
                            ind2 = soilGroup == np.int(ferrasol)
                            fsbins = sbins[ind2]
                            fsvarfreq = svarfreq[ind2]
                            fscountSoilGroup = fsvarfreq.sum()
#####################

        fraction=0.042
        pad=0.06

        plt.subplot(221)
        plt.imshow( self.TcarbMassMod,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' A*** Full Mod calc TSoilCarb (0-0.3m)')
        plt.subplot(222)
        plt.imshow( self.ScarbMassMod,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' Full calc Mod SSoilCarb (0.3-1m) kg/m3')
        plt.subplot(223)
        plt.imshow( self.TcarbMass ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' Full calc TSoilCarb (0-0.3m) kg/m3')
        plt.subplot(224)
        plt.imshow( self.ScarbMass ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' Full calc SSoilCarb (0.3-1m) kg/m3')
        plt.savefig( "ferrasols1.png", orientation = "landscape")
        plt.clf()
        plt.close()

##
    def readLookUpGroups(self):
# Routine to find and label the major soil groups

        fid = open( "/project/jchmr/hadng/hwsd/HWSD_DATA.txt","r")  ## HWSD lookup table
        header = fid.readline()
        colMu = header.replace('"',"").split(",").index("MU_GLOBAL")        
        colsym74 = header.replace('"',"").split(",").index("SU_SYM74")
        colsym90 = header.replace('"',"").split(",").index("SU_SYM90")
        colcode90 = header.replace('"',"").split(",").index("SU_CODE90")

        colSa = header.replace('"',"").split(",").index("T_SAND")
        colSi = header.replace('"',"").split(",").index("T_SILT")
        colCl = header.replace('"',"").split(",").index("T_CLAY")
        colTOC = header.replace('"',"").split(",").index("T_OC")        
        colSOC = header.replace('"',"").split(",").index("S_OC")
        colCEC = header.replace('"',"").split(",").index("T_CEC_SOIL")
        colpH = header.replace('"',"").split(",").index("T_PH_H2O")
        colTBD = header.replace('"',"").split(",").index("T_BULK_DENSITY")
        colSBD = header.replace('"',"").split(",").index("S_BULK_DENSITY")
        colTREFBD = header.replace('"',"").split(",").index("T_REF_BULK_DENSITY")
        colSREFBD = header.replace('"',"").split(",").index("S_REF_BULK_DENSITY")
        colTGr = header.replace('"',"").split(",").index("T_GRAVEL")
        colSGr = header.replace('"',"").split(",").index("S_GRAVEL")

        print( "sym90 ", colsym90)
        lookup = fid.readlines()

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
        sym74=[ "G", "R", "I", "Q", "E", "U", "T", "V", "Z", "S", "Y" \
                   , "X", "K", "C", "H", "M", "B", "L", "D", "P", "W" \
                   , "A", "N", "F", "O", "J"]
# LP=lithic leptosols = very limited by depth over hard rock
        sym90=["GL","RG","LP","AR","-1","-2", "AN","VR","SC","SN","-3" \
                   ,"-4","KS","CH","PH","GR","CB","LV","PD","PZ", "PL"\
                   ,"AC","NT","FR","HS","FL"]

        print("len(lookup)",len(lookup))
        fao90_list=[]
        allMu=[]
        soilCode = {} #setting up dictionary
        soilInd = {}
        textures = {}
        oldmu=''
        for line in lookup:
# line.replace= replaces first string with second.
# This is effectively made into an array:

            line = line.replace('"',"").split(",") 
            fao90_line = line[colsym90]
            thismu=line[colMu]
            csum = int(line[colSa]) + int(line[colSi]) + int(line[colCl]) 
            if csum != 100:
                if csum != 0:
                    if csum != -29997:
                        print("check 100%",int(line[colSa]),  int(line[colSi]),  int(line[colCl]) )
                        print("check 100%",csum)
#            print("check 100%",[colSa] + [colSi] + [colCl] )
            if thismu != oldmu :  # only consider first find of soil label
                allMu.append(int(thismu))

#OC is in % kgC kg^-1
#TBD is in kg dm^-3
                textures[line[colMu]] = ( line[colSa], line[colSi], \
                                              line[colCl], line[colTOC], \
                                              line[colCEC],line[colpH], \
                                              line[colTBD], line[colSOC], \
                                              line[colSBD],line[colTREFBD],\
                                              line[colSREFBD], line[colTGr], line[colSGr] )
#                print("checking soil carbon",line[colTOC],line[colSOC])

## There is a bug in the text file at number 20941 (from the FAO90 code
## and comparison with others, it should be an RG$ soil symbol.
                if fao90_line == "1.45" :
                    fao90_line = "RGe"

## if there is no code in FAO90, take appropriate code in FAO74
## and do conversion
                if fao90_line == "-9999" :
                    fao74_line=line[colsym74]
## just take first letter as this is all that's used in FAO74 soil groups:
                    fao74_line=fao74_line[0]  
                    ind=sym74.index(fao74_line)
                    fao90_line=sym90[ind]

# WARNING even after the above procedure some soils have been re-classed using FAO90 therefore have numbers < 0.

# ensure only have first two characters:
                if len(fao90_line) > 2 :
                    fao90_line = fao90_line[0] + fao90_line[1]

                if fao90_line == 'HS' :
                    print("checking soil carbon in histosols",line[colTOC],line[colSOC])
                old = fao90_line in fao90_list
                if old == False :
                    fao90_list.append(fao90_line)
                    print(fao90_list.index(fao90_line),fao90_line)
                    if fao90_line == 'FR' :
                        ferrasol = fao90_list.index(fao90_line)
                        print("ferrasol",ferrasol,fao90_line)
                        self.ferrasol = ferrasol
                    if fao90_line == 'HS' :
                        histosol = fao90_list.index(fao90_line)
                        self.histosol = histosol
#                        print("checking soil carbon in histosols",line[colTOC],line[colSOC])
                    if fao90_line == 'AN' :
                        andosol = fao90_list.index(fao90_line)
                        self.andosol = andosol
                soilCode[int(line[colMu])] = fao90_line
                soilInd[int(line[colMu])] = ( fao90_list.index(fao90_line) )
                oldmu = line[colMu]
#            if int(line[colMu]) < 7000:
#                print("textures[line[colMu]]",line[colMu],textures[line[colMu]][0])
        fid.close()
        self.Textures = textures
        self.soilInd =soilInd
        self.soilCode = soilCode
        self.fao90_index = fao90_list.index
        self.colMu = allMu
#        quit()



    def readLandMask(self):
        fn = \
           "WFD-EI-LandFraction2d.nc"
        fid = Dataset( fn, "r")
        self.mask = np.flipud(fid.variables["lsmask"][:])
        fid.close()
##
##
    def fillHoles(self,sa,si,cl):
#        dxy = 5
        dxy = self.dxy

        plt.figure()
        plt.subplot(211)
        plt.imshow(sa   )
        plt.title('sand a')

        for it in xrange(1000):

            ii = ( self.mask > 0) & ( sa <= 0) & ( si <= 0) & ( cl <= 0) 
            if ii.sum() < 1:
                print( "Filled all the holes!",it)
                break
            print( it, dxy, ii.sum())
            for y in xrange( self.nyOut):
                for x in xrange( self.nxOut):
                    if ii[y,x] > 0:
                        ind =sa[y-dxy:y+dxy,x-dxy:x+dxy] \
                            + si[y-dxy:y+dxy,x-dxy:x+dxy] \
                            + cl[y-dxy:y+dxy,x-dxy:x+dxy] > 0
                        vsa = self.sand[y-dxy:y+dxy,x-dxy:x+dxy][ind]
                        vsi = self.silt[y-dxy:y+dxy,x-dxy:x+dxy][ind]
                        vcl = self.clay[y-dxy:y+dxy,x-dxy:x+dxy][ind]

                        if len(ind) > 0 :
                            ind2 = (vsa >= 0) & (vsi >= 0) & (vcl >= 0)
                            if len(ind2) > 0 :
                                sa[y,x] = np.mean(vsa[ind2])
                                si[y,x] = np.mean(vsi[ind2])
                                cl[y,x] = np.mean(vcl[ind2])
                                checkS=sa[y,x]+si[y,x]+cl[y,x]

#            dxy += 10
            dxy += 2*self.dxy

        return sa, si, cl
##
##
##
    def fillHoles_1fld(self,fld):
#        dxy = 5
        dxy = self.dxy

        temp_fld = deepcopy(fld)
        for it in xrange(1000):

            ii = ( self.mask > 0) & ( fld < 0) 
            if ii.sum() < 1:
                print( "Filled all the holes for fld!",it)
                break
            print( it, dxy, ii.sum())
            for y in xrange( self.nyOut):
                for x in xrange( self.nxOut):
                    if ii[y,x] > 0:
                        ind =temp_fld[y-dxy:y+dxy,x-dxy:x+dxy] >= 0
                        if ind.sum() > 0 :
                            vfld = temp_fld[y-dxy:y+dxy,x-dxy:x+dxy][ind]
                            fld[y,x] = np.mean(vfld)

#            dxy += 10
            dxy += 2*self.dxy


        return fld
##
##

    def doStuff(self):
##
        fraction=0.042
        pad=0.06
        

        self.readLookUpGroups()
### Read Emma's lookup table

        dy = self.nyIn / self.nyOut
        dx = self.nxIn / self.nxOut

## read in the raw 16 bit integer file with the soil codes
        self.readBIL()
        self.regridHSWD(self.nxIn / self.nxOut, self.nyIn / self.nyOut)



##
## Now we have to check this against the land sea mask we are using to see if
##   there are any points with no paramaters.        
        self.readLandMask()

# soil groups:
        ferra = deepcopy(self.newGrid)
        histo = deepcopy(self.newGrid)
        ando = deepcopy(self.newGrid)
        ferra = np.ma.masked_array( self.newGrid, self.newGrid != np.int(self.ferrasol))
        histo = np.ma.masked_array( self.newGrid, self.newGrid != np.int(self.histosol))
        ando = np.ma.masked_array( self.newGrid, self.newGrid != np.int(self.andosol))

################################################################
# Fill gaps:

        
        if self.fill:
            self.TcarbPerc= self.fillHoles_1fld(self.TcarbPerc) 
            self.ScarbPerc= self.fillHoles_1fld(self.ScarbPerc) 

        self.TcarbPerc = np.ma.masked_array( self.TcarbPerc, self.TcarbPerc < 0.0)
        self.ScarbPerc = np.ma.masked_array( self.ScarbPerc, self.ScarbPerc < 0.0)
        self.TcarbPerc = np.ma.masked_array( self.TcarbPerc, self.mask <= 0.0)
        self.ScarbPerc = np.ma.masked_array( self.ScarbPerc, self.mask <= 0.0)


        if self.fill:
            self.TcarbMass= self.fillHoles_1fld(self.TcarbMass) 
            self.ScarbMass= self.fillHoles_1fld(self.ScarbMass) 
            self.TcarbMassRef= self.fillHoles_1fld(self.TcarbMassRef) 
            self.ScarbMassRef= self.fillHoles_1fld(self.ScarbMassRef) 
            self.TcarbMassMod= self.fillHoles_1fld(self.TcarbMassMod) 
            self.ScarbMassMod= self.fillHoles_1fld(self.ScarbMassMod) 

        self.TcarbMass = np.ma.masked_array( self.TcarbMass, self.TcarbMass < 0.0)
        self.ScarbMass = np.ma.masked_array( self.ScarbMass, self.ScarbMass < 0.0)
        self.TcarbMassRef = np.ma.masked_array( self.TcarbMassRef, self.TcarbMassRef < 0.0)
        self.ScarbMassRef = np.ma.masked_array( self.ScarbMassRef, self.ScarbMassRef < 0.0)
        self.TcarbMassMod = np.ma.masked_array( self.TcarbMassMod, self.TcarbMassMod < 0.0)
        self.ScarbMassMod = np.ma.masked_array( self.ScarbMassMod, self.ScarbMassMod < 0.0)

        self.TcarbMass = np.ma.masked_array( self.TcarbMass, self.mask <= 0.0)
        self.ScarbMass = np.ma.masked_array( self.ScarbMass, self.mask <= 0.0)
        self.TcarbMassRef = np.ma.masked_array( self.TcarbMassRef, self.mask <= 0.0)
        self.ScarbMassRef = np.ma.masked_array( self.ScarbMassRef, self.mask <= 0.0)
        self.TcarbMassMod = np.ma.masked_array( self.TcarbMassMod, self.mask <= 0.0)
        self.ScarbMassMod = np.ma.masked_array( self.ScarbMassMod, self.mask <= 0.0)



        self.carbMass = 0.3*self.TcarbMass + 0.7*self.ScarbMass
        self.carbMassRef = 0.3*self.TcarbMassRef + 0.7*self.ScarbMassRef
        self.carbMassMod = 0.3*self.TcarbMassMod + 0.7*self.ScarbMassMod
        self.carbMass = np.ma.masked_array( self.carbMass, self.carbMass < 0.0)
        self.carbMassRef = np.ma.masked_array( self.carbMassRef, self.carbMassRef < 0.0)
        self.carbMassMod = np.ma.masked_array( self.carbMassMod, self.carbMassMod < 0.0)
        self.carbMass = np.ma.masked_array( self.carbMass, self.mask <= 0.0)
        self.carbMassRef = np.ma.masked_array( self.carbMassRef, self.mask <= 0.0)
        self.carbMassMod = np.ma.masked_array( self.carbMassMod, self.mask <= 0.0)

#        self.carbMass= self.fillHoles_1fld(self.carbMass) 
#        self.carbMassRef= self.fillHoles_1fld(self.carbMassRef) 
#        self.carbMassMod= self.fillHoles_1fld(self.carbMassMod) 


################################################################
# Now calculate 0->1m mean bulk densities and soil carbon
        RefBD = 0.3*self.TRefBD+0.7*self.SRefBD # ref bulk density
        BD = 0.3*self.TBD+0.7*self.SBD    # bulk density
        

        carbPerc = 0.3*self.TcarbPerc+0.7*self.ScarbPerc

#        nlat=len(self.sand)
#        nlon=len(self.sand[0])
        nlat=self.nyOut
        nlon=self.nxOut
        lat=np.zeros(nlat)
        lon=np.zeros(nlon)
        for j in xrange(0,nlat):
            lat[j] = -90.0+(j+0.5)*0.5
        for i in xrange(0,nlon):
            lon[i] = -180.0+(i+0.5)*0.5
 

        RefBD = np.ma.masked_array( RefBD, self.mask <= 0.0)
        BD = np.ma.masked_array( BD, self.mask <= 0.0)

        self.TRefBD = np.ma.masked_array( self.TRefBD, self.mask <= 0.0)
        self.SRefBD = np.ma.masked_array( self.SRefBD, self.mask <= 0.0)
        self.TBD = np.ma.masked_array( self.TBD, self.mask <= 0.0)
        self.SBD = np.ma.masked_array( self.SBD, self.mask <= 0.0)


###
        sgroups = self.newGrid
        sand=self.sand
        silt=self.silt
        clay=self.clay
        TcarbPerc=self.TcarbPerc
        ScarbPerc=self.ScarbPerc
        TBD=self.TBD
        SBD=self.SBD
        TRefBD=self.TRefBD
        SRefBD=self.SRefBD

# Check soil totals add to 100%:
        checkSum = sand + silt + clay
        checkSum = np.ma.masked_array( checkSum, sand < 0.0)       
        vv =np.logical_and(checkSum > 0, checkSum < 99.9) 
        print("min max",checkSum.min(),checkSum.max())
        print("sand,silt,clay",sand[vv],silt[vv],clay[vv])


        sgroups = np.ma.masked_array( sgroups, sgroups <= 0.0)
        sand = np.ma.masked_array( sand, sand < 0.0)
        silt = np.ma.masked_array( silt, silt < 0.0)
        clay = np.ma.masked_array( clay, clay < 0.0)
        

        TcarbPerc = np.ma.masked_array( TcarbPerc, self.mask <= 0.0)
        ScarbPerc = np.ma.masked_array( ScarbPerc, self.mask <= 0.0)
        TBD = np.ma.masked_array( TBD, TBD < 0.0)
        SBD = np.ma.masked_array( SBD, SBD < 0.0)
        TRefBD = np.ma.masked_array( TRefBD, TRefBD < 0.0)
        SRefBD = np.ma.masked_array( SRefBD, SRefBD < 0.0)

        
        plt.subplot(211)
        plt.imshow( self.TcarbMassMod,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('TSoilCarbMod (0-0.3m)')
        plt.subplot(212)
        plt.imshow( self.ScarbMassMod,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('SSoilCarbMod (0.3-1m) kg/m3')
        plt.savefig( "ferrasols2.png", orientation = "landscape")
        plt.clf()
        plt.close()

        plt.subplot(211)
        plt.imshow( self.TcarbMass ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' TSoilCarb (0-0.3m) kg/m3')
        plt.subplot(212)
        plt.imshow( self.ScarbMass ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' SSoilCarb (0.3-1m) kg/m3')
        plt.savefig( "ferrasols3.png", orientation = "landscape")
        plt.clf()
        plt.close()

        plt.subplot(211)
        plt.imshow( self.TcarbMassRef ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' TSoilCarbRef (0-0.3m) kg/m3')
        plt.subplot(212)
        plt.imshow( self.ScarbMassRef ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' SSoilCarbRef (0.3-1m) kg/m3')
        plt.savefig( "ferrasols4.png", orientation = "landscape")
        plt.clf()
        plt.close()

        plt.subplot(221)
        plt.imshow( self.TcarbPerc ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' TSoilCarbPerc (0-0.3m) %')
        plt.subplot(222)
        plt.imshow( self.ScarbPerc ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' SSoilCarbPerc (0.3-1m) %')
        plt.subplot(223)
        plt.imshow( carbPerc ,vmin=0,vmax=40)            
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' SoilCarbPerc (0-1m) %')
        plt.savefig( "ferrasols5.png", orientation = "landscape")
        plt.clf()
        plt.close()


        plt.figure()
        plt.subplot(221)
        plt.imshow( self.TBD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('T Bulk Density')
        plt.subplot(222)
        plt.imshow( self.SBD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('S Bulk Density')
        plt.subplot(223)
        plt.imshow( BD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title(' Bulk Density')
        plt.savefig( "ferrasols9.png", orientation = "landscape")
        plt.clf
        plt.close()


        plt.figure()
        plt.subplot(221)
        plt.imshow( self.TRefBD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('Ref T Bulk Density')
        plt.subplot(222)
        plt.imshow( self.SRefBD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('Ref S Bulk Density')
        plt.subplot(223)
        plt.imshow( RefBD, vmin=0.2, vmax=1.8)
        plt.colorbar(fraction=fraction, pad=pad)
        plt.title('Ref  Bulk Density')
        plt.savefig( "ferrasols10.png", orientation = "landscape")
        plt.clf
        plt.close()


        fillval=-999

        fid = Dataset( 'soil_cs.nc', "w", \
                       format = "NETCDF3_CLASSIC")
        fid.description = 'soil Carbon fields'
        fid.history = 'Created from ~hadng/WFDEI/make_ancils/soil_class_luts/HWSD2WFDEI.py'
        fid.source = 'HWSD file: /project/jchmr/hadng/hwsd/hwsd_rjel.bil'
        fid.createDimension( 'lon', nlon)
        fid.createDimension( 'lat', nlat)

        outvar=fid.createVariable('lon', np.float, ('lon') )
        outvar[:]=lon
        outvar=fid.createVariable('lat', np.float, ('lat') )
        outvar[:]=lat

        outvar=fid.createVariable('soilcarbon_Perc', np.float, ('lat','lon'),fill_value=fillval )
        outvar.units='% kg/kg (0->1m avg)'
        outvar[:]=carbPerc[::-1,:]  # need to reverse lat dimension



        outvar=fid.createVariable('TcarbMass', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->30cm avg) (uses Bulk Density)'
        outvar[:]=self.TcarbMass[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('ScarbMass', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (30->100cm avg) (uses Bulk Density)'
        outvar[:]=self.ScarbMass[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('carbMass', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->100cm avg) (uses Bulk Density)'
        outvar[:]=self.carbMass[::-1,:]  # need to reverse lat dimension

        outvar=fid.createVariable('TcarbMassRef', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->30cm avg) (uses Ref Bulk Density)'
        outvar[:]=self.TcarbMassRef[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('ScarbMassRef', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (30->100cm avg) (uses Ref Bulk Density)'
        outvar[:]=self.ScarbMassRef[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('carbMassRef', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->100cm avg) (uses Ref Bulk Density)'
        outvar[:]=self.carbMassRef[::-1,:]  # need to reverse lat dimension

        outvar=fid.createVariable('TcarbMassMod', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->30cm avg) (uses Bulk Density for Histosols and Andosols, otherwise Ref Bulk Density)'
        outvar[:]=self.TcarbMassMod[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('ScarbMassMod', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (30->100cm avg) (uses Bulk Density for Histosols and Andosols, otherwise Ref Bulk Density)'
        outvar[:]=self.ScarbMassMod[::-1,:]  # need to reverse lat dimension
        outvar=fid.createVariable('carbMassMod', np.float, ('lat','lon') ,fill_value=fillval )
        outvar.units='kgC/m3 (0->100cm avg) (uses Bulk Density for Histosols and Andosols, otherwise Ref Bulk Density)'
        outvar[:]=self.carbMassMod[::-1,:]  # need to reverse lat dimension

        fid.close()

if __name__ == "__main__":
    t=test()
    t.doStuff()

##
##
