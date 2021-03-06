#
# Python module of plotting routines
#
# Edward Comyn-Platt
# Centre for Ecology and Hydrology
# 2015
#
#
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.style.use('classic')
#
def custom_div_cmap(numcolors=6, name='custom_div_cmap',\
        colors=['lime','Green','Orange', 'Aqua', 'white', 'k']):
    cmap = LinearSegmentedColormap.from_list(name=name,
            colors=colors, N=numcolors)
    return cmap


import matplotlib.colors as mcolors
def custom_colour_list(
    numcolors=6, name='custom_div_cmap',
    colors=['lime','Green','Orange', 'Aqua', 'white', 'k']
):
    cmap = mcolors.LinearSegmentedColormap.from_list(name=name,
            colors=colors, N=numcolors)
    clist = [mcolors.rgb2hex(cmap(i)) for i in range(cmap.N)] 
    return clist

def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []   
    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []
        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
            reverse.append(sorted(data))    
    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r

#def reverse_colourmap(cmap):
#    reverse = []
#    for channel in cmap:
#        data = []
#        for t in channel:
#            data.append((1 - t[0], t[1], t[2]))
#        reverse.append(sorted(data))
#        
#    return reverse

###############################################################################
# Function: plot_map
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Plot geographic data as a map plot with colorbar.
#
# Required Inputs: 
#    DATA - numpy array or list of data to be plotted, 
#             can be of any dimensions
#    LONS - numpy array or list of longitudes of DATA,
#             dimensions must be the same as DATA
#    LATS - numpy array or list of latitudes of DATA,
#             dimensions must be the same as DATA
# Optional Inputs:
#    DATA2 - numpy array or list of data to be plotted on top of DATA, 
#             can be of any dimensions
#    LONS2 - numpy array or list of longitudes of DATA2,
#             dimensions must be the same as DATA2
#    LATS2 - numpy array or list of latitudes of DATA2,
#             dimensions must be the same as DATA2
#    DATA_RANGE - Range of color bar scale [Min,Max].
#                  Data outside this range is set to Min/Max,
#                  unless SET_UNDER/SET_OVER are set.
#    LON_RANGE - longitude range to plot [West_Lon,East_Lon]
#    LAT_RANGE - latitude range to plot [South_Lat,North_Lat]
#    MAP_TYPE  - Method of plotting data, corresponding to MPL plotting techniques.
#                 Currently accepted 'Mesh' (Default), 'Contour' and 'Map'
#    COLOURS   - List of colours to use for plotting data.
#                 Colours must be listed in min to max order.
#    MPL_CBAR  - MPL colorbar to use for plotting
#    CMAP      - Colormap/Palette to use for plotting
#    NLEVELS   - Number of levels for plotting, Default sets to, in priority order:
#                   Number of CLEVELS, Number of COLOURS, or 5
#    CLEVELS   - List of intervals to split the colorbar into.
#    CBAR_ORIENTATION - 'horizontal' at the bottom or 'vertical' to the right
#    WIDTH/HEIGHT - in cm of the polotting window
#    CBAR_LABEL   - Label for Colorbar
#    INTERPOLATE_COLOURS - logical to interpolate the list  of colours (COLOURS) to match NLEVELS
#    TICK_FORMAT  - Format for the CBAR tick labels, default is '%0.3f' (3 decimal places)
#    CBAR_TICK_LENGTH=10 - length of colorbar ticks, to manually adjust the seperator tick
#    PLOT_TITLE   - Plot title, at top of plot
#    FONTSIZES    - List of FONTSIZES - [ plot_default, lat/lon labels, cbar, plot_title]
#    iDISPLAY     - Display? 'Y' or 'N'. default is 'N', save as file.
#    iCLOSE       - Close? 'Y' or 'N'. default is 'N', save as file.
#    FILE_PLOT    - Output filename
#    SET_OVER/SET_UNDER  - Colour for over/under bounds values
#    LATDEL/LONDEL - Lat/Lon deliminater for gridlines
#    RIVERS  - Logical to set Whether or not to overplot rivers 
#    RESOLUTION - Resolution of coastlines/countries/rivers.
#                   'h' high, 'i' intermdiate (default), 'l' low
#    PROJECTION - Projection of plot, corresponding to MPL.Basemap projections
#    LEFT_FRAC/RIGHT_FRAC/BOTTOM_FRAC/TOP_FRAC - fractions to cut out of npstere/spstere projections
#    BOUNDINGLAT - Bounding latitude for npstere/spstere projections
#    LON_0  - Central longitude for required projections
#    LAT_0  - Central latitude for required projections
#
def plot_map(DATA,LONS,LATS, 
             DATA2=None, LONS2=None, LATS2=None,
             GREYMASK=None, MASKCOLOR='grey', 
             DATA_RANGE=None, LON_RANGE=None, LAT_RANGE=None, 
             MAP_TYPE='Mesh', 
             COLOURS=None,MPL_CBAR=None,CMAP=None,NLEVELS=None,CLEVELS=None,
             TickLEVELS=None,NTICKS=None,TickLABELS=None, 
             CBAR_ORIENTATION='horizontal',CBAR_SIZE='6%',CBAR_PAD=0.3,INTERPOLATE_COLOURS=False,
             TICK_FORMAT='%.2f', CBAR_TICK_LENGTH=10, 
             WIDTH=12, HEIGHT=8, CBAR_LABEL=None,PLOT_TITLE=None, FONTSIZES=[10,10,12,12],
             iDISPLAY='N', iCLOSE='N', FILE_PLOT=None, iORIENTATION='landscape', 
             FIGURE=None, AXIS=None, 
             SET_OVER=None,SET_UNDER=None,LATDEL=None,LONDEL=None,
             RIVERS=False, COASTLINES=True, COUNTRIES=False, 
             RESOLUTION='i',PROJECTION='cyl', 
             LEFT_FRAC=0.,RIGHT_FRAC=0.,BOTTOM_FRAC=0.,TOP_FRAC=0., 
             BOUNDINGLAT=35,LON_0=-32,LAT_0=90.,RSPHERE=[6378137.00, 6356752.3142], 
             ):  #Fraction of polar-steroegraphic to use
    
    from matplotlib.pylab import rc, cm
    from mpl_toolkits.basemap import Basemap
    import matplotlib.colors as col
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #
    # If plotting ranges not set set to limits of data 
    if (DATA_RANGE==None):
        DATA_RANGE=[np.amin(DATA),np.amax(DATA)]
    # If plotting geo limits not set, set to integer limits of lats/lons
    if (LON_RANGE==None):
        if (PROJECTION in ['stere']):
            LON_RANGE= [ LONS[0,0], LONS[-1,-1] ]
        else:
            LON_RANGE=[np.floor(np.amin(LONS)),np.ceil(np.amax(LONS))]
    if (LAT_RANGE==None):
        if (PROJECTION in ['stere']):
            LAT_RANGE= [ LATS[0,0], LATS[-1,-1] ]
        else:
            LAT_RANGE=[np.floor(np.amin(LATS)),np.ceil(np.amax(LATS))]
    #
    # If NLEVELS not set, set to number of CLEVELS if set, else to number of COLOURS if set, else to 5
    if (NLEVELS==None):
        if (CLEVELS!=None):
            NLEVELS=len(CLEVELS)
        elif (COLOURS!=None):
            NLEVELS=len(COLOURS)
        else:
            NLEVELS=5
    # If NLEVELS set but so too is CLEVELS over-ride NLEVELS
    elif (CLEVELS!=None):
        print('plot_map: Number of CLEVELS superceding NLEVELS')
        NLEVELS=len(CLEVELS)
    #
    # If CLEVELS is not set construct CLEVELS as equidistant levels from min(data) to max(data) for NLEVELS
    if (CLEVELS==None):      
        CLEVELS=np.linspace(DATA_RANGE[0],DATA_RANGE[1],num=NLEVELS)   
    #
    # If TickLEVELS not set, set to CLEVELS
    if (TickLEVELS==None):
        if (NTICKS==None):
            TickLEVELS=CLEVELS
        else:
            # NTICKS defined, create Ticklevels based on interpolation of
            # DATA_RANGE and NTICKS
            TickLEVELS = np.linspace(DATA_RANGE[0],DATA_RANGE[1],num=NTICKS)
    #
    # If no Figure and/or Axis provided create provided create figure and/or axis
    if FIGURE==None and AXIS==None:
        FIGURE = plt.figure(figsize=(WIDTH,HEIGHT))
        AXIS   = FIGURE.add_subplot(1,1,1)
    elif AXIS==None:
        AXIS   = FIGURE.add_subplot(1,1,1)
    
    
    #
    FONT     = {'family' : 'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : FONTSIZES[0] }
    rc('font', **FONT) 
    rc('text',usetex=True)
    #
    # Create cylindrical Basemap instance.
    if (PROJECTION in ['npstere','spstere']):
        M         = Basemap(projection=PROJECTION,resolution=RESOLUTION, \
                            boundinglat=BOUNDINGLAT,lon_0=LON_0, \
                            ax=AXIS) 
        M.llcrnrx = M.urcrnrx*LEFT_FRAC
        M.urcrnrx = M.urcrnrx*(1.-RIGHT_FRAC)
        M.llcrnry = M.urcrnry*BOTTOM_FRAC
        M.urcrnry = M.urcrnry*(1.-TOP_FRAC)
        xx,yy     = M(LONS,LATS)
        if (DATA2!=None):
            xx2,yy2 = M(LONS2,LATS2)
    elif (PROJECTION.lower() in ['stere','emep','emep4uk']):
        if (PROJECTION.lower() in ['emep','emep4uk']):
            LON_0=-32.
            LAT_0=90.

        M         = Basemap(projection='stere',resolution=RESOLUTION, \
                            lon_0=LON_0, lat_0=LAT_0, \
                            llcrnrlat =LAT_RANGE[0],urcrnrlat=LAT_RANGE[1], \
                            llcrnrlon =LON_RANGE[0],urcrnrlon=LON_RANGE[1], \
                            ax=AXIS ) 
        xx,yy     = M(LONS,LATS)
        if (DATA2!=None):
            xx2,yy2 = M(LONS2,LATS2)
    elif (PROJECTION.lower() in ['chess','tmerc']):
        if (PROJECTION.lower() in ['chess']):
            LON_0=-2.
            LAT_0=49.
            RSPHERE=[6378137.00, 6356752.3142]

        M         = Basemap(projection='tmerc',resolution=RESOLUTION, \
                            lon_0=LON_0, lat_0=LAT_0,rsphere=RSPHERE, \
                            llcrnrlat =LAT_RANGE[0],urcrnrlat=LAT_RANGE[1], \
                            llcrnrlon =LON_RANGE[0],urcrnrlon=LON_RANGE[1], \
                            ax=AXIS ) 
        xx,yy     = M(LONS,LATS)
        if (DATA2!=None):
            xx2,yy2 = M(LONS2,LATS2)
    else:
        M         = Basemap(projection=PROJECTION, \
                            llcrnrlat =LAT_RANGE[0],urcrnrlat=LAT_RANGE[1], \
                            llcrnrlon =LON_RANGE[0],urcrnrlon=LON_RANGE[1], \
                            resolution=RESOLUTION,ax=AXIS )
        xx,yy     = M(LONS,LATS)
        if (DATA2!=None):
            xx2,yy2 = (LONS2,LATS2)
    #
    #,fix_aspect=ASPECT)
    ##AX        = plt.gca()            # Extract axis from 
    #
    # Set up a colormap:
    #
    if CMAP!=None:
        PALETTE=CMAP
        if (NLEVELS<255):
            NORM      = col.BoundaryNorm(CLEVELS, NLEVELS, clip=False)
    elif (MPL_CBAR==None):
        if (COLOURS==None):
            # Default set of colours to use if no colour bar is selected
            COLOURS   = [ '#00a0ff', '#00c8c8', '#00d28c' , '#a0e632', '#e6dc32', '#f08228' , '#fa3c3c' ]

        if not INTERPOLATE_COLOURS:
            NCOLOURS=len(COLOURS)
            COL_MAP   = col.ListedColormap(COLOURS[0:NCOLOURS],'indexed')
            cm.register_cmap(cmap=COL_MAP)
            PALETTE   = COL_MAP
            if (NCOLOURS<255):
                NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
        else:
            NCOLOURS=NLEVELS
            PALETTE=custom_div_cmap(numcolors=NCOLOURS, colors=COLOURS)
            if (NCOLOURS<255):
                NORM      = col.BoundaryNorm(CLEVELS, NCOLOURS, clip=False)
    else:
        PALETTE   = cm.get_cmap(name=MPL_CBAR,lut=NLEVELS)
        if (NLEVELS<255):
            NORM      = col.BoundaryNorm(CLEVELS, NLEVELS, clip=False)
                

    if SET_UNDER:
        PALETTE.set_under(SET_UNDER, 1.0)
    if SET_OVER:
        PALETTE.set_over(SET_OVER,  1.0)
    
    if (MAP_TYPE=='Map'):
        IMAGE    = M.imshow(DATA,cmap=PALETTE,interpolation='nearest',norm=NORM)
    elif (MAP_TYPE=='Mesh'):
        IMAGE    = M.pcolormesh(xx,yy,DATA,cmap=PALETTE,norm=NORM)
    elif (MAP_TYPE=='Contour'):
        IMAGE    = M.contourf(xx,yy,DATA,CLEVELS,cmap=PALETTE,norm=NORM,extend='both')
    
    if (DATA2!=None):
        if (MAP_TYPE=='Map'):
            M.imshow(DATA2,cmap=PALETTE,interpolation='nearest',norm=NORM)
        elif (MAP_TYPE=='Mesh'):
            M.pcolormesh(xx2,yy2,np.zeros_like(xx2),cmap='binary',norm=NORM)
            M.pcolormesh(xx2,yy2,DATA2,cmap=PALETTE,norm=NORM)
        elif (MAP_TYPE=='Contour'):
            M.contourf(xx2,yy2,np.zeros_like(xx2),cmap='binary',norm=NORM)
            M.contourf(xx2,yy2,DATA2,CLEVELS,cmap=PALETTE,norm=NORM)
    #
    
    #if (GREYMASK!=None):
    if GREYMASK is not None:
        mask_cmap=custom_div_cmap(numcolors=2,colors=[MASKCOLOR,MASKCOLOR])
        if (MAP_TYPE=='Mesh'):
            M.pcolormesh(xx,yy,GREYMASK,cmap=mask_cmap,norm=NORM)
        elif (MAP_TYPE=='Contour'):
            M.contourf(xx,yy,GREYMASK,cmap=mask_cmap,norm=NORM)

    
    if (PLOT_TITLE!=None):
        AXIS.set_title(PLOT_TITLE,fontsize = FONTSIZES[3])
    #
    # draw coastlines.
    if COUNTRIES:
        M.drawcountries(linewidth =0.4)
    if COASTLINES:
        M.drawcoastlines(linewidth =0.4)
    # draw rivers if desired
    if RIVERS:
        M.drawrivers(linewidth =0.1,color='grey')
    #
    # Draw a line around the map region.
    M.drawmapboundary()
    #
    # draw parallels and meridians.
    if (LATDEL):
        parallels = np.arange(-90.,90.,LATDEL)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5,fontsize=FONTSIZES[1])
    if (LONDEL):
        meridians = np.arange(-180.,180.,LONDEL)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5,fontsize=FONTSIZES[1])
    
    #
    # Colorbar options
    if (CBAR_ORIENTATION!='off'):
        DIVIDER   = make_axes_locatable(AXIS)
        # Horizontal at the bottom or Vertical to the right
        if CBAR_ORIENTATION=='horizontal':
            CB_AXIS   = DIVIDER.append_axes('bottom', size=CBAR_SIZE, pad=CBAR_PAD)
        elif CBAR_ORIENTATION=='vertical':
            CB_AXIS   = DIVIDER.append_axes('right', size=CBAR_SIZE, pad=CBAR_PAD)
        # pointed colour bar at one or both ends?
        # selection based on whether SET_UNDER/SET_OVER are set.
        # this method should change
        if SET_OVER:
            if SET_UNDER:
                COLORBAR  = plt.colorbar(IMAGE,cax=CB_AXIS,extend='both',orientation=CBAR_ORIENTATION,\
                                          ticks=TickLEVELS,format=TICK_FORMAT)
            else:
                COLORBAR  = plt.colorbar(IMAGE,cax=CB_AXIS,extend='max',orientation=CBAR_ORIENTATION,\
                                         ticks=TickLEVELS,format=TICK_FORMAT)
        elif SET_UNDER:
            COLORBAR  = plt.colorbar(IMAGE,cax=CB_AXIS,extend='min',orientation=CBAR_ORIENTATION,\
                                     ticks=TickLEVELS,format=TICK_FORMAT)
        else:
            COLORBAR  = plt.colorbar(IMAGE,cax=CB_AXIS,orientation=CBAR_ORIENTATION,\
                                     ticks=TickLEVELS,format=TICK_FORMAT)
        #
        if (CBAR_LABEL!=None):
            COLORBAR.set_label(CBAR_LABEL,fontsize=FONTSIZES[2])
        #
        # Set ticks in line 
        if CBAR_ORIENTATION=='horizontal':
            #COLORBAR.ax.axhline(linewidth=2,color='black')
            #COLORBAR.ax.axhline(y=1,linewidth=2,color='black')
            COLORBAR.ax.xaxis.set_tick_params(length=CBAR_TICK_LENGTH,bottom=True,top=True,\
                                              width=1.5,labelsize=FONTSIZES[0])
            if TickLABELS!=None:
                COLORBAR.ax.xaxis.set_ticklabels(TickLABELS,fontsize=FONTSIZES[0])

        elif CBAR_ORIENTATION=='vertical':
            COLORBAR.ax.axvline(linewidth=2,color='black')
            COLORBAR.ax.axvline(x=1,linewidth=2,color='black')
            COLORBAR.ax.yaxis.set_tick_params(length=CBAR_TICK_LENGTH,bottom=True,top=True,\
                                              width=1.5,labelsize=FONTSIZES[0])
            if TickLABELS!=None:
                COLORBAR.ax.yaxis.set_ticklabels(TickLABELS,fontsize=FONTSIZES[0])
        #

    #
    if iDISPLAY =='Y':
        # Display onscreen.
        plt.show() # display onscreen.
    elif (FILE_PLOT!=None):
        # Write to file
        plt.savefig(FILE_PLOT, orientation=iORIENTATION, bbox_inches='tight')
    
    if iCLOSE=='Y':
        plt.close()
    #
    # Return to calling routine
    #
    return IMAGE

###############################################################################
# Function: plot_map_multi
# Author: Edward Comyn-Platt, Feb 2016
# Purpose: Plot multi-map. 
#
# Required Inputs: 
#    DATA - A list of numpy arrays to be plotted, 
#             can be of any dimensions
#    LONS - A list of numpy arrays of longitudes of DATA,
#             dimensions must be the same as DATA
#    LATS - A list of numpy arrays of latitudes of DATA,
#             dimensions must be the same as DATA
# Optional Inputs:
#    Ncols(=None) - Number of columns. Default is to use number of maps to be plotted
#    Nrows(=None) - Number of rows. Default is 1
#    FIGSIZE(=(20,10)) - Figure size.
#    FILEPLOT - Filename to save plot
#    lDISPLAY - Logical to display the plot or not.
#    lCLOSE   - Logical to close the plot or not.
#
def plot_map_multi(DATA,LONS,LATS, \
                    Ncols=None,Nrows=None,FIGSIZE=(20,10),FONTSIZE=50, \
                    FILEPLOT=None, lDISPLAY=False, lCLOSE=False, SUPTITLE=None, \
                    COMMON_CBAR=False, NTICKS_COM=11, PLOT_TITLES=None, \
                    pad=0.06,fraction=0.15, \
                    **kwargs ):
    
    # get number of maps to be plotted
    Nmaps=len(DATA)
    # calculate ncols and nrows
    if (Ncols==None)&(Nrows==None):
        Nrows=1
        Ncols=Nmaps
    elif (Ncols!=None)&(Nrows==None):
        Nrows=int(np.ceil(float(Nmaps)/float(Ncols)))
    elif (Ncols==None)&(Nrows!=None):
        Ncols=int(np.ceil(float(Nmaps)/float(Nrows)))
    
    if COMMON_CBAR:
        kwargs['CBAR_ORIENTATION']='off'

    # Create figure and axes:
    FIG,AXES = plt.subplots(ncols=Ncols,nrows=Nrows,figsize=FIGSIZE)
    # Loop round each map
    for iMAP,AX in zip(range(Nmaps),AXES.flat):
        IMAGE=plot_map(DATA[iMAP],LONS[iMAP],LATS[iMAP],AXIS=AX,\
                        PLOT_TITLE=PLOT_TITLES[iMAP], \
                     **kwargs )
    if COMMON_CBAR:
        if 'DATA_RANGE' in kwargs:
            DATA_RANGE=kwargs['DATA_RANGE']
            TickLEVELS = np.linspace(DATA_RANGE[0],DATA_RANGE[1],num=NTICKS_COM)
            CBAR=plt.colorbar(IMAGE,ax=AXES.flatten().tolist(),orientation='horizontal',\
                                    ticks=TickLEVELS,pad=pad,fraction=fraction)
        else:
            CBAR=plt.colorbar(IMAGE,ax=AXES.flatten().tolist(),orientation='horizontal',\
                               pad=pad,fraction=fraction)
        CBAR.ax.axhline(linewidth=2,color='black')
        CBAR.ax.axhline(y=1,linewidth=2,color='black')
        CBAR.ax.xaxis.set_tick_params(length=30,bottom=True,top=True,\
                                           width=1.5,labelsize=FONTSIZE*0.5)
        if 'CBAR_LABEL' in kwargs:
            CBAR.set_label(kwargs['CBAR_LABEL'],fontsize=FONTSIZE*0.7)
        
    if SUPTITLE!=None:
        FIG.suptitle(SUPTITLE,fontsize=FONTSIZE) #,y=0.92)
    
    if FILEPLOT!=None:
        FIG.savefig(FILEPLOT, bbox_inches='tight')

    if lDISPLAY:
        plt.show()
    
    if lCLOSE:
        plt.close()

#
#
###############################################################################
# Function: plot_timeseries
# Author: Edward Comyn-Platt, Feb 2015
# Purpose: Plot time-series data. Can plot multiple series to one plot
#
# Required Inputs: 
#    DATA - single numpy array or list or list of numpy arrays or lists of data to be plotted.
#    TIME - datetime object or list of datetime objects corresponding to each DATA.
#           one datetime object can be replicated for for each DATA
#    DATA_RANGE - [min,max] 2 element list of max and min for data (y) axis
#    TIME_RANGE - [min,max] 2 element list of max and min for time (x) axis
#    PLOT_TITLE - Title for plot
#    WIDTH/HEIGHT - in cm of the polotting window
#    FONTSIZES    - Font sizes for virous texts:
#                      [Axis labels, tick labels, Title, Legend]
#    COLOURS - Colour or list of Colours for plot lines and points,
#                     if List, should correspond to number of DATA series 
#                     unallocated colours will be final colour repeated
#                     accepts any matplotlib type colour
#    MARKERS - mpl marker string, or list of markers,
#                   if List, should correspond to number of DATA series 
#                   unallocated markers will be final marker repeated
#                   accepts any matplotlib type marker
#    LINESTYLE  - linestyle, or list of linestyles for plot lines,
#                     if List, should correspond to number of DATA series 
#                     unallocated linestyles will be final linestyle repeated
#                     accepts any matplotlib type linestyle
#    LEGEND - Set indicate location of LEGEND if desired.
#              Accepts any MPL legend location code or string
#              If not set, no legend is drawn
#    LEGEND_DATANAMES - Datanames for legend, list of strings corresponding to DATA series order
#    PLOT_AXIS    - Set to plot predfined plot axis if exist, 
#                    useful for incorporating this routine in a multiplot routine
#    iDISPLAY     - Display? 'Y' or 'N'. default is 'N', save as file.
#    FILE_PLOT    - Output filename, set to save plot within routine.
#    Y_LABEL      - label for y-axis 
#    
#    
#
def plot_timeseries(DATA,TIME, \
                        DATA_RANGE=None, TIME_RANGE=None, \
                        PLOT_TITLE=None, Y_LABEL = None, \
                        WIDTH=12, HEIGHT=6, FONTSIZES=[10,10,12,10], \
                        COLOURS=None, MARKERS=None, LINESTYLES=None, \
                        LEGEND=None, LEGEND_DATANAMES=None, \
                        PLOT_AXIS=None, iDISPLAY='N', FILE_PLOT=None):
    # 
    # import Function specific modulesImport neccessary modules
    import numpy as np
    import copy
    import matplotlib.colors as col
    import matplotlib.pyplot as plt
    from matplotlib import dates as mdates
    #
    # Find out format of DATA and how many time-series we are plotting
    # Convert to list of numpy arrays for mathematical convenience
    if (type(DATA).__name__=='list'):
        if (type(DATA[0]).__name__=='list'):
            #Conver list of lists to list of np.arrays
            DATA=[np.array(DAT) for DAT in DATA]
            Nseries = len(DATA)
        elif ('rray' in type(DATA[0]).__name__):
            # Good, no conversion required
            Nseries = len(DATA)
        elif (type(DATA[0]).__name__=='int')      | \
                (type(DATA[0]).__name__=='float') | \
                (type(DATA[0]).__name__=='long')  | \
                (type(DATA[0]).__name__=='complex'):
            # convert list to list of np.array
            DATA=[np.array(DATA)]
            Nseries = 1
    elif ('rray' in type(DATA).__name__):
        # convert np.array to list of np.array
        DATA=[DATA]
        Nseries=1
    else:
        print('ERROR in plot_timeseries: Unrecognised DATA format, exiting plotting procedure')
        return
    #
    # Find out format of TIME and how many time-series we are plotting
    # Convert to list of numpy arrays for mathematical convenience
    if (type(TIME).__name__=='list'):
        if (len(TIME)==1) & (Nseries>1):
            TIME=[TIME for i in range(Nseries)]
        elif (len(TIME)!=Nseries):
            print('ERROR in plot_timeseries: TIME not compatible with DATA')
            return   
    elif ('array' in type(TIME).__name__):
        TIME=[TIME for i in range(Nseries)]
    
    if not all([(len(DATA[i])==len(TIME[i])) for i in range(Nseries)]):
        print([(len(DATA[i])==len(TIME[i])) for i in range(Nseries)])
        print('ERROR in plot_timeseries: Lengths of time and data series do not match')
        return   
    # 
    # 
    # Check for line colours list, if doesn't exist set to list of red for each series
    if (COLOURS==None):
        COLOURS=['red' for i in range(Nseries)]
    # else, if a single string, convert to list
    elif (type(COLOURS)=='str'):
        COLOURS=[COLOURS for i in range(Nseries)]
    # else, if len of COLOURS list is less than Nseries, fill remaining colours with last colour
    elif (type(COLOURS)=='list'):
        # if list check if equal to or longer than Nseries
        if(len(COLOURS)<Nseries):
            for i in range(Nseries-len(COLOURS)):
                COLOURS.append(COLOURS[-1])
    
    # Check for LINESTYLES
    if (LINESTYLES==None):
        # if None, set to list of blank strings
        LINESTYLES=['' for i in range(Nseries)]
    elif (type(LINESTYLES)=='str'):
        # if single string convert to list
        LINESTYLES=[LINESTYLES for i in range(Nseries)]
    elif (type(LINESTYLES)=='list'):
        # if list check if equal to or longer than Nseries
        if (len(LINESTYLES)<Nseries):
            # fill unfilled with final linestyle
            for i in range(Nseries-len(LINESTYLES)):
                LINESTYLES.append(LINESTYLES[-1])
    
    # Check for MARKERS
    if (MARKERS==None):
        # if None, set to list of blank strings
        MARKERS=['' for i in range(Nseries)]
    elif (type(MARKERS)=='str') | (type(MARKERS)=='tuple'):
        # if single string or single tuple convert to list
        MARKERS=[MARKERS for i in range(Nseries)]
    elif (type(MARKERS)=='list'):
        # if list check if equal to or longer than Nseries
        if (len(MARKERS)<Nseries):
            # fill unfilled with final linestyle
            for i in range(Nseries-len(MARKERS)):
                MARKERS.append(MARKERS[-1])
    
    # If PLOT_AXIS not set create new axis as subplot:
    if (PLOT_AXIS==None):
        FIG,PLOT_AXIS = plt.subplots()
    else:
        # If Plotting to existing axis overide any display or save commands
        if (iDISPLAY=='Y'):
            print('WARNING in plot_timeseries: Plotting to existing AXIS, overwriting iDISPLAY command')
            iDISPLAY='N'
        if (FILE_PLOT!=None):
            print('WARNING in plot_timeseries: Plotting to existing AXIS, overwriting FILE_PLOT command')
            FILE_PLOT=None
    
    # Plot data series
    for DAT,TIM,COL,LS,MK in zip(DATA,TIME,COLOURS,LINESTYLES,MARKERS):  
        PLOT_AXIS.plot(TIM,DAT,color=COL,ls=LS,marker=MK)
    
    # If DATA_RANGE set, set y-axis accordingly
    if (DATA_RANGE!=None):
        PLOT_AXIS.axes.set_ylim(bottom=DATA_RANGE[0],top=DATA_RANGE[1])
        
    # If TIME_RANGE set, set x-axis accordingly
    if (TIME_RANGE!=None):
        PLOT_AXIS.axes.set_xlim(left=TIME_RANGE[0],right=TIME_RANGE[1])
    
    if (Y_LABEL!=None):
        PLOT_AXIS.set_ylabel(Y_LABEL,fontsize=FONTSIZES[1])

    if (PLOT_TITLE!=None):
        PLOT_AXIS.set_title(PLOT_TITLE,fontsize = FONTSIZES[3])
    
    if (LEGEND!=None):
        if (LEGEND_DATANAMES==None):
            LEGEND_DATANAMES=[str(i) for i in range(Nseries)]
        PLOT_AXIS.legend(LEGEND_DATANAMES,loc=LEGEND)
    
    FIG.set_size_inches(WIDTH,HEIGHT)
    
    # Display if desired
    if (iDISPLAY=='Y'):
        FIG.show()
    # Save to file if desired
    elif (FILE_PLOT!=None):
        FIG.savefig(FILE_PLOT,dpi=100, bbox_inches='tight')
    
    return
        
#

###################################################################################
# plot_taylor_diagram
#
# Produce a Taylor plot for one or more series of standard deviations 
#     and correlation values
#
# Input: STD - a numpy array or a list of numpy arrays containing standard deviations
#        CORR- a numpy array or a list of numpy arrays containing correlation
#
# Options: smax - maximum standard deviation, default is the max value in STD
#          corr_range - range of correlation values, default is [0,1], only other usable
#                         option is [-1,1]
#          FIG - Figure to plot to, if not defined a figure is created and returned
#          AX - Axis to plot to, if not defined Axis is created 
#          FIG_subplot - if creating axis choose location based on standard mpl syntax
#                           default is (111) = entire figure space
#          iDISPLAY - display the figure as part of function call, default is 'N'
#          COLORS - list colors to plot the series. If number of colors is less than the 
#                      number of series COLORS are repeated.
#          MARKERS - as COLORS but for plotting marker shape
#          SYMSIZE - as COLORS but for plotting marker size
#          LINEWIDTHS - as COLORS but for plotting marker thickness
#          FILE_PLOT - filename to save the figure to
#          PLOT_LEGEND - boolean to control whether or not to plot legend, default is False
#          LABELS - list of legend labels, if not set labels are just a sequence of numbers
#          FONTSIZES - [ AX_title, axis labels, FIG_title, unused ]
#          AX_TITLE - title for Axis, Default is None
#          FIG_TITLE - title for Figure, default is None
#
# 
####################################################################################
def plot_taylor_diagram( STD, CORR ,\
                         smax=None, corr_range=[0,1] , \
                         FIG=None, AX=None, FIG_subplot=(111), iDISPLAY='N', \
                         COLORS=['b'],MARKERS=['x'],SYMSIZES=[30],LINEWIDTHS=[1.5], \
                         FILE_PLOT=None, PLOT_LEGEND=False, LABELS=None, \
                         FONTSIZES=[20,15,25,15], \
                         AX_TITLE=None, FIG_TITLE=None, \
                         ):
    
    from matplotlib.projections import PolarAxes
    import mpl_toolkits.axisartist.floating_axes as FA
    import mpl_toolkits.axisartist.grid_finder as GF
    import types

    if type(STD)!=list:
        STD=[STD]
    if type(CORR)!=list:
        CORR=[CORR]
    
    nSERIES=len(STD)
    if len(STD)!=len(CORR):
        print('ERROR in plot_taylor_diagram: ')
        print('STD and CORR are different lengths')
   
    # Popoulate plotting option lists where neccessary
    if type(COLORS)!=list:
        COLORS=[COLORS]
    while len(COLORS)<len(STD):
        COLORS=COLORS+COLORS
    if type(MARKERS)!=list:
        MARKERS=[MARKERS]
    while len(MARKERS)<len(STD):
        MARKERS=MARKERS+MARKERS
    if type(SYMSIZES)!=list:
        SYMSIZES=[SYMSIZES]
    while len(SYMSIZES)<len(STD):
        SYMSIZES=SYMSIZES+SYMSIZES
    if type(LINEWIDTHS)!=list:
        LINEWIDTHS=[LINEWIDTHS]
    while len(LINEWIDTHS)<len(STD):
        LINEWIDTHS=LINEWIDTHS+LINEWIDTHS

    # If labels not defined create
    if (LABELS is None):
        LABELS=[str(i) for i in range(nSERIES)]

    tr = PolarAxes.PolarTransform()
    corr_range_rad = [ (1+corr_r)*(np.pi/2) for corr_r in corr_range ] 
    # Correlation labels
    if (corr_range[0]==0):
        rlocs = np.concatenate((np.arange(0.0,0.91,0.1),[0.95,0.99]))
    elif (corr_range[0]==-1):
        rlocs = np.concatenate(([-0.99,-0.95,-0.9],np.arange(-0.8,0.9,0.2),[0.9,0.95,0.99]))
        
    rlocs[np.absolute(rlocs)<0.01] = 0.0

    tlocs = np.arccos(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs)    # Positions
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str,rlocs))))
    
    # Standard deviation axis extent
    smin = 0
    if smax==None:
        smaxes = [ np.max(np.array(std)) for std in STD ]
        smax = int(np.ceil(max(smaxes)*10))/10.
    
    ghelper = FA.GridHelperCurveLinear(tr,
                                       extremes=(corr_range_rad[0],corr_range_rad[1],
                                                 smin,smax),
                                       grid_locator1=gl1,
                                       tick_formatter1=tf1,)
    
    #if FIG is None:
    if (FIG is None) & (AX is None):
        FIG = plt.figure()

    if (AX is None):
        AX = FA.FloatingSubplot(FIG,FIG_subplot, grid_helper=ghelper)
        FIG.add_subplot(AX)
        # Adjust axes
        # Correlation Axis (top and bottom)
        AX.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        AX.axis["top"].toggle(ticklabels=True, label=True)
        AX.axis["top"].major_ticklabels.set_axis_direction("top")
        AX.axis["top"].label.set_axis_direction("top")
        AX.axis["top"].label.set_text("Correlation") #,fontsize=FONTSIZES[1])
        AX.axis["bottom"].set_visible(False)         # Useless
    
        # Standard Deviation Axis - left and right
        AX.axis["left"].set_axis_direction("bottom") # "X axis"
        AX.axis["left"].label.set_text("Standard deviation") #,fontsize=FONTSIZES[1])
        AX.axis["right"].set_axis_direction("top")   # "Y axis"
    
        # Contours along standard deviations
        AX.grid(True)
        
    # Get Plotting Axis in polar coordinates
    AX_plt = AX.get_aux_axes(tr)

    i=0
    for std,corr in zip(STD,CORR):
        plt_std=np.array(std)
        plt_corr=np.arccos(np.array(corr))
        AX_plt.scatter(plt_corr,plt_std, \
                    s=SYMSIZES[i],marker=MARKERS[i],\
                    c=COLORS[i],lw=LINEWIDTHS[i],   \
                    label=LABELS[i] )
        i+=1
   
    if (AX_TITLE!=None):
        AX_plt.set_title(AX_title,fontsize=FONTSIZES[0])

    if (PLOT_LEGEND==True):
        handles,labels=AX_plt.get_legend_handles_labels()
        FIG.legend(handles,labels,pos='upper right')
    
    if (FIG_TITLE!=None):
        FIG_plt.set_title(FIG_title,fontsize=FONTSIZES[0])

    if FILE_PLOT!=None:
        FIG.savefig(FILE_PLOT,dpi=100, bbox_inches='tight')

    if iDISPLAY=='Y':
        plt.show()
    else:
        return FIG, AX
  
