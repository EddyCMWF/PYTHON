# Python module containing a number of plotting routine
#
# Garry Hayman
# Centre for Ecology and Hydrology
# July 2012
#
# Contains:
#
# Plot_Histogram_Single                     : single histogram
# Plot_Regression                           : linear regression
# Plot_Zonal                                : zonal plot (against latitude)
#
import math
import numpy as np
from numpy import arange
import matplotlib.pyplot as plt
from matplotlib import rc,cm,colorbar
#
# *****************************************************************************************
# Histograms
# *****************************************************************************************
#
def Plot_Histogram_Single(DATA,DATA_MIN,DATA_MAX,RANGE, \
	PLOT_TITLE,FACECOLOUR,XLABEL,YLABEL,XTICKS,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(6.0,6.0))
#
	NUMBER,BINS,PATCHES=plt.hist(DATA,bins=RANGE,normed=0,facecolor=FACECOLOUR,alpha=0.75)
#
	if DEBUG == 'Y':
		print(RANGE)
		print(DATA)
#
	MAX_COUNT      = NUMBER.max()
	POWER          = 0.1
#
	while MAX_COUNT > 1:
		POWER       = 10.0*POWER
		MAX_COUNT   = MAX_COUNT/POWER
#
	NYTICKS        = int(1+MAX_COUNTS.max()/POWER)
	MAX_COUNT      = POWER*NYTICKS
	YTICKS         = int(POWER)*arange(NYTICKS+1)
#
	TOTAL          = NUMBER.sum()
	MEAN           = DATA[DATA > 0.0].mean()
	STD            = np.std(DATA[DATA > 0.0])
	print(len(DATA),TOTAL,MEAN,DATA.max(),MAX_COUNT)
	plt.axis([DATA_MIN,DATA_MAX,0.0,MAX_COUNT])
#
	TEXT           = 'Counts = '+('%d' % TOTAL)+'; Mean = '+('%6.1f' % MEAN)+' ppb; std. dev = '+('%6.1f' % STD)+' ppb'
	plt.title(TEXT,fontsize=12)
	plt.suptitle(PLOT_TITLE,fontsize=8)
	plt.xlabel(XLABEL,fontsize=10)
	plt.ylabel(YLABEL,fontsize=10)
	plt.xticks(XTICKS,XTICKS)
	plt.yticks(YTICKS,YTICKS)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_Histogram_Multi(NDATA,DATA,DATA_MIN,DATA_MAX,RANGE,ORIENTATION, \
	PLOT_TITLE,FACECOLOUR,XLABEL,YLABEL,XTICKS,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 10 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
		SUB_PLOT0    = '1'+str(NDATA)
	else:
		plt.figure(figsize=(8.0,12.0))
		SUB_PLOT0    = str(NDATA)+'1'
#
# Need to cycle through each dataset to define a common maximum for the histogram plots
#
	MAX_COUNTS   = np.zeros(NDATA)
	DATA         = np.array(DATA)
#
	for iDATA in range(NDATA):
		NUMBER,BINS,PATCHES = plt.hist(DATA[iDATA,:],bins=RANGE,normed=0,facecolor=FACECOLOUR[iDATA],alpha=0.75)
		MAX_COUNTS[iDATA]   = NUMBER.max()
	plt.close()
#
	if DEBUG == 'Y':
		print(RANGE)
		print(DATA)
#
	MAX_COUNT      = MAX_COUNTS.max()
	POWER          = 0.1
#
	while MAX_COUNT > 1:
		POWER       = 10.0*POWER
		MAX_COUNT   = MAX_COUNT/POWER
#
	NYTICKS        = int(1+MAX_COUNTS.max()/POWER)
	MAX_COUNT      = POWER*NYTICKS
	YTICKS         = int(POWER)*arange(NYTICKS+1)
#
	for iDATA in range(NDATA):
#
		SUB_PLOT     = int(SUB_PLOT0+str(iDATA+1))
		plt.subplot(SUB_PLOT)
		DATA_HIST    = DATA[iDATA,:]
		NUMBER,BINS,PATCHES = plt.hist(DATA_HIST,bins=RANGE,normed=0,facecolor=FACECOLOUR[iDATA],alpha=0.75)
#
		TOTAL          = NUMBER.sum()
		MEAN           = DATA_HIST[DATA_HIST > 0.0].mean()
		STD            = np.std(DATA_HIST[DATA_HIST > 0.0])
		print(SUB_PLOT,len(DATA_HIST),TOTAL,MEAN,DATA_HIST.max(),MAX_COUNT,PLOT_OPT)
		plt.axis([DATA_MIN,DATA_MAX,0.0,MAX_COUNT])
#
		TEXT           = 'Counts = '+('%d' % TOTAL)+'; Mean = '+('%6.1f' % MEAN)+' ppb; std. dev = '+('%6.1f' % STD)+' ppb'
		plt.title(TEXT,fontsize=8)
		plt.xlabel(XLABEL,fontsize=10)
		plt.ylabel(YLABEL,fontsize=10)
		plt.xticks(XTICKS,XTICKS)
		plt.yticks(YTICKS,YTICKS)
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Linear regression
# *****************************************************************************************
#
def Plot_Regression(X,Y,DATA_MIN,DATA_MAX,XLABEL,YLABEL,SLOPE,INTERCEPT,CORREL, \
	XTICKS,YTICKS,PLOT_TITLE,PLOT_OPT,FILE_PLOT,DEBUG):
# 
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 10 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(8.0,8.0))
	plt.suptitle(PLOT_TITLE,fontsize=12)
	plt.xlim(DATA_MIN,DATA_MAX)
	plt.xlabel(XLABEL,fontsize=10)
	plt.xticks(XTICKS,XTICKS)
	plt.ylim(DATA_MIN,DATA_MAX)
	plt.ylabel(YLABEL,fontsize=10)
	plt.yticks(YTICKS,YTICKS)
	plt.plot(X,Y,'ro',markersize=2,markeredgecolor='r')
#
# Derive end points for best-fit line
#
	YMIN         = SLOPE*DATA_MIN+INTERCEPT
	XMIN         = DATA_MIN
#
	if YMIN < DATA_MIN:
		XMIN         = (DATA_MIN-INTERCEPT)/SLOPE
		YMIN         = DATA_MIN
#
	YMAX         = SLOPE*DATA_MAX+INTERCEPT
	XMAX         = DATA_MAX
#
	if YMAX > DATA_MAX:
		XMAX         = (DATA_MAX-INTERCEPT)/SLOPE
		YMAX         = DATA_MAX
#
	plt.plot([XMIN,XMAX],[YMIN,YMAX],'k-')
#
	ONE2ONE      = [DATA_MIN,DATA_MAX]
	plt.plot(ONE2ONE,ONE2ONE,'k--')
#
	if INTERCEPT >= 0.0:
		FIT          = ('Fit y = %.4f x + %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
	else:
		FIT          = ('Fit y = %.4f x %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
#
	plt.legend(('Data',FIT,'1:1 line'),loc=2)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Time series
# *****************************************************************************************
#
def Plot_TimeSeries(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	PLOT_TITLE,LEGEND,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(6.0,6.0))
#
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	if YMIN < 0:
		plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	plt.legend((LEGEND),loc=2)
	plt.suptitle(PLOT_TITLE,fontsize=10)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries2(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(6.0,6.0))
#
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print iDATA
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS)
	plt.suptitle(PLOT_TITLE,fontsize=10)
#
	if YMIN < 0:
		plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Date(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL,POWER, \
	PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	FIGURE      = plt.figure(figsize=(6.0,6.0))
	TIME_MAJOR  = dates.DayLocator(interval=7)
	TIME_MINOR  = dates.DayLocator()
	TIME_FMT    = dates.DateFormatter('%d/%m')
	AXIS        = FIGURE.add_subplot(111)
	AXIS.xaxis.set_major_locator(TIME_MAJOR)
	AXIS.xaxis.set_minor_locator(TIME_MINOR)
	AXIS.xaxis.set_major_formatter(TIME_FMT)
#
	plt.xlim(XMIN,XMAX)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS*POWER,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
	        if DEBUG == 'Y':
	                print(X[iDATA][:])
	                print(Y[iDATA,:])
#
	        plt.plot_date(X[iDATA][:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	plt.suptitle(PLOT_TITLE,fontsize=10)
#
	if PLOT_OPT == '1':
	        plt.show()
	elif PLOT_OPT == '2':
	        print(FILE_PLOT)
	        plt.savefig(FILE_PLOT)
	        plt.close()
#
	return
#
def Plot_TimeSeries_Multi(NPLOTS,NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	ORIENTATION,PLOT_TITLE,SUBTITLES,LEGEND,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NPLOTS != 8:
		print('Incorrect number of plots')
		quit()
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
		NUM_PER_ROW  = 4
		SUB_PLOT0    = str(NUM_PER_ROW)+str(int(NPLOTS/NUM_PER_ROW))
		NPLOT_XLABEL = 5
	else:
		plt.figure(figsize=(8.0,12.0))
		NUM_PER_ROW  = 2
		SUB_PLOT0    = str(int(NPLOTS/NUM_PER_ROW))+str(NUM_PER_ROW)
		NPLOT_XLABEL = 7
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		SUB_PLOT     = int(SUB_PLOT0+str(iPLOT+1))
		plt.subplot(SUB_PLOT)
#
		plt.title(SUB_TITLE,fontsize=10)
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICKS)
		if iPLOT+1 >= NPLOT_XLABEL: plt.xlabel(XLABEL,fontsize=10)
		plt.ylim(YMIN,YMAX)
		plt.yticks(YTICKS,YTICKS)
		if NUM_PER_ROW*int(iPLOT/NUM_PER_ROW) ==  iPLOT: plt.ylabel(YLABEL,fontsize=10)
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
#
		plt.legend((LEGEND),loc=2)
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Multi2(NPLOTS,NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	ORIENTATION,PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NPLOTS != 8:
		print('Incorrect number of plots')
		quit()
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
		NUM_PER_ROW  = 4
		SUB_PLOT0    = str(NUM_PER_ROW)+str(int(NPLOTS/NUM_PER_ROW))
		NPLOT_XLABEL = 5
	else:
		plt.figure(figsize=(8.0,12.0))
		NUM_PER_ROW  = 2
		SUB_PLOT0    = str(int(NPLOTS/NUM_PER_ROW))+str(NUM_PER_ROW)
		NPLOT_XLABEL = 7
#
	for iPLOT in range(NPLOTS):
#
		LINES       = []
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		SUB_PLOT     = int(SUB_PLOT0+str(iPLOT+1))
		plt.subplot(SUB_PLOT)
#
		plt.title(SUB_TITLE,fontsize=10)
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICKS)
		if iPLOT+1 >= NPLOT_XLABEL: plt.xlabel(XLABEL,fontsize=8)
		plt.ylim(YMIN,YMAX)
		plt.yticks(YTICKS,YTICKS)
		if NUM_PER_ROW*int(iPLOT/NUM_PER_ROW) ==  iPLOT: plt.ylabel(YLABEL,fontsize=8)
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA],markersize=2)
			LINES.append(L1)
#
		plt.plot([XMIN,XMAX],[0.0,0.0],'k-')
#		plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	plt.figlegend((LINES),(LEGEND),'lower center',prop=FONT)
#
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Regression(NDATA,X,Y,XREGRESS,YREGRESS,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	SLOPE,INTERCEPT,CORREL,PLOT_TITLE,LEGEND,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(12.0,6.0))
#
# Time series
#
	plt.subplot(121)
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=2)
#
# Regression
#
	DATA_MIN     = YMIN
	DATA_MAX     = YMAX
#
	plt.subplot(122)
	plt.xlabel(YLABEL+': '+LEGEND[0],fontsize=10)
	plt.ylabel(YLABEL+': '+LEGEND[1],fontsize=10)
	plt.xlim(DATA_MIN,DATA_MAX)
	plt.ylim(DATA_MIN,DATA_MAX)
	plt.xticks(YTICKS,YTICKS)
	plt.yticks(YTICKS,YTICKS)
	plt.plot(XREGRESS,YREGRESS,'ro',markeredgecolor='r')
#
# Derive end points for best-fit line
#
	YMIN         = SLOPE*DATA_MIN+INTERCEPT
	XMIN         = DATA_MIN
#
	if YMIN < DATA_MIN:
		XMIN         = (DATA_MIN-INTERCEPT)/SLOPE
		YMIN         = DATA_MIN
#
	YMAX         = SLOPE*DATA_MAX+INTERCEPT
	XMAX         = DATA_MAX
#
	if YMAX > DATA_MAX:
		XMAX         = (DATA_MAX-INTERCEPT)/SLOPE
		YMAX         = DATA_MAX
#
	plt.plot([XMIN,XMAX],[YMIN,YMAX],'k-')
#
	ONE2ONE      = [DATA_MIN,DATA_MAX]
	plt.plot(ONE2ONE,ONE2ONE,'k--')
#
	if INTERCEPT >= 0.0:
		FIT          = ('Fit y = %.4f x + %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
	else:
		FIT          = ('Fit y = %.4f x %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
#
	plt.legend(('Data',FIT,'1:1 line'),loc=2)
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Regression2(NDATA,X,Y,XREGRESS,YREGRESS,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	SLOPE,INTERCEPT,CORREL,PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(12.0,6.0))
#
# Time series
#
	plt.subplot(121)
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS)
#
# Regression
#
	DATA_MIN     = YMIN
	DATA_MAX     = YMAX
#
	plt.subplot(122)
	plt.xlabel(YLABEL+': '+LEGEND[0],fontsize=10)
	plt.ylabel(YLABEL+': '+LEGEND[1],fontsize=10)
	plt.xlim(DATA_MIN,DATA_MAX)
	plt.ylim(DATA_MIN,DATA_MAX)
	plt.xticks(YTICKS,YTICKS)
	plt.yticks(YTICKS,YTICKS)
	plt.plot(XREGRESS,YREGRESS,'ro',markeredgecolor='r')
#
# Derive end points for best-fit line
#
	YMIN         = SLOPE*DATA_MIN+INTERCEPT
	XMIN         = DATA_MIN
#
	if YMIN < DATA_MIN:
		XMIN         = (DATA_MIN-INTERCEPT)/SLOPE
		YMIN         = DATA_MIN
#
	YMAX         = SLOPE*DATA_MAX+INTERCEPT
	XMAX         = DATA_MAX
#
	if YMAX > DATA_MAX:
		XMAX         = (DATA_MAX-INTERCEPT)/SLOPE
		YMAX         = DATA_MAX
#
	plt.plot([XMIN,XMAX],[YMIN,YMAX],'k-')
#
	ONE2ONE      = [DATA_MIN,DATA_MAX]
	plt.plot(ONE2ONE,ONE2ONE,'k--')
#
	if INTERCEPT >= 0.0:
		FIT          = ('Fit y = %.4f x + %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
	else:
		FIT          = ('Fit y = %.4f x %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
#
	plt.legend(('Data',FIT,'1:1 line'),loc=LEGEND_POS)
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Power(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL,LEGEND,PLOT_CODE, \
	NDATA2,X2,Y2,XMIN2,XMAX2,XTICKS2,XLABEL2,YMIN2,YMAX2,YTICKS2,YLABEL2,LEGEND2,PLOT_CODE2,
	LEGEND_POS,PLOT_TITLE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(12.0,6.0))
#
# Time series
#
	plt.subplot(121)
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
# Power spectrum/Growth Rate
#
	plt.subplot(122)
	plt.xlim(XMIN2,XMAX2)
	plt.xticks(XTICKS2,XTICKS2)
	plt.xlabel(XLABEL2,fontsize=10)
#
	plt.ylim(YMIN2,YMAX2)
	plt.yticks(YTICKS2,YTICKS2)
	plt.ylabel(YLABEL2,fontsize=10)
#
	for iDATA in range(NDATA2):
		if DEBUG == 'Y':
			print(X2[iDATA,:])
			print(Y2[iDATA,:])
#
		plt.plot(X2[iDATA,:],Y2[iDATA,:],PLOT_CODE2[iDATA])
#
	plt.plot([XMIN,XMAX],[0.0,0.0],'k-')
#
	plt.legend((LEGEND2),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	plt.suptitle(PLOT_TITLE,fontsize=12)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_TimeSeries_Power_Multi(NPLOTS,NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL,
	NDATA2,X2,Y2,XMIN2,XMAX2,XTICKS2,XLABEL2,YMIN2,YMAX2,YTICKS2,YLABEL2, \
	SITE_CODES,ORIENTATION,FONTSIZES,SUBTITLES,LEGEND,PLOT_CODE,LEGEND2,PLOT_CODE2, \
	LEGEND_POS,PLOT_TITLE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if 2*NPLOTS != 8:
		print('Incorrect number of plots')
		quit()
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
		NUM_PER_ROW  = 4
		SUB_PLOT0    = str(NUM_PER_ROW)+str(int(2*NPLOTS/NUM_PER_ROW))
		NPLOT_XLABEL = 5
	else:
		plt.figure(figsize=(8.0,12.0))
		NUM_PER_ROW  = 2
		SUB_PLOT0    = str(int(2*NPLOTS/NUM_PER_ROW))+str(NUM_PER_ROW)
		NPLOT_XLABEL = 7
#
	for iPLOT in range(NPLOTS):
#
		iSITE        = SITE_CODES[iPLOT]
		SUB_TITLE    = SUBTITLES[iPLOT]
		SUB_PLOT     = int(SUB_PLOT0+str(2*iPLOT+1))
		plt.subplot(SUB_PLOT)
# Time series
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICKS)
		if 2*iPLOT+1 >= NPLOT_XLABEL: plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#
		plt.ylim(YMIN,YMAX)
		plt.yticks(YTICKS,YTICKS)
		plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#
		for iDATA in range(NDATA):
			if DEBUG == 'Y':
				print(X[iSITE,iDATA,:])
				print(Y[iSITE,iDATA,:])
#
			plt.plot(X[iSITE,iDATA,:],Y[iSITE,iDATA,:],PLOT_CODE[iDATA])
#
		plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
# Power spectrum/Growth Rate
#
		SUB_PLOT     = int(SUB_PLOT0+str(2*iPLOT+2))
		plt.subplot(SUB_PLOT)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN2,XMAX2)
		plt.xticks(XTICKS2,XTICKS2)
		if 2*iPLOT+1 >= NPLOT_XLABEL: plt.xlabel(XLABEL2,fontsize=FONTSIZES[1])
#
		plt.ylim(YMIN2,YMAX2)
		plt.yticks(YTICKS2,YTICKS2)
		plt.ylabel(YLABEL2,fontsize=FONTSIZES[1])
#
		for iDATA in range(NDATA2):
			if DEBUG == 'Y':
				print(X2[iSITE,iDATA,:])
				print(Y2[iSITE,iDATA,:])
#
			plt.plot(X2[iSITE,iDATA,:],Y2[iSITE,iDATA,:],PLOT_CODE2[iDATA])
#
		plt.plot([XMIN,XMAX],[0.0,0.0],'k-')
#
		plt.legend((LEGEND2),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Taylor plot
# *****************************************************************************************
#
def Plot_Taylor(TAYLOR_R,TAYLOR_A,RADII,PLOT_TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT     = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	PI       = np.pi
	XMAX     = RADII[-1]
#
	if DEBUG == 'Y':
		print(TAYLOR_A)
		print(TAYLOR_R)
#
	if TAYLOR_A.max() > PI/2:
		plt.figure(figsize=(12.0,6.0))
		plt.xlim(-XMAX,XMAX)
		ANGLEE    = 180
		CORR      = [ -0.99,-0.95,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, \
			       0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
#
		NRADII    = len(RADII)
		XRADII    = arange(2*NRADII)
		LRADII    = arange(2*NRADII)
#
		for i in range(NRADII):
			XRADII[i]        = -int(RADII[i])
			XRADII[i+NRADII] =  int(RADII[i])
			LRADII[i]        =  int(RADII[i])
			LRADII[i+NRADII] =  int(RADII[i])
#
	else:
		plt.figure(figsize=(6.0,6.0))
		plt.xlim(0.0,XMAX)
		ANGLEE    =  90
		CORR      = [ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
		XRADII    = RADII
		LRADII    = RADII
#
	ANGLES   = 0
#
	for iRAD   in range(1,XMAX+1):
		X       = []
		Y       = []
		for iANGLE in range(ANGLES,ANGLEE+1):
			X.append(iRAD*math.cos(math.radians(iANGLE)))
			Y.append(iRAD*math.sin(math.radians(iANGLE)))
		X       = np.array(X)
		Y       = np.array(Y)
		plt.plot(X,Y,'k:')
#
	YCORR   = np.sin(np.arccos(CORR))
#
	for iDATA in range(len(CORR)):
		plt.plot([0.0,CORR[iDATA]*XMAX],[0.0,YCORR[iDATA]*XMAX],'k:')
		plt.text(0.95*CORR[iDATA]*XMAX,0.95*YCORR[iDATA]*XMAX,str(CORR[iDATA]),fontsize=FONTSIZES[1])
#
	plt.plot([0.0,0.0],[0.0,XMAX],'k-')
#
	NDATA    = len(TAYLOR_R)
	TAYLOR_X = np.zeros(NDATA)
	TAYLOR_Y = np.zeros(NDATA)
# 
	TAYLOR_X = TAYLOR_R*np.cos(TAYLOR_A)
	TAYLOR_Y = TAYLOR_R*np.sin(TAYLOR_A)
#
	plt.ylim(0.0,XMAX)
	plt.yticks(RADII,RADII)
	plt.xticks(XRADII,LRADII)
	plt.plot(TAYLOR_X,TAYLOR_Y,PLOT_CODE)
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_Taylor_Multi(NRUNS,TAYLOR_R,TAYLOR_A,RADII,PLOT_TITLE,LEGEND,LEGEND_POS,FONTSIZES,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT     = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	PI       = np.pi
	XMAX     = RADII[-1]
#
	if DEBUG == 'Y':
		print(TAYLOR_A)
		print(TAYLOR_R)
#
	NDATA    = TAYLOR_R.shape[0]
	AMAX     = TAYLOR_A[~np.isnan(TAYLOR_A)].max()
#
	print TAYLOR_R.shape,TAYLOR_A.shape,AMAX
#
	if AMAX > PI/2:
		plt.figure(figsize=(12.0,6.0))
		plt.xlim(-XMAX,XMAX)
		ANGLEE    = 180
		CORR      = [ -0.99,-0.95,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, \
			       0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
#
		NRADII    = len(RADII)
		XRADII    = arange(2*NRADII)
		LRADII    = arange(2*NRADII)
#
		for i in range(NRADII):
			XRADII[i]        = -int(RADII[i])
			XRADII[i+NRADII] =  int(RADII[i])
			LRADII[i]        =  int(RADII[i])
			LRADII[i+NRADII] =  int(RADII[i])
#
	else:
		plt.figure(figsize=(6.0,6.0))
		plt.xlim(0.0,XMAX)
		ANGLEE    =  90
		CORR      = [ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
		XRADII    = RADII
		LRADII    = RADII
#
	for iRUN in range(NRUNS):
		TAYLOR_X = TAYLOR_R[:,iRUN]*np.cos(TAYLOR_A[:,iRUN])
		TAYLOR_Y = TAYLOR_R[:,iRUN]*np.sin(TAYLOR_A[:,iRUN])
		plt.plot(TAYLOR_X,TAYLOR_Y,PLOT_CODE[iRUN])
#
	print LEGEND
	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	ANGLES   = 0
#
	for iRAD   in range(1,XMAX+1):
		X       = []
		Y       = []
		for iANGLE in range(ANGLES,ANGLEE+1):
			X.append(iRAD*math.cos(math.radians(iANGLE)))
			Y.append(iRAD*math.sin(math.radians(iANGLE)))
		X       = np.array(X)
		Y       = np.array(Y)
		plt.plot(X,Y,'k:')
#
	YCORR   = np.sin(np.arccos(CORR))
#
	for iDATA in range(len(CORR)):
		plt.plot([0.0,CORR[iDATA]*XMAX],[0.0,YCORR[iDATA]*XMAX],'k:')
		plt.text(0.95*CORR[iDATA]*XMAX,0.95*YCORR[iDATA]*XMAX,str(CORR[iDATA]),fontsize=FONTSIZES[1])
#
	plt.plot([0.0,0.0],[0.0,XMAX],'k-')
#
	plt.ylim(0.0,XMAX)
	plt.yticks(RADII,RADII)
	plt.xticks(XRADII,LRADII)
#
	TAYLOR_X = np.zeros(NDATA)
	TAYLOR_Y = np.zeros(NDATA)
# 
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_Taylor_Polar(TAYLOR_R,TAYLOR_A,RADII,PLOT_TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT     = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	PI       = np.pi
#
	if DEBUG == 'Y':
		print(TAYLOR_A)
		print(TAYLOR_R)
#
	if TAYLOR_A.max() > PI/2:
		CORR      = [ -0.99,-0.95,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, \
			       0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
#
	else:
		CORR      = [ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99 ]
#
	plt.figure(figsize=(6.0,6.0))
	plt.polar(TAYLOR_A,TAYLOR_R,PLOT_CODE)
	plt.rgrids(RADII,RADII,angle=0.0,fontsize=FONTSIZES[1])
	plt.thetagrids(np.degrees((np.arccos(CORR))),labels=CORR,frac=0.95,fontsize=FONTSIZES[1])
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[12])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Zonal plot
# *****************************************************************************************
#
def Plot_Zonal(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	PLOT_TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG,WIDTH=6,HEIGHT=6, \
	LEGEND_POS=0):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(WIDTH,HEIGHT))
#
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#
	plt.title(PLOT_TITLE,fontsize=FONTSIZES[2])
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_Zonal_Regression(NDATA,X,Y,XREGRESS,YREGRESS, \
	XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	SLOPE,INTERCEPT,CORREL,PLOT_TITLE,LEGEND,FONTSIZES,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(12.0,6.0))
#
	plt.subplot(121)
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICKS)
	plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=2)
#
# Regression
#
	DATA_MIN     = YMIN
	DATA_MAX     = YMAX
#
	plt.subplot(122)
	plt.xlabel(YLABEL+': '+LEGEND[0],fontsize=FONTSIZES[1])
	plt.ylabel(YLABEL+': '+LEGEND[1],fontsize=FONTSIZES[1])
	plt.xlim(DATA_MIN,DATA_MAX)
	plt.ylim(DATA_MIN,DATA_MAX)
	plt.xticks(YTICKS,YTICKS)
	plt.yticks(YTICKS,YTICKS)
	plt.plot(XREGRESS,YREGRESS,'ro',markeredgecolor='r')
#
# Derive end points for best-fit line
#
	YMIN         = SLOPE*DATA_MIN+INTERCEPT
	XMIN         = DATA_MIN
#
	if YMIN < DATA_MIN:
		XMIN         = (DATA_MIN-INTERCEPT)/SLOPE
		YMIN         = DATA_MIN
#
	YMAX         = SLOPE*DATA_MAX+INTERCEPT
	XMAX         = DATA_MAX
#
	if YMAX > DATA_MAX:
		XMAX         = (DATA_MAX-INTERCEPT)/SLOPE
		YMAX         = DATA_MAX
#
	plt.plot([XMIN,XMAX],[YMIN,YMAX],'k-')
#
	ONE2ONE      = [DATA_MIN,DATA_MAX]
	plt.plot(ONE2ONE,ONE2ONE,'k--')
#
	if INTERCEPT >= 0.0:
		FIT          = ('Fit y = %.4f x + %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
	else:
		FIT          = ('Fit y = %.4f x %.4f, R = %.4f' % (SLOPE,INTERCEPT,CORREL))
#
	plt.legend(('Data',FIT,'1:1 line'),loc=2)
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Climatologies
# *****************************************************************************************
#
def Plot_Climatology(NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	PLOT_TITLE,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : 8 }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	plt.figure(figsize=(6.0,6.0))
#
	XTICK_LABS  = ['J','F','M','A','M','J','J','A','S','O','N','D' ]
#
	plt.xlim(XMIN,XMAX)
	plt.xticks(XTICKS,XTICK_LABS)
	plt.xlabel(XLABEL,fontsize=10)
#
	plt.ylim(YMIN,YMAX)
	plt.yticks(YTICKS,YTICKS)
	plt.ylabel(YLABEL,fontsize=10)
#
	for iDATA in range(NDATA):
		if DEBUG == 'Y':
			print(X[iDATA,:])
			print(Y[iDATA,:])
#
		plt.plot(X[iDATA,:],Y[iDATA,:],PLOT_CODE[iDATA])
#
	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	plt.suptitle(PLOT_TITLE,fontsize=10)
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_Climatology_Multi(NPLOTS,NDATA,X,Y,XMIN,XMAX,XTICKS,XLABEL,YMIN,YMAX,YTICKS,YLABEL, \
	ORIENTATION,FONTSIZES,PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	XTICK_LABS  = ['J','F','M','A','M','J','J','A','S','O','N','D' ]
#
	if NPLOTS != 8:
		print('Incorrect number of plots')
		quit()
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
		NUM_PER_ROW  = 4
		SUB_PLOT0    = str(NUM_PER_ROW)+str(int(NPLOTS/NUM_PER_ROW))
		NPLOT_XLABEL = 5
	else:
		plt.figure(figsize=(8.0,12.0))
		NUM_PER_ROW  = 2
		SUB_PLOT0    = str(int(NPLOTS/NUM_PER_ROW))+str(NUM_PER_ROW)
		NPLOT_XLABEL = 7
#
#
	for iPLOT in range(NPLOTS):
#
		LINES       = []
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		SUB_PLOT     = int(SUB_PLOT0+str(iPLOT+1))
		plt.subplot(SUB_PLOT)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT+1 >= NPLOT_XLABEL: plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
		plt.ylim(YMIN,YMAX)
		plt.yticks(YTICKS,YTICKS)
		if NUM_PER_ROW*int(iPLOT/NUM_PER_ROW) ==  iPLOT: plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
#		plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
#
	plt.figlegend((LINES),(LEGEND),'lower center',prop=FONT)
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
# Generic multi-subplot multi-run plot
# *****************************************************************************************
#
def Plot_General_MultiPlot(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,ORIENTATION,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG,NSIZE=5.0):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
#	DEBUG       = 'Y'
#
	if ORIENTATION == 'Landscape':
		plt.figure(figsize=(NROWS*NSIZE,NCOLUMNS*NSIZE))
	else:
		plt.figure(figsize=(NCOLUMNS*NSIZE,NROWS*NSIZE))
#
#
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN,YMAX)
		plt.yticks(YTICKS,YTICK_LABS)
		if NCOLUMNS*int(iPLOT/NCOLUMNS) == iPLOT:
			plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,YLABEL
#
		LINES       = []
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
#
			if DEBUG == 'Y':
				print
				print iPLOT,iDATA
				print(XPLOT)
				print(YPLOT)
#
			if YMAX > YMIN:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMAX)] = YMAX
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMIN)] = YMIN
			else:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMAX)] = YMAX
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMIN)] = YMIN
#
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_General_MultiPlot_VarScale(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,YSCALE,ORIENTATION,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NROWS*NCOLUMNS == 1:
		plt.figure(figsize=(6.0,6.0))
	elif ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
	else:
		plt.figure(figsize=(8.0,12.0))
#
#	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN*YSCALE[iPLOT],YMAX*YSCALE[iPLOT])
		plt.yticks(YTICKS*YSCALE[iPLOT],YTICK_LABS*YSCALE[iPLOT])
		if NCOLUMNS*int(iPLOT/NCOLUMNS) == iPLOT:
			plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,YLABEL
#
		LINES       = []
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
#
			if YMAX > YMIN:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMAX*YSCALE[iPLOT])] = YMAX*YSCALE[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMIN*YSCALE[iPLOT])] = YMIN*YSCALE[iPLOT]
			else:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMAX*YSCALE[iPLOT])] = YMAX*YSCALE[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMIN*YSCALE[iPLOT])] = YMIN*YSCALE[iPLOT]
#
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_General_MultiPlot_VarScale2(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,ORIENTATION,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NROWS*NCOLUMNS == 1:
		plt.figure(figsize=(6.0,6.0))
	elif ORIENTATION == 'Landscape':
		plt.figure(figsize=(12.0,8.0))
	else:
		plt.figure(figsize=(8.0,12.0))
#
#	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN[iPLOT],YMAX[iPLOT])
		plt.yticks(YTICKS[iPLOT],YTICK_LABS[iPLOT])
		plt.ylabel(YLABEL[iPLOT],fontsize=FONTSIZES[1])
#
		LINES       = []
#
		for iDATA in range(NDATA[iPLOT]):
			XPLOT       = np.array(X[iPLOT][iDATA])
			YPLOT       = np.array(Y[iPLOT][iDATA])
#
			if YMAX[iPLOT] > YMIN[iPLOT]:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMIN[iPLOT])] = YMIN[iPLOT]
			else:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMIN[iPLOT])] = YMIN[iPLOT]
#
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN[iPLOT] < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_General_MultiPlot_VarScale2a(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,WIDTH,HEIGHT,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NROWS*NCOLUMNS == 1:
		plt.figure(figsize=(6.0,6.0))
	else:
		plt.figure(figsize=(WIDTH,HEIGHT))
#
#	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN[iPLOT],YMAX[iPLOT])
		plt.yticks(YTICKS[iPLOT],YTICK_LABS[iPLOT])
		plt.ylabel(YLABEL[iPLOT],fontsize=FONTSIZES[1])
#
		LINES       = []
#
		for iDATA in range(NDATA[iPLOT]):
			XPLOT       = np.array(X[iPLOT][iDATA])
			YPLOT       = np.array(Y[iPLOT][iDATA])
#
			if DEBUG == 'Y':
				print iPLOT,iDATA
				print(XPLOT)
				print(YPLOT)
#
			if YMAX[iPLOT] > YMIN[iPLOT]:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMIN[iPLOT])] = YMIN[iPLOT]
			else:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMIN[iPLOT])] = YMIN[iPLOT]
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN[iPLOT] < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_General_MultiPlot_VarScale2b(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,WIDTH,HEIGHT,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,PLOT_TEXT,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NROWS*NCOLUMNS == 1:
		plt.figure(figsize=(6.0,6.0))
	else:
		plt.figure(figsize=(WIDTH,HEIGHT))
#
#	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN[iPLOT],YMAX[iPLOT])
		plt.yticks(YTICKS[iPLOT],YTICK_LABS[iPLOT])
		plt.ylabel(YLABEL[iPLOT],fontsize=FONTSIZES[1])
#
		if len(PLOT_TEXT) != 0:
			plt.annotate(PLOT_TEXT[iPLOT],xy=(0.02,0.92),xycoords='axes fraction')

		LINES       = []
#
		for iDATA in range(NDATA[iPLOT]):
			XPLOT       = np.array(X[iPLOT][iDATA])
			YPLOT       = np.array(Y[iPLOT][iDATA])
#
			if YMAX[iPLOT] > YMIN[iPLOT]:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMIN[iPLOT])] = YMIN[iPLOT]
			else:
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT < YMAX[iPLOT])] = YMAX[iPLOT]
				YPLOT[(~np.isnan(YPLOT)) & (YPLOT > YMIN[iPLOT])] = YMIN[iPLOT]
#
			if DEBUG == 'Y':
				print iPLOT,iDATA
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN[iPLOT] < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
def Plot_General_MultiPlot_VarScale3(NROWS,NCOLUMNS,NDATA,X,Y,XMIN,XMAX,XTICKS,XTICK_LABS,XLABEL, \
	YMIN,YMAX,YTICKS,YTICK_LABS,YLABEL,YSCALE,WIDTH,HEIGHT,FONTSIZES, \
	PLOT_TITLE,SUBTITLES,PLOT_TEXT,LEGEND,LEGEND_POS,PLOT_CODE,PLOT_OPT,FILE_PLOT,DEBUG):
#
	FONT        = {'family' : 'sans-serif', 'weight' : 'normal', 'size' : FONTSIZES[0] }
	rc('font', **FONT)
	rc('text',usetex=True)
#
	if NROWS*NCOLUMNS == 1:
		plt.figure(figsize=(6.0,6.0))
	else:
		plt.figure(figsize=(WIDTH,HEIGHT))
#
#	plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
	NPLOTS      = NROWS*NCOLUMNS
	NPLOT_XLABS = NPLOTS-NCOLUMNS
#	print NPLOTS,NPLOT_XLABS
#
	for iPLOT in range(NPLOTS):
#
		SUB_TITLE    = SUBTITLES[iPLOT]
		plt.subplot(NROWS,NCOLUMNS,iPLOT+1)
#
		plt.title(SUB_TITLE,fontsize=FONTSIZES[2])
		plt.xlim(XMIN,XMAX)
		plt.xticks(XTICKS,XTICK_LABS)
		if iPLOT >= NPLOT_XLABS:
			plt.xlabel(XLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,XLABEL
		plt.ylim(YMIN*YSCALE[iPLOT],YMAX*YSCALE[iPLOT])
		plt.yticks(YTICKS*YSCALE[iPLOT],YTICK_LABS*YSCALE[iPLOT])
		if NCOLUMNS*int(iPLOT/NCOLUMNS) == iPLOT:
			plt.ylabel(YLABEL,fontsize=FONTSIZES[1])
#			print iPLOT,YLABEL
#
		if len(PLOT_TEXT) != 0:
			plt.annotate(PLOT_TEXT[iPLOT],xy=(5,5),xycoords='axes points',fontsize=FONTSIZES[0])

		LINES       = []
#
		for iDATA in range(NDATA):
			XPLOT       = np.array(X[NDATA*iPLOT+iDATA])
			YPLOT       = np.array(Y[NDATA*iPLOT+iDATA])
#
			if YMAX > YMIN:
				YPLOT[YPLOT > YMAX*YSCALE[iPLOT]] = YMAX*YSCALE[iPLOT]
				YPLOT[YPLOT < YMIN*YSCALE[iPLOT]] = YMIN*YSCALE[iPLOT]
			else:
				YPLOT[YPLOT < YMAX*YSCALE[iPLOT]] = YMAX*YSCALE[iPLOT]
				YPLOT[YPLOT > YMIN*YSCALE[iPLOT]] = YMIN*YSCALE[iPLOT]
#
			if DEBUG == 'Y':
				print(XPLOT)
				print(YPLOT)
#
			L1,         = plt.plot(XPLOT,YPLOT,PLOT_CODE[iDATA])
			LINES.append(L1)
#
		if YMIN < 0:
			plt.plot([XMIN,XMAX],[0.0,0.0],'k--')
#
	if LEGEND_POS >= 0:
		if NPLOTS == 1:
			plt.legend((LEGEND),loc=LEGEND_POS,frameon='False',prop=FONT)
		else:
			plt.figlegend((LINES),(LEGEND),'lower center',frameon='False',prop=FONT,ncol=min(2,NDATA))
#
	plt.suptitle(PLOT_TITLE,fontsize=FONTSIZES[3])
#
	if PLOT_OPT == '1':
		plt.show()
	elif PLOT_OPT == '2':
		print(FILE_PLOT)
		plt.savefig(FILE_PLOT)
		plt.close()
#
	return
#
# *****************************************************************************************
#
def get_ticks(MAX0,MIN0,NTICKS,DEBUG):
#
	MAX         = MAX0
#
	if MAX >= 1.0:
#
		POWER       =  1.0
		while MAX > 1.0:
			POWER       = 10.0*POWER
			MAX         = MAX/POWER
#
	elif MAX < 1.0:
#
		POWER       = 10.0
		while MAX < 1.0:
			POWER       = POWER/10.0
			MAX         = MAX/POWER

#
	if DEBUG == 'Y': print POWER,MAX,MAX0 % POWER
#
	if MAX0 % POWER == 0:
		MAX         = POWER*int(MAX0/POWER)
	else:
		MAX         = POWER*(1+int(MAX0/POWER))
#
	MIN         = POWER*int(MIN0/POWER)
#		
	if MAX0 > 1:
		MAX         = int(MAX)
		MIN         = int(MIN)
#
	INC         = (MAX-MIN)/(NTICKS-1)
#
	if DEBUG == 'Y': print POWER,MIN,MAX,INC
	print POWER,MIN,MAX,INC
#
	TICKS       = np.arange(MIN,MAX+INC,INC)
#
	return TICKS
#
def get_ticks2(MAX0,MIN0,NTICKS,DEBUG):
#
	INC         = (MAX0-MIN0)/(NTICKS-1)
	TICKS       = MIN0 + INC*np.arange(NTICKS)
#
	return TICKS
#
# End of Program
