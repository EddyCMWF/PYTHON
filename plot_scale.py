#
# Python module to derive ticks and labels for plots
#
# Garry Hayman
# Centre for Ecology and Hydrology
# October 2013
#
import numpy as np
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
	if MIN0 != 0.0:
		MIN         = -MAX
	else:
		MIN         = MIN0
#
	if MAX0 > 1:
		MAX         = int(MAX)
		MIN         = int(MIN)
#
	INC         = (MAX-MIN)/(NTICKS-1)
#
	if DEBUG == 'Y': print POWER,MIN,MAX,INC
#
	TICKS       = np.arange(MIN,MAX+INC,INC)
#
	return TICKS
