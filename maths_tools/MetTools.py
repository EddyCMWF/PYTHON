########################################
#  ECP - Python Meterological Tools
########################################
import numpy as np

Mw=18.0160 # molecular weight of water
Md=28.9660 # molecular weight of dry air
Mrat = Mw/Md  # molecular mass ratios
R =  8.31432E3 # gas constant
Rd = R/Md # specific gas constant for dry air
Rv = R/Mw # specific gas constant for vapour
Tref = 273.16
Pref = 101325.0

# air density
def rhoair(p,T):
    '''calulate air density from pressure (Pa) and Temperature (K) using ideal gas law
    '''
    return p/(T*Rd*1e3)

#saturation pressure
def esat(T):
    ''' get saturation pressure (units [Pa]) for a given air temperature (units [K])'''
    TK = Tref
    e1 = Pref
    logTTK = np.log10(T/TK)
    esat =  e1*10**(10.79586*(1-TK/T)-5.02808*logTTK+ \
                1.50474*1e-4*(1.-10**(-8.29692*(T/TK-1)))+ \
                0.42873*1e-3*(10**(4.76955*(1-TK/T))-1)-2.2195983)
    return esat

#def esat(T):
#    ''' get saturation pressure (units [Pa]) for a given air temperature (units [K])'''
#    if T>150: T=T-273.15
#    esat =  610.8 * np.exp(17.27 * T / (T+237.3) )
#    return esat

# conversion mass mixing ratio (mixr) to volume mixing ratio
def mixr2vmr(mixr):
    '''conversion from volume mixing ratio (units [n/n]) to mass mixing ratio (units [kg/kg])
    '''
    return mixr*Mrat

# conversion mixing ratio to specific humidity (always kg/kg)
def mixr2sh(W):
    '''conversion from mixing ratio (units [kg/kg]) to specific humidity (units also [kg/kg])
    '''
    return W/(1.+W)

# conversion specific humidity to mass mixing ratio (always kg/kg)
def sh2mixr(qv):
    '''conversion from specific humidity (units [kg/kg]) to mixing ratio (units also [kg/kg])
    '''
    return qv/(1.-qv)

def vmr2ea(vmr,P=Pref):
    '''conversion from volume mixing ratio (units [n/n]) to water vapour pressure (units [Pa])
        Requires Pressure, if not given then set to reference pressure of 101325 Pa
    '''
    return vmr*P

def ea2vmr(ea,P=Pref):
    '''conversion from volume mixing ratio (units [n/n]) to water vapour pressure (units [Pa])
        Requires Pressure, if not given then set to reference pressure of 101325 Pa
    '''
    return ea/P 

def ea2mixr(ea,P=Pref):
    '''conversion from vapour pressure (both units [Pa]) to mixing ratio (units [kg/kg]) 
        Requires Pressure, if not given then set to reference pressure of 101325 Pa
    '''
    return Mrat*ea/(P-ea)

def mixr2ea(mixr,P=Pref):
    '''conversion from mixing ratio (units [kg/kg] to vapour pressure (units [Pa])
        Requires Pressure, if not given then set to reference pressure of 101325 Pa
    '''
    return P*mixr/(Mrat+mixr)

# conversion from mixing ratio to relative humidity
def mixr2rh(mixr,P=Pref,T=Tref):
    '''conversion from mixing ratio (units [kg/kg]) to relative humidity (units [%])
        Requires Pressure and temperature, if not given then set to 
           reference pressure of 101325 Pa and temperature of 273.15
    '''
    es = esat(T)
    ea = mixr2ea(mixr,P=P)
    RH = (ea/es) * 100.
    return RH 

# conversion from mixing ratio to relative humidity
def rh2mixr(RH,P=Pref,T=Tref):
    '''conversion from relative humidity (units [%]) to mixing ratio (units [kg/kg])
        Requires Pressure and temperature, if not given then set to 
           reference pressure of 101325 Pa and temperature of 273.15
    '''
    es   = esat(T)
    ea   = (RH/100.)*es
    mixr = ea2mixr(ea,P=P)
    return mixr

def rh2sh(RH,P=Pref,T=Tref):
    '''conversion from relative humidity (units %) to specific humidity (units [kg/kg])  
        Requires Pressure and temperature, if not given then set to 
           reference pressure of 101325 Pa and temperature of 273.15
    '''
    mixr=rh2mixr(RH,P=P,T=T)
    SH=mixr2sh(mixr)
    return SH

def sh2rh(SH,P=Pref,T=Tref):
    '''conversion from relative humidity (units %) to specific humidity (units [kg/kg])  
        Requires Pressure and temperature, if not given then set to 
           reference pressure of 101325 Pa and temperature of 273.15
    '''
    mixr=sh2mixr(SH)
    RH=mixr2rh(mixr,P=P,T=T)
    return RH


def rh2vpd(RH, T=Tref):
    '''conversion from relative humidity (%) to vapour pressure deficit (Pa)
        Requires temperature. 
    '''
    es =  esat(T)  
    ea = (RH/100)*es
    ##print('rh2vpd')
    ##print('RH,T,es,ea = ', RH,T,es,ea)
    return ( es - ea )

def vpd2rh(VPD, T=Tref):
    '''conversion from vapour pressure deficit (Pa) to relative humidity (%) 
        Requires temperature, if not given then set to 
           reference temperature of 273.15
    '''
    es =  esat(T)  
    ea = es - VPD
    ##print('vpd2rh')
    ##print('VPD,T,es,ea = ', VPD,T,es,ea)
    return 100.*ea/es


def sh2vpd(SH, P=Pref, T=Tref):
    '''conversion from specifc humidity ([kg/kg]) to vapour pressure deficit (Pa)
       Requires temperaturei and pressure, if not given then set to 
       reference temperature of 273.15 and pressure of 101325 Pa
    '''
    RH = sh2rh(SH, P=P, T=T)
    ##print('sh2vpd')
    ##print('SH, RH,T,P,esat(T) = ', SH, RH,T,P,esat(T))
    return rh2vpd(RH, T=T)

def vpd2sh(VPD, P=Pref, T=Tref):
    '''conversion from specifc humidity ([kg/kg]) to vapour pressure deficit (Pa)
       Requires temperaturei and pressure, if not given then set to 
       reference temperature of 273.15 and pressure of 101325 Pa
    '''
    RH = vpd2rh(VPD, T=T)
    ##print('vpd2sh')
    ##print('VPD,RH,T,P,esat(T)=',VPD,RH,T,P,esat(T))
    return rh2sh(RH, P=P, T=T)



##################################################################################
# The following is scrap left over from stuff i inherited that I have not checked and simplified

# not tested
def ah2mixr (rhov,p,T):
    '''conversion from absolute humidity (units [kg/m**3]) to mixing ratio (units also [kg/kg])
    '''
    return (Rd * T)/(p/rhov-Rv*T)

# not tested
def wvap2sh(e,p):
    '''conversion from water vapour pressure (units [Pa]) to specific humidity (units [kg/kg]) 
    '''
    return Mrat*e/(p-(1.-Mrat)*e)

def conc2mixr(conc,p=Pref,T=Tref,Msp=Mw):
    '''conversion from concentration (units mol/m**3) to mass mixing ratio (units [kg/kg])
       Msp is the molecular mass of the gas species in question
    '''
    # density of air:
    rho_a = rhoair(p,T)
    return conc * Msp * 1e-3 / rho_a

def conc2sh(conc,p=Pref,T=Tref,Msp=Mw):
    '''conversion from concentration (units mol/m**3) to mass specific humidity (units [kg/kg])
       Msp is the molecular mass of the gas species in question
    '''
    # density of air:
    rho_a = rhoair(p,T)
    return mixr2sh(conc * Msp * 1e-3 / rho_a)





def SVPfromT(T):
    # Calculation from Murray (1967):
    #   "On the computation of saturation vapour pressure"
    #   T must be in Celcius
    SVP = 610.7 * 10**(7.5*T/(237.3+T))  # Units = Pa
    return SVP

def VPDfromSVPnRHum(SVP,RHum):
    # Calculate Vapour Pressure Deficit from Relative Humidity
    # SVP Units = Pa
    # RHum Units = %
    VPD = ( 1. - (RHum/100.))*SVP
    return VPD

def VPDfromTnRHum(T,RHum):
    # T units = Celcius
    # RHum units = %
    VPD = ( 1. - (RHum/100.))*SVPfromT(T)
    return VPD

def RHumfromVPDnSVP(VPD,SVP):
    # VPD and SVP Units = The same as each other (Pa)
    RHum = (1. - (VPD/SVP))*100.
    return RHum

def RHumfromVPDnT(VPD,T):
    # SPD Units = Pa
    # T Units = Celcius
    RHum = (1. - ( VPD/ SVPfromT(T) ))*100.
    return RHum


def vpd2spechum( VPD, T, P, \
                  Mdair=Md, \
                  Mwv=Mw ):
    # VPD Units = Pa
    # T Units = Celcius
    # P Units = Pa
    # Mdair,Mwv Units = g/mol (Molar mass of dry air and water vapour, respectively
    
    # First calculate SVP
    RHum  = RHumfromVPDnT(VPD,T)
    
    

    
