#!/bin/env python2.7


from maths_tools import MetTools as MTs

vpd=1.687*100
T=296.06
P=100431
sh=0.01642
rh=93.9

print('SH = ', MTs.vpd2sh(vpd,T=T,P=P))

print('RH = ', MTs.vpd2rh(vpd,T=T))

print('VPD = ', MTs.sh2vpd(sh,T=T,P=P))

print('VPD = ', MTs.rh2vpd(rh,T=T))

