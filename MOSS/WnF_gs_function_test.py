#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np


def WnF_gc(W):
   gc = (-0.195 + (0.134*W) - (0.0256 * W**2) + (0.00228*W**3)
            - (0.0000984*W**4) + (0.00000168*W**5) )
   return gc

def ECP_gc(W):
   gc = (-0.0943325 + (0.0725185*W) - (0.01171817 * W**2) + (0.0008136875*W**3)
            - (0.000026223156*W**4) + (0.000000326075727*W**5) )
   return gc


data_points = np.loadtxt('WnF_data.dat')
W_pts = data_points[:,0]
gc_pts = data_points[:,1]

W=np.arange(0,17,0.1)
Gc_WnF=WnF_gc(W)
Gc_ECP=ECP_gc(W)

np.polyfit(W_pts,gc_pts,5)

plt.plot(W,Gc_WnF,label='WnF')
plt.plot(W,Gc_ECP,label='ECP')
plt.scatter(W_pts,gc_pts)
plt.ylim(0,0.08)
plt.xlim(0,20)
plt.show()


