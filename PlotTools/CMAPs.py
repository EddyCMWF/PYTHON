#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:19:53 2019

@author: edwcom
"""
from matplotlib import colors

def ESA_CCI_cmap(infile='/prj/CLIFFTOP/ESA_CCI_LC/ESACCI-LC-Legend_sname.csv'):
    inlines = open(infile).readlines()
    header = inlines.pop(0)
    indices, labels, R, G, B = [],[],[],[],[]
    for line in inlines:
        split=line.replace('\r\n','').split(';')
        indices.append(int(split[0]))
        labels.append(split[1])
        R.append(float(split[2])/255.)
        G.append(float(split[3])/255.)
        B.append(float(split[4])/255.)
    
    cmap = colors.ListedColormap([RGB for RGB in zip(R,G,B)])
    norm = colors.BoundaryNorm(indices,cmap.N)
    
    return indices, labels, cmap, norm


def ESA_CCI_cmap_reduced(red_index, infile='/prj/CLIFFTOP/ESA_CCI_LC/ESACCI-LC-Legend_sname.csv'):
    inlines = open(infile).readlines()
    header = inlines.pop(0)
    indices_in, labels_in, R_in, G_in, B_in = [],[],[],[],[]
    for line in inlines:
        split=line.replace('\r\n','').split(';')
        indices_in.append(int(split[0]))
        labels_in.append(split[1])
        R_in.append(float(split[2])/255.)
        G_in.append(float(split[3])/255.)
        B_in.append(float(split[4])/255.)
    
    indices, labels, R, G, B = [],[],[],[],[]
    for ri in red_index:
        index = indices_in.index(ri)
        indices.append(indices_in[index])
        labels.append(labels_in[index])
        R.append(R_in[index])
        G.append(G_in[index])
        B.append(B_in[index])
    
    indices[-1]=250
    
    cmap = colors.ListedColormap([RGB for RGB in zip(R,G,B)])
    norm = colors.BoundaryNorm(indices,cmap.N)
    
    return indices, labels, cmap, norm


def CiFOR_1973_cmap():
    
    R = [0.1, 0.,  0.5, 0.,  0.5]
    G = [0.1, 0.5, 0.5, 0.3, 0.5]
    B = [0.1, 0.,  0.,  0.8, 0.5]
    indices = [0,1,2,3,4]
    labels = ['No Data', 'Forest', 'Non-forest', 'Water', 'Cloud' ]
    cmap = colors.ListedColormap([RGB for RGB in zip(R,G,B)])
    norm = colors.BoundaryNorm(indices,cmap.N)
    
    return indices, labels, cmap, norm
