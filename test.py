#!/bin/env python2.7

import numpy as np
a=np.arange(10)

names1=['a','b','c']
names2=['1','2']


DICT={ name1: { name2: {} for name2 in names2} for name1 in names1 }

for name1 in names1: 
    for name2 in names2: 
        DICT[name1][name2]=a

print(DICT)


