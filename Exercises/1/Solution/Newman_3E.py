# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:30 2026

@author: jfcte
"""

from pylab import imshow, show, jet, colorbar
from numpy import loadtxt
datos = loadtxt("circular.txt", float)
imshow(datos), jet(), colorbar()
show()