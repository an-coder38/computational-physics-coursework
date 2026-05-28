# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:41:58 2026

@author: jfcte
"""

import pylab as py
from numpy import linspace, sin, exp, pi

# Data
x1 = linspace(0.0, 2.0, 20)
x2 = linspace(0.0, 2.0, 200)

y1 = exp(-x1)
y2 = exp(-x2)
y3 = sin(2 * pi * x2)
y4 = y2 * y3

# Graphs
l1 = py.plot(x1, y1, "bD-",
             markeredgecolor="black", # adds a black border to the graphed figures
             markeredgewidth=0.5) # width of the borders
l3 = py.plot(x2, y3, "go-", markeredgecolor="black", markeredgewidth=0.5)
l4 = py.plot(x2, y4, "rs-", markeredgecolor="black", markeredgewidth=0.5)

# Labeling
py.ylim(-1.1, 1.1)
py.xlabel("Segundos")
py.ylabel("Voltios")

py.legend((l3[0], l4[0]),
          ("Oscilatorio", "Amortiguado"),
          shadow=True)

py.title("Movimiento Oscilatorio Amortiguado")
py.show()
