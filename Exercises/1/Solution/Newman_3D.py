# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:21 2026

@author: jfcte
"""

from pylab import scatter, xlabel, ylabel, xlim, ylim, title, show, text
from numpy import loadtxt
datos = loadtxt("stars.txt", float)
x, y = datos[: , 0], datos[: , 1]
scatter(x, y, c = x, cmap = "hsv", s = 75,
        edgecolors="black",   # black border of the circles graphed
        linewidths=0.5)        # width of the border
# both lines commented are needed so the graph obtained with the code is the 
# one in the slides "Tema1"

xlabel("Temperatura")
ylabel("Magnitud")
xlim(13000, 0)
ylim(20, -5)
text(4500, 4.5, "Secuencia Principal", fontsize = 11, style = "italic")
text(11000, 15, "Enanas Blancas", fontsize = 11, style = "italic")
title("Diagrama Hertzsprung-Russell")
show()