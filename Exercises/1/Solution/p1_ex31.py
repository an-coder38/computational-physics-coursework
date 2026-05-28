# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:30 2026

@author: jfcte
"""
import numpy as np
import matplotlib.pyplot as plt

# In the on-line resources you will find a file called sunspots.txt, which contains the observed
# number of sunspots on the Sun for each month since January 1749. The file contains two
# columns of numbers, the first being the month and the second being the sunspot number.

# a) Write a program that reads in the data and makes a graph of sunspots as a function of time.

data = np.loadtxt("sunspots.txt") # loads the data from sunspots.txt
# time = data[:, 0] # first column contains the month
# sunspots = data[:, 1] # second column contains the sunspot number

# could be done at the same time using 
time, sunspots = data[: , 0], data[: , 1]
# time contains the month starting in Jan-1749
# sunspots contains the sunspot number for that month

# designing the graph
plt.figure()
plt.plot(time, sunspots)
plt.xlabel("Time (months since January 1749)")
plt.ylabel("Sunspot number")
plt.title("a) Sunspots as a function of time")
plt.show()


# b) Modify your program to display only the first 1000 data points on the graph.

N = 1000 #first 1000 points
time_1000 = time[:N] #only loads from 0 to N-1
sunspots_1000 = sunspots[:N] #only loads from 0 to N-1

plt.figure()
plt.plot(time_1000, sunspots_1000)
plt.xlabel("Time (months since January 1749)")
plt.ylabel("Sunspot number")
plt.title("b) 1000 first sunspots as a function of time")
plt.show()


# c) Modify your program further to calculate and plot the running average of the data, defined by --
# where r = 5 in this case (and the yk are the sunspot numbers). Have the program plot
# both the original data and the running average on the same graph, again over the range
# covered by the first 1000 data points.

l=len(time)
r = 5
running_avg = []

# Work out the mobile mean value
for k in range(0, l):
    if k<r:
        avg = np.mean(sunspots[0: k + r + 1]) #using np.mean saves to 
        running_avg.append(avg) #creates a list of the averages and adds each new result
    elif k>l-r:
        avg = np.mean(sunspots[k-r: l + 1]) #using np.mean saves to 
        running_avg.append(avg) #creates a list of the averages and adds each new result
    else:
        avg = np.mean(sunspots[k - r : k + r + 1])
        running_avg.append(avg)

#could be done with one for and if, elif, elif but the code runs fast and it is clearer this way
running_avg = np.array(running_avg)
# the whole plot
plt.figure()
plt.plot(time, sunspots, label="Sunspots", alpha=0.6)
plt.plot(time, running_avg, label="Average sunspots (r = 5)", linewidth=2)
plt.xlabel("Time (months since January 1749)")
plt.ylabel("Sunspot number")
plt.title("Sunspots as a function of time")
plt.legend()
plt.show()

# c) only 1000 data points
avg_1000 = running_avg[:N]
#change time to time_1000
plt.figure()
plt.plot(time_1000, sunspots_1000, label="Sunspots", alpha=0.6)
plt.plot(time_1000, avg_1000, label="Average sunspots (r = 5)", linewidth=2)
plt.xlabel("Time (months since January 1749)")
plt.ylabel("Sunspot number")
plt.title("1000 first sunspots as a function of time")
plt.legend()
plt.show()
