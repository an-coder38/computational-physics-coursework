# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:30 2026

@author: jfcte
"""
# Write a program that prints in increasing order all Catalan numbers less than or equal to one
# billion.


c=1 #first catalan number
n=0
lim = 1*10**9 #the limit is set to one billion = one thousand million, if it were the spanish billion it would be 10**12
while c<lim: #notice that the last calculated catalan number may be calculated but not printed
    print(int(c)) #prints c as an integer as the catalan numbers are integers
    c=(4*n+2)*c/(n+2) #defining c(n+1) with c(n)
    n+=1

#limiting at 10**9, after the simulation the numbers are:
# 1
# 1
# 2
# 5
# 14
# 42
# 132
# 429
# 1430
# 4862
# 16796
# 58786
# 208012
# 742900
# 2674440
# 9694845
# 35357670
# 129644790
# 477638700
#the next number would be 1767263190 (bigger than 10**9)

# limiting at 10**12 after the simulation the numbers are:
# 1
# 1
# 2
# 5
# 14
# 42
# 132
# 429
# 1430
# 4862
# 16796
# 58786
# 208012
# 742900
# 2674440
# 9694845
# 35357670
# 129644790
# 477638700
# 1767263190
# 6564120420
# 24466267020
# 91482563640
# 343059613650