# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:41:26 2018

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt
import math

#dont worry about this. 
'''
    for i in range(0,931):
            file1 = str(int(i)).zfill(5)
            #loop through stellar structure files
            filepath = 'C:/Users/Grace/Documents/AnacondaProjects/stellar/structure_'+file1+'.txt'
            exists = os.path.isfile(filepath)
            if exists:
                with open(filepath, 'r') as f:
                    f1 = f.readlines()
                    for i in range(len(f1)):
                        li = f1[i].strip() #Removes whitespace from string ends
                        if not li.startswith('#'):  #Ignore all of the column headers
                            line_array1 = li.split() #Breaks the strings into elements
                            stellar_file.append(line_array1)
    #end for loop
'''
def stellar_star(file):

    stellar_file = []

    #enter in the path name of the structure files    
    with open('C:/Users/Grace/Documents/AnacondaProjects/stellar/structure_'+file+'.txt', 'r') as f1:
       classif = f1.readlines()
       for i1 in range(len(classif)):
           li1 = classif[i1].strip() #Removes whitespace from string ends
           if not li1.startswith('#'):  #Ignore all of the column headers
               line_array = li1.split() #Breaks the strings into elements
               stellar_file.append(line_array)
        
    #create whatever types of lists that you need to store your data (ex: mass, radius, etc)    
    mass = []
    radius = []
    temp = []  
    ad = []  
    rad = []     
    
    #Just append the data to the list. each column should be the same length. Be careful about indexing               
    for i0 in range(len(stellar_file)):
        mass.append(float(stellar_file[i0][0]))
        radius.append(float(stellar_file[i0][1]))
        temp.append(float(stellar_file[i0][5]))
        ad.append(float(stellar_file[i0][10]))
        rad.append(float(stellar_file[i0][15]))
    
    #convert to log scale if needed
    rad_log = []
    ad_log = []
    for i1 in range(len(rad)):
        rad_log.append(math.log(rad[i1]))
        ad_log.append(math.log(ad[i1]))
        
    #return all lists so you can plot them below
    return(mass, radius, temp, ad_log, rad_log)

#the wacky formatting entered as the parameter just creates the leading zeros that the files are labeled with
#int(NUMBER OF FILE).zfill(NUMBER OF ZEROS) (number of zeros should not change)
m, r, t, a, rad = stellar_star(file=str(int(0)).zfill(5))
m1, r1, t1, a1, rad1 = stellar_star(file=str(int(66)).zfill(5))
m2, r2, t2, a2, rad2 = stellar_star(file=str(int(128)).zfill(5))
m3, r3, t3, a3, rad3 = stellar_star(file=str(int(65)).zfill(5))

#plot statements. To plot on the same graph, get rid of the plt.figure() line
plt.plot(r,m, label='0 Gyr')
plt.plot(r1,m1, label='4.7 Gyr')
plt.plot(r2,m2, label='10 Gyr')
plt.xlabel("Radius")
plt.ylabel("Mass")
plt.title("Mass vs. Radius through Time")
plt.legend()

plt.figure()
plt.plot(t,m, label='0 Gyr')
plt.plot(t1,m1, label='4.7 Gyr')
plt.plot(t2,m2, label='10 Gyr')
plt.xlabel("Temperature")
plt.ylabel("Mass")
plt.title("Mass vs. Temperature through Time")
plt.legend()

plt.figure()
plt.plot(r3,a3, label='Adiabatic')
plt.plot(r3,rad3, label='Radiative')
plt.xlabel("Radius")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of the Present Day Sun")
plt.legend()

plt.figure()
plt.plot(m3,a3, label='Adiabatic')
plt.plot(m3,rad3, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of the Present Day Sun")
plt.legend()