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
def stellar_star(star, file):

    stellar_file = []

    #enter in the path name of the structure files    
    with open('C:/Users/Grace/Documents/AnacondaProjects/'+star+'/structure_'+file+'.txt', 'r') as f1:
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
    conv_vel = []
    
    #Just append the data to the list. each column should be the same length. Be careful about indexing               
    for i0 in range(len(stellar_file)):
        mass.append(float(stellar_file[i0][0]))
        radius.append(float(stellar_file[i0][1]))
        temp.append(float(stellar_file[i0][5]))
        ad.append(float(stellar_file[i0][10]))
        rad.append(float(stellar_file[i0][15]))
        conv_vel.append(float(stellar_file[i0][17]))
    
    #convert to log scale if needed
    rad_log = []
    ad_log = []
    for i1 in range(len(rad)):
        rad_log.append(math.log(rad[i1]))
        ad_log.append(math.log(ad[i1]))
        
    #return all lists so you can plot them below
    return(mass, radius, temp, ad_log, rad_log, conv_vel)

#the wacky formatting entered as the parameter just creates the leading zeros that the files are labeled with
#int(NUMBER OF FILE).zfill(NUMBER OF ZEROS) (number of zeros should not change)
m, r, t, a, rad, conv_vel = stellar_star(star='10M Z0.0001', file=str(int(0)).zfill(5))
m1, r1, t1, a1, rad1, conv_vel1 = stellar_star(star='10M Z0.0001', file=str(int(20)).zfill(5))
m2, r2, t2, a2, rad2, conv_vel2 = stellar_star(star='10M Z0.0001', file=str(int(40)).zfill(5))
m3, r3, t3, a3, rad3, conv_vel3 = stellar_star(star='10M Z0.0001', file=str(int(60)).zfill(5))
#m4, r4, t4, a4, rad4, conv_vel4 = stellar_star(star='10M Z0.0001', file=str(int(800)).zfill(5))
#m5, r5, t5, a5, rad5 = stellar_star(star='100M Z0.0001', file=str(int(750)).zfill(5))
#m6, r6, t6, a6, rad6 = stellar_star(star='1M Z0.004', file=str(int(850)).zfill(5))
#m7, r7, t7, a7, rad7 = stellar_star(star='1M Z0.004', file=str(int(1000)).zfill(5))


#plot statements. To plot on the same graph, get rid of the plt.figure() line

plt.figure()
plt.plot(r,conv_vel, label='Adiabatic')
#plt.plot(m,rad, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of 100M Star: 0 Yrs")
plt.legend()

plt.figure()
plt.plot(r1,conv_vel1, label='Adiabatic')
#plt.plot(m1,rad1, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of 100M Star: 2.97x10$^{6}$ yrs")
plt.legend()

plt.figure()
plt.plot(r2,conv_vel2, label='Adiabatic')
#plt.plot(m2,rad2, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of 100M Star: 3.01x10$^{6}$ yrs")
plt.legend()

plt.figure()
plt.plot(r3,conv_vel3, label='Adiabatic')
#plt.plot(m3,rad3, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of 100M Star: 3.05x10$^{6}$ yrs")
plt.legend()

plt.figure()
plt.plot(r4,conv_vel4, label='Adiabatic')
#plt.plot(m4,rad4, label='Radiative')
plt.xlabel("Mass")
plt.ylabel("Temperature Gradient")
plt.title("Temperature Gradients of 100M Star: 3.23x10$^{6}$ yrs")
plt.legend()

#plt.figure()
#plt.plot(r5,a5, label='Adiabatic')
#plt.plot(r5,rad5, label='Radiative')
#plt.xlabel("Radius")
#plt.ylabel("Temperature Gradient")
#plt.title("Temperature Gradients of 100M Star: 2.63x10$^{7}$ yrs")
#plt.legend()

#plt.figure()
#plt.plot(m6,a6, label='Adiabatic')
#plt.plot(m6,rad6, label='Radiative')
#plt.xlabel("Mass")
#plt.ylabel("Temperature Gradient")
#plt.title("Temperature Gradients of 1M Star: 2.67x10$^{7}$ yrs")
#plt.legend()

#plt.figure()
#plt.plot(m7,a7, label='Adiabatic')
#plt.plot(m7,rad7, label='Radiative')
#plt.xlabel("Mass")
#plt.ylabel("Temperature Gradient")
#plt.title("Temperature Gradients of 1M Star: 2.7x10$^{7}$ yrs")
#plt.legend()