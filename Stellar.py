# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:09:31 2018

@author: Grace
"""
#imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#Start of code
#Inputs solar mass, solar luminosity, effective temperature, H fraction, metal fraction
def star():
    
    #define variables
    #ALL UNITS IN CGS
    Nstarti = 10 #number of steps for starting equations
    Nstop = 999 #max number of steps
    Igoof = -1 #final model condition flag
    ierr = 0 
    P0 = 0 #surface pressure
    T0 = 0 #surface temp
    dlPlim = 99.9 
    Rsun = 6.9599E10 #radius of sun
    Msun = 1.989E33 #mass of sun
    Lsun = 3.826E33 #sun luminosity
    sigma = 5.67E-5 #stefan boltzmann
    c = 3E8 #speed of light
    a = (4*sigma)/c #rad pressure constant
    G = 6.67259E-8 #grav constant
    k_B = 1.380658E-16 #boltzman constants
    m_H = 1.673534E-24 #mass hydrogen
    pi = np.pi
    gamma=5/3
    gamrat=gamma/(gamma-1)
    tog_bf=0.01 #bound free opacity constant
    g_ff=1.0 #free-free opacity gaunt factor
    
    #Next open up a file but we don't have one: should ask about it
    
    Msolar=input("Enter star mass (solar units):")
    Lsolar=input("Enter star luminosity (solar units):")
    Teff=input("Enter effective temperature:")

    ###FIGURE THIS OUT###
    #user inputs mass fraction values, runs through loop to calculate Y
    X= input("Enter the mass fraction of hydrogen:")
    Z=input("Enter the mass fraction of metals:")
    Y = 1-X-Z #helium mass fraction
    while Y < 0:
        print("\n")
        X= input("Enter the mass fraction of hydrogen:")
        Z=input("Enter the mass fraction of metals:")
        Y = 1-X-Z #helium mass fraction
    
    XCNO = Z/2 #mass frac of CNO = 50% fo Z
    #Star mass, luminosity, radius
    Ms = Msolar*Msun
    Ls = Lsolar*Lsun
    Rs = np.sqrt(Ls/(4*pi*sigma))/(Teff**2)
    Rsolar = Rs/Rsun
    
    #start with small step size
    delstar = -Rs/1000
    #idrflg = size flag. 0 is initial surface step size
    idrflg = 0
    
    #mean molecular weight for complete ionization
    mu = 1/((2*X)+(0.75*Y)+(0.5*Z))
    
    #delimiter between adiabatic convection and radiation
    gamrat = gamma/(gamma-1)
    
    #initialize values at the surface. Outermost zone = 1
    r1 = Rs
    M_r1 = Ms
    L_r1 = Ls
    T1 = T0
    P1 = P0
    if P0 <= 0 or T0 <= 0:
        rho1 = 0
        kappa1 = 0
        epsilon1 = 0
    else:
        #Is a subroutine necessary for python?
        if ierr != 0:
            break #Is this right?
        
        
        
        
