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
def star(Msolar, Lsolar, Teff, X, Z):
    
    #define variables
    #ALL UNITS IN CGS
    Nstart=10 #number of steps for starting equations
    Nstop=999 #max number of steps
    Igoof=-1 #final model condition flag
    ierr=0 
    P0=0 #surface pressure
    T0=0 #surface temp
    dlPlim=99.9 
    Rsun=6.9599E10 #radius of sun
    Msun=1.989E33 #mass of sun
    Lsun=3.826E33 #sun luminosity
    sigma= 5.67E-5 #stefan boltzmann
    c=3E8 #speed of light
    a=(4*sigma)/c #rad pressure constant