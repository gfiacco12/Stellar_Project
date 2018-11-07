# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:09:31 2018

@author: Grace
"""
#imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc        
        
class Model:
    
    #these are class variables, shared by all instances
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
    kpad = 0.3 #adiabatic gas law constant
        
    #constuctor, called when a new instance of a class is created
    def __init__(self):
        #these are instance variables, different for each instance
        self.exampleB = 2
        self.m_H = 1.673534E-24 #mass hydrogen

        star()
        return
    
    #makes sure all terminal inputs are converted to integers
    def getInt(s):
        return(int(input(s)))
    
    #EOS sub-function 
    def EOS(X, Z, XCNO, mu, P, T, rho, kappa, epslon, tog_bf, izone, ierr):
        oneo3 = 0.33333
        twoo3 = 0.66667
        
        if T<0 or P<0:
            Model.ierr = 1
            print("Negative T or P, redo initial conditions. Problem in zone "+izone+" with conditions:"+T+'/n'+P)
        
        Prad = Model.a * T ** (4/3)
        Pgas = P - Prad
        rho = ((mu * Model.m_H) / Model.k_B) * (Pgas / T)
        
        if rho < 0:
            Model.ierr = 1
            print("Negative density, probably too high radiation pressure. Problem in zone "+izone+" with conditions:"+T+'/n'+P+'/n'+Prad+'/n'+Pgas+'/n'+rho)
        return()

    #Starmodel sub-function
    def Starmodel(deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, r, P_ip1, M_rip1, L_rip1, T_ip1, tog_bf, irc): 
        ###WHY RENAME THESE??###
        r = r_i + deltar
        M_rip1 = M_ri
        L_rip1 = L_ri
        
        #Radiative approximation
        if irc == 0:
            T_ip1 = ((Model.G * M_rip1 * Model.mu * Model.m_H) / (4.25 * Model.k_B)) * ((1/r) - (1/Rs))
            A_bf = 4.34 * Z * (1 + X) / tog_bf
            A_ff = 3.68 * Model.g_ff * (1 - Z) * (1 + X)
            Afac = A_bf + A_ff
            P_ip1 = np.sqrt((1/4.25) * (16 / (3 * Model.pi * Model.a * Model.c)) * ((Model.G * M_rip1) / L_rip1) * (Model.k_B / (Afac * Model.mu * Model.m_H))) * T_ip1**4.25
        else:
            #convective approximation
            T_ip1 = ((Model.G * M_rip1 * Model.mu * Model.m_H) / Model.k_B) * ((1/r) - (1/Rs)) / Model.gamrat
            P_ip1 = Model.kpad * T_ip1**Model.gamrat
        return()
    
    #Beginning of actual stellar model function
    def star():        
        #Next open up a file but we don't have one: should ask about it        
        Msolar = Model.getInt("Enter star mass (solar units):")
        Lsolar = Model.getInt("Enter star luminosity (solar units):")
        Teff = Model.getInt("Enter effective temperature:")
    
        #user inputs mass fraction values, runs through loop to calculate Y
        while True:
            X = Model.getInt("Enter the mass fraction of hydrogen:")
            Z = Model.getInt("Enter the mass fraction of metals:")
            Y = 1-X-Z #helium mass fraction
            if Y < 0:
                break
        
        XCNO = Z/2 #mass frac of CNO = 50% fo Z
        #Star mass, luminosity, radius
        Ms = Msolar * Model.Msun
        Ls = Lsolar * Model.Lsun
        Rs = np.sqrt(Ls / (4 * Model.pi * Model.sigma)) / (Teff**2)
        Rsolar = Rs / Model.Rsun
        
        #start with small step size
        delstar = -Rs / 1000
        #idrflg = size flag. 0 is initial surface step size
        idrflg = 0
        
        #mean molecular weight for complete ionization
        mu = 1 / ((2 * X) + (0.75 * Y) + (0.5 * Z))
          
        r = [Rs]
        M_r = [Ms]
        L_r = [Ls]
        T = [Model.T0]
        P = [Model.P0]
        rho = []
        kappa = []
        epslon = []
        if Model.P0 <= 0 or Model.T0 <= 0:
            rho = 0
            kappa = 0
            epslon = 0
        else:
            Model.EOS(X, Z, XCNO, mu, P, T, rho, kappa, epslon, Model.tog_bf, 1, Model.ierr)
        
        #Apply surface conditions to begin integration. Radiation transport in outermost zone: irc = 0.
        #Arbitrary initial values for kPad and dlPdlT
        ##dlPdlT = dlnP/dlnT##
        
        irc = 0 
        dlPdlT = [4.25]
        
        ###Learn about DO LOOPS###
        
        return()
    
model = Model()