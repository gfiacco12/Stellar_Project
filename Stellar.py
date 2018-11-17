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
    Nstart = 10 #number of steps for starting equations
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
    c = 3E10 #speed of light
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
       
        #calculate opacity
        Model.tog_bf = 2.8 * (rho * (1 + X)**0.20)
        ##Double check this with the book##
        #bound-free, free-free, e scattering opacities
        k_bf = 4.34 / Model.tog_bf * Z * (1 + X) * rho /(T**3.5)
        k_ff = 3.68 * Model.g_ff * (1 - Z)*(1 + X)*rho/(T**3.5)
        k_e = 0.2 * (1 + X)
        kappa = k_bf + k_ff + k_e
        
        #energy generation from pp chain and CNO
        T6 = T * 1E-06
        fx = 0.133 * X * np.sqrt((3 + X) * rho) / T6**1.5
        fpp = 1 + fx*X
        psipp = 1 + 1.412E8 * (1/X - 1)*np.exp(-49.98 * T6**(-oneo3))
        Cpp = 1 + 0.0123 * (T6**oneo3) + 0.0109 * (T6**twoo3) + 0.000938 * T6
        epspp = 2.38E6 * rho * X**2 * fpp * psipp * Cpp * T6**(-twoo3) * np.exp(-33.80 * T6**(-oneo3))
        CCNO = 1 + 0.0027 * T6**(-oneo3) - 0.00778 * T6**(-twoo3) - 0.000149 * T6
        epsCNO = 8.67E27 * rho * X * XCNO * CCNO * T6**(-twoo3) * np.exp(-152.28 * T6**(-oneo3))
        ##HELP
        epslon = epspp + epsCNO
        
        return()

    #Starmodel sub-function
    def Starmodel(deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, r, P_ip1, M_rip1, L_rip1, T_ip1, tog_bf, irc): 

        r = r_i + Model.deltar
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
    
    #Functions to calculate pressure, mass, luminosity, an temp gradients 
    def dPdr(r, M_r, rho):
        dPdr = -Model.G * rho * M_r / r**2 
        return()
    
    def dMdr(r, rho):
        dMdr = 4 * Model.pi * rho * r**2
        return()
    
    def dLdr(r, rho, epslon):
        dLdr = 4 * Model.pi * rho * epslon * r**2
        return()
    
    def dTdr(r, M_r, L_r, T, rho, kappa, mu, irc):
        #radiative temp gradient
        if irc == 0:
            dTdr = -((3 / (16 * Model.pi * Model.a * Model.c)) * kappa * (rho / (T**3)) * (L_r / r**2))
        else:
            #adiabatic convective gradient
            dTdr = -(1 / Model.gamrat) * (Model.G * M_r /r**2) * (mu * Model.m_H / Model.k_B)
        return()
    
        #FUNDEQ sub-function 
    def FUNDEQ(r, f, dfdr, irc, X, Z, XCNO, mu, izone, ierr):
        
        P = f[0]
        M_r = f[1]
        L_r = f[2]
        T = f[3]
        
        Model.EOS(X, Z, XCNO, mu, P, T, Model.rho, Model.kappa, Model.epslon, Model.tog_bf, izone, ierr)
        
        dfdr[0] = Model.dPdr(r, M_r, Model.rho)
        dfdr[1] = Model.dMdr(r, Model.rho)
        dfdr[2] = Model.dLdr(r, Model.rho, Model.epslon)
        dfdr[3] = Model.dTdr(r, M_r, L_r, T, Model.rho, Model.kappa, mu, irc)
        
        return()
        
    #RUNGE sub-function 
    def RUNGE(f_im1, dfdr, f_i, r_im1, deltar, irc, X, Z, XCNO, mu, izone, ierr):
        
        dr12 = deltar / 2
        dr16 = deltar / 6
        r12 = r_im1 + dr12
        r_i = r_im1 +deltar
        
        f_temp = []
        df1 = []
        df2 = []
        df3 = []
        f_i = []
        
        for i in range(0,4):
            f_temp[i] = f_im1[i] + dr12 * dfdr[i]
            
        Model.FUNDEQ(r12, f_temp, df1, irc, X, Z, XCNO, mu, izone, ierr)
        
        if ierr != 0:
            return()
                
        for i in range(0,4):
            f_temp[i] = f_im1[i] + dr12 * df1[i]
            
        Model.FUNDEQ(r12, f_temp, df2, irc, X, Z, XCNO, mu, izone, ierr)
        
        if ierr != 0:
            return()
            
        for i in range(0,4):
            f_temp[i] = f_im1[i] + deltar*df2[i]

        Model.FUNDEQ(r12, f_temp, df3, irc, X, Z, XCNO, mu, izone, ierr)
        if ierr != 0:
            return()
            
        for i in range(0,4):
            f_i[i] = f_im1[i] + dr16 * (dfdr[i] + 2 * df1[i] + 2 * df2[i] + df3[i])
            
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
                print('X+Z must be <= 1. Reenter composition.')
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
            rho[0] = 0
            kappa[0] = 0
            epslon[0] = 0
        else:
            Model.EOS(X, Z, XCNO, mu, P[0], T[0], rho[0], kappa[0], epslon[0], Model.tog_bf, 1, Model.ierr)
        
        #Apply surface conditions to begin integration. Radiation transport in outermost zone: irc = 0.
        #Arbitrary initial values for kPad and dlPdlT
        ##dlPdlT = dlnP/dlnT##
        
        irc = 0 
        dlPdlT = [4.25]
        
        for i in range(Model.Nstart):
            ip1 = i + 1
            Model.Starmodel(deltar, X, Z, mu, Rs, r[i], M_r[i], L_r[i], r[ip1], P[ip1], M_r[ip1], L_r[ip1], T[ip1], Model.tog_bf, irc)
            Model.EOS(X, Z, XCNO, mu, P[ip1], T[ip1], rho[ip1], kappa[ip1], epslon[ip1], Model.tog_bf, ip1, Model.ierr)
        
            if Model.ierr != 0:
                print('Values from the previous zone are:', r[i]/Rs, rho[1], M_r[i]/Ms, kappa[i], T[i], epslon[i], P[i], L_r[i]/Ls)

            #Determine whether conduction happens in next zone
            if i > 1:
                dlPdlT[ip1] = np.log10(P[ip1]/P[i]) / np.log10(T[ip1]/T[i])
            else:
               dlPdlT[ip1] = dlPdlT[i] 
            
            if dlPdlT[ip1] < Model.gamrat:
                irc = 1
            else:
                irc = 0
                Model.kpad = (P[ip1]/T[ip1])**Model.gamrat
            
            #Check surface assumption of constant mass
            deltaM = deltar * Model.dMdr(r[ip1], rho[ip1])
            M_r[ip1] = M_r[i] + deltaM
            
            if np.abs(deltaM) > 0.001*Ms:
                print("The variation in mass has become larger than 0.001Ms. Leaving loop before Nstart was reached.")
                if ip1 > 2:
                    ip1 = ip1 - 1
                    break #exits the loop early
                else:
                    continue #continues with the loop
            
        #Main integration loop
        Nsrtp1 = ip1 + 1
        #start RK4 routine
        for i in range(Model.Nstart, Model.Nstop):
            im1 = i - 1
            f_im1 = [P[im1], M_r[im1], L_r[im1], T[im1]]
            dfdr = [Model.dPdr(r[im1], M_r[im1], rho[im1]), Model.dMdr(r[im1], rho[im1]), Model.dLdr(r[im1], rho[im1], epslon[im1]), Model.dTdr(r[im1], M_r[im1], L_r[im1], T[im1], rho[im1], kappa[im1], mu, irc)]
            
            Model.RUNGE(f_im1, dfdr, Model.f_i, r[im1], deltar, irc, X, Z, XCNO, mu, izone, Model.ierr)
           
            if Model.ierr != 0:
                print("The Problem occurred in the Runge-Kutta routine")
                print("Values from the previous zone are: ", r/Rs, rho, M_r/Ms, kappa, T, epslon, P, L_r/Ls)
                
            r[i] = r[im1] + deltar
            P[i] = Model.f_i[0]
            M_r[i] = Model.f_i[1]
            L_r[i] = Model.f_i[2]
            T[i] = Model.f_i[3]
            
            #calculate density, opacity, energy generation rate
            Model.EOS(X, Z, XCNO, mu, P, T, rho, kappa, epslon, Model.tog_bf, izone, Model.ierr)
            
            if Model.ierr != 0:
                print("Values from the previous zone are: ", r[im1]/Rs, rho[im1], M_r[im1]/Ms, kappa[im1], T[im1], epslon[im1], P[im1], L_r[im1]/Ls)
            
            #check if convection is operating in next zone
            dlPdlT[i] = np.log(P[i] / P[im1]) / np.log(T[i] / T[im1])
            
            if dlPdlT[i] < Model.gamrat:
                irc = 1
                
            else:
                irc = 0
                
            #checks for center
            if r[i] <= np.abs(deltar) and np.logical_or(L_r[i] >= 0.1*Ls, M_r >= 0.1*Ms):
                Model.Igoof = 6
                
            elif L_r[i] <= 0:
                Model.Igoof = 5
                rhocor = M_r[i] / (4 / (3 * Model.pi * r[i]**3))
                
                if M_r[i] != 0:
                    epscor = L_r[i] / M_r[i]
                    
                else:
                    epscor = 0
                    
                Pcore = P[i] + 2/(3 * Model.pi * Model.G * rhocor**2 * r[i]**2)
                Tcore = Pcore * mu * Model.m_H / (rhocor * Model.k_B)
                
            elif M_r[i] <= 0:
                Model.Igoof = 4
                rhocor = 0
                epscor = 0
                Pcore = 0
                Tcore = 0
                
            elif r[i] < 0.02*Rs and M_r[i] < 0.01*Ms and L_r[i] < 0.1*Ls:
                rhocor = M_r[i] / (4 / (3 * Model.pi * r[i]**3))
                rhomax = 10 * (rho[i] / rho[im1]) * rho[i]
                epscor = L_r[i] / M_r[i]
                Pcore = P[i] + 2 / (3 * Model.pi * Model.G *rhocor**2 * r[i]**2)
                Tcore = Pcore * mu * Model.m_H / (rhocor * Model.k_B)
                
                if  rhocor < rho[i] or rhocor > rhomax:
                    Model.Igoof = 1
                    
                elif epscor < epslon[i]:
                    Model.Igoof = 2
                
                elif Tcore < T[i]:
                    Model.Igoof = 3
                    
                else:
                    Model.Igoof = 0
            
            if Model.Igoof != -1:
                istop = i
            
        #Change step sizes 
            if idrflg == 0 and M_r[i] < 0.99*Ms:
                deltar = -Rs / 100
                idrflg = 1
            
            if idrflg == 1 and deltar >= 0.5*r[i]:
                deltar = -Rs / 5000
                idrflg = 2
            
            istop = i
            
        #generate warnings for central conditions:
        
            rhocor = 
            
        return()
    
model = Model()