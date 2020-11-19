# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 17:18:03 2020

@author: tropitac
"""

import numpy as np


def train_motion(t,y,params):
    #initialize train parameters
    Ls=params[0]    #stroke length of piston
    rw=params[1]    #wheel radius
    rg=params[2]    #gear radius
    Pistonradius=params[3]  #piston radius
    g=params[4]     #acceleration due to gravity
    mw=params[5]    #mass of a wheel
    pguage=params[6]    #initial guage pressure
    rho=params[7]   #density of air
    m=params[8]     #mass of train
    A=params[9]     #cross-sectional area of train
    Mus=params[10]  
    Cd=params[11]   #coefficient of drag
    Crr=params[12]  #coefficient of friction for rolling resistance
    Patm=params[13] #atmospheric pressure
    tlen=params[14] #length of train tank
    
    #train formulas
    Ap=Pistonradius**2*np.pi     #area of piston
    
    Frr=m*g*Crr         #rolling resistance force
    
    La=Ls*rw/rg         #Length  of acceleration
    
    v0=(tlen-Ls)*Ap
    
    #check for acceleration or deceleration
    if y[0]<La:
        dxdt=y[1]
        Fd=rho*Cd*A*dxdt**2/2
        Fteq=((rg*Ap)/rw)*(((pguage*v0)/(v0+Ap*(rg/rw)*y[0]))-Patm)
        dvdt=(Fteq-Fd-Frr)/(m+mw)
        Ft=((rg*Ap)/rw)*((pguage*v0)/(v0+Ap*(rg/rw)*y[0]))-mw*dvdt
        
        #slip check
        if Ft>Mus*(m/2)*g:
            print('Ft exceedingly high, slip occurs')
            raise ValueError
            
    #deceleration term               
    else:
        dxdt=y[1]
        Fd=rho*Cd*A*dxdt**2/2
        dvdt=(-Fd-Frr)/m
    
    dydt=np.array((dxdt,dvdt))
    
    return dydt









