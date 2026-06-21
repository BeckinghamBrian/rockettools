# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 18:25:27 2020

@author: brian

Code and method from https://www.benjaminmunro.com/moc, modified to plot 
nozzle out to design exit radius
"""

import numpy as np
from scipy import optimize

# calculating the Prandtl-Meyer function for a given number
def prandtlMeyer(M, gam = 1.4):
    return np.sqrt((gam+1)/(gam-1))*np.arctan(np.sqrt((gam-1)/(gam+1)*(M**2 - 1))) - np.arctan(np.sqrt(M**2 - 1))

# finds the mach number from the Prandtl-Meyer function
def machPM(M, nu, gam=1.4):
    return (np.sqrt((gam+1)/(gam-1))*np.arctan(np.sqrt((gam-1)/(gam+1)*(M**2 - 1))) - np.arctan(np.sqrt(M**2 - 1)))-nu

def nozzle(exitMach, n, Rt, Re):

    thetaInital = np.deg2rad(0.865) # angle of the first characteristic
    gam = 1.4 # ratio of specific heats
    x0 = 0 # location of throat and point where MOC begins
    y0 = Rt/1000 # location of throat and point where MOC begins, throat radius in mm
    Re = Re/1000 # back to m
    
    # Initialize arrays
    nu = np.zeros([n,n])  # Prandtl-Meyer function
    mu = np.zeros([n,n])  # Mach angle
    theta = np.zeros([n,n]) # Angle
    M = np.zeros([n,n])  # Mach number
    Kp = np.zeros([n,n])  # 
    Kn = np.zeros([n,n])  #
    x = np.zeros([n,n])  # x location of intersection points
    y = np.zeros([n,n])  # y location of intersection points
    xWall = np.zeros([n+1])  # Wall x location
    yWall = np.zeros([n+1])  # Wall y location
    
    
    
    # Find maximum angle of expansion corner
    thetaMax = prandtlMeyer(exitMach, gam)/2
    thetaDel = (thetaMax-thetaInital)/(n-1)
    
    for j in range(n):
            theta[0,j] = thetaInital +j*thetaDel
            nu[0,j] = theta[0,j]
            M[0,j] = optimize.newton(machPM, 1.1, args=(nu[0,j], gam))
            mu[0,j] = np.arcsin(1/M[0,j])
            Kp[0,j] = theta[0,j] -  nu[0,j]
            Kn[0,j] = theta[0,j]+nu[0,j]
            
    for j in range(n):
        for i in range(1,n-j):  
                Kp[i,j] = -Kn[0,i]
                theta[i,j] = j*thetaDel
                nu[i,j] =  theta[i,j] - Kp[i,j]
                M[i,j] = optimize.newton(machPM, 2, args=(nu[i,j], gam))
                mu[i,j] = np.arcsin(1/M[i,j])
                Kn[i,j] = theta[i,j]+nu[i,j]
                
    # need to calulate the first characteristic first. 
    for j in range(n):
        if j == 0:
            x[0,j] = x0 - y0/np.tan(theta[0,j]-mu[0,j])
            y[0,j] = 0
        else:
            alpha = np.tan(theta[0,j]-mu[0,j])
            beta = np.tan((theta[0,j-1] + mu[0,j-1] + theta[0,j] + mu[0,j])/2)
            x[0,j] = (y0-alpha*x0 - y[0,j-1]+beta*x[0,j-1])/(beta-alpha)
            y[0,j] = beta*(x[0,j]-x[0,j-1]) + y[0,j-1]
             
            #change out xo and yo vals
            
    # now to calulate the locations of the rest of the characteristic net. 
    for i in range(1,n):
        for j in range(0,n-i):
            if j == 0:
                x[i,j] = x[i-1,j+1] - y[i-1,j+1]/np.tan(theta[i,j]-mu[i,j])
                y[i,j] = 0
            else:
                alpha = np.tan((theta[i-1,j+1]-mu[i-1,j+1]+theta[i,j]-mu[i,j])/2)
                beta = np.tan((theta[i,j-1] + mu[i,j-1] + theta[i,j] + mu[i,j])/2)
                x[i,j] = (y[i-1,j+1]-alpha*x[i-1,j+1] - y[i,j-1]+beta*x[i,j-1])/(beta-alpha)
                y[i,j] = beta*(x[i,j]-x[i,j-1]) + y[i,j-1]
                
    # Calculate Nozzle wall
    xWall[0] = x0
    yWall[0] = y0
    for i in range(0,n):
        if i == 0:
            thetaWall = np.tan(thetaMax)
            beta = np.tan(theta[0,-1]+mu[0,-1])
            xWall[1] = (yWall[0]-thetaWall*xWall[0] - y[0,-1]+beta*x[0,-1])/(beta-thetaWall)
            yWall[1] = thetaWall*(xWall[1]-xWall[0]) + yWall[0]  
        else:
            thetaWall = np.tan((theta[i-1,-i]+theta[i,-(i+1)])/2)
            beta = np.tan(theta[i,-(i+1)]+mu[i,-(i+1)])
            xWall[i+1] = (yWall[i]-thetaWall*xWall[i] - y[i,-(i+1)]+beta*x[i,-(i+1)])/(beta-thetaWall)
            yWall[i+1] = thetaWall*(xWall[i+1]-xWall[i]) + yWall[i]
            
        # # This matches desired exit angle, but not expansion ratio
        # u = (yWall[i] - yWall[i-1]) / (xWall[i] - xWall[i-1])
        
        # # When the nozzle gets to the exit angle, remove the rest to the right
        # if 8 >= math.degrees(math.atan(u)):
        #     yWall = yWall[:len(yWall)-(n-i)]
        #     xWall = xWall[:len(xWall)-(n-i)]
        #     # z replaces n for the remaining calcs
        #     z = len(xWall)-1
        #     break
        
        # When the nozzle gets to the exit radius, remove the rest to the right
        if yWall[i] >= Re:
            yWall = yWall[:len(yWall)-(n-i)]
            xWall = xWall[:len(xWall)-(n-i)]
            # z replaces n for the remaining calcs
            z = len(xWall)-1
            break
        
    xChar = np.zeros([z,z+2]) 
    yChar = np.zeros([z,z+2]) 
    
    xChar[:,0] = x0
    xChar[:,-1] = xWall[1:]
    yChar[:,0] = y0
    yChar[:,-1] = yWall[1:]
    
    for j in range(z):
        xChar[j, j+2:-1] = x[j, 1:z-j]
        yChar[j, j+2:-1] = y[j, 1:z-j]
        w = 0
        for i in range(j, -1, -1):
            xChar[j, w+1] = x[w,i]
            yChar[j, w+1] = y[w,i]
            w += 1

    step = 0.00001 # step size for nozzle geo (100th of a mm)
    # fitting the wall coordinates to a 4th degree polynomial to space 
    # points equidistant apart horizontally
    coefs = np.polyfit(xWall,yWall,4)
    poly = np.poly1d(coefs)
    xContour = np.arange(xWall[0], xWall[-1], step)
    yContour = poly(xContour)
            
    xContour = np.delete(xContour,0)
    yContour = np.delete(yContour,0)
    
    # theta = ((yContour[-1] - yContour[-2])/(xContour[-1] - xContour[-2]))
    # print(theta)
    # tan = math.degrees(math.atan(theta))
    # print(tan)

    return xContour*1000, yContour*1000 # returns the x and y points in mm