#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:41:33 2023

@author: tempchiao
"""

import numpy as np
import csv
import multipletau               #Install this packaged
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd

def get(path):
    geta = []
    getb = []
    readData = pd.read_csv(path, header = None, sep = '\\t', engine='python')
    a = readData[0]
    b = readData[1]
    for i in a:
        geta.append(i)
    for i in b:
        getb.append(i)
    return geta, getb


########################################### Run all scripts ############################################

k_green=-4.59832392481725  
w_green=2.4793999921704524e-07  
k_red=2.439674219185456 
w_red=2.886072421052287e-07

'''def runall():
    toloadsample = r"/Volumes/Tianxiao/20231117_Lysate/Beads/1in10_1in10_FCS" #C:/Users/Admin User/Desktop/Tianxiao/FCS/TDP-43/1in10/Zoe
    # This is the FCS file to load that has the data from the dye. This is to determine the beam waist. '''
                                        
                                             

########################################### Autocorrelate ############################################
    
def autocorrelate():
   global c
   global new_c
   global new_d
  
   file_path = '/Volumes/Elements/20240419_FCS/EV/'
   root = 'FCS'
   file_number = 3
   dataA = []
   dataB = []
   for i in range(1,file_number+1):
        
       if i == 1:
           path = file_path + root
           m,n = get(path)
           for p in m:
               dataA.append(p)
           for q in n:
                dataB.append(q)
            
       if i > 1:
           if i < 10:
                path = file_path + root + '_0' + str(i)
                m,n = get(path)
                for p in m:
                    dataA.append(p)
                for q in n:
                    dataB.append(q)
       
       if i > 9:
            path = file_path + root + '_' + str(i)
            m,n = get(path)
            for p in m:
                dataA.append(p)
            for q in n:
                dataB.append(q)        
     
        
   x=np.array(dataA,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   c=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_c=np.delete(c, (0), axis=0)                                                  # This deletes the first row which is at t = 0. 
   
       
   
   x=np.array(dataB,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   d=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_d=np.delete(d, (0), axis=0) 
   
   return file_path, root




def fungreen(x,n,td):
    k=k_green                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitgreen(path,root):
    xdata=new_c[:, 0]  
    ydata=new_c[:,1]
    guesses=np.array([20,6e-5])
    (n_, td_), _ = opt.curve_fit(fungreen, xdata, ydata,guesses)
    params= opt.curve_fit(fungreen, xdata, ydata,guesses)
    
   
    
    y_fit = fungreen(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "--k",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_c[:,1]))
    ax1.plot(xdata, y_fit, '-',color='green')
    
    
    print ("Green_N = %r \r") %params[0][0]
    print ("Green_td = %r \r") %params[0][1]

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print ("Green_D = %r \r") %Diff
    
    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print ("Green_r = %r \r") %Rh
    plt.savefig('{}{}_Green.png'.format(path, root)) 
    
def funred(x,n,td):
    k=k_red                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitred(path,root):
    xdata=new_d[:, 0]  
    ydata=new_d[:,1]
    guesses=np.array([25,5e-5])
    (n_, td_), _ = opt.curve_fit(funred, xdata, ydata,guesses)
    params= opt.curve_fit(funred, xdata, ydata,guesses)
    
   
    
    y_fit = funred(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "--k",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_d[:,1]))
    ax1.plot(xdata, y_fit, '-',color='red')
    
    
    print ("Red_N = %r \r") %params[0][0]
    print ("Red_td = %r \r") %params[0][1]

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print ("Red_D = %r \r") %Diff

    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print ("Red_r = %r \r") %Rh
    plt.savefig('{}{}_Green.png'.format(path, root)) 
    

p,r = autocorrelate()
fitgreen(p,r)
fitred(p,r)
    