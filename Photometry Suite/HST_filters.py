#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:05:17 2021

@author: brianmerino

This code opens all of the HST throughput data that my project uses so that my photometry code
can call on this code and overplot the data over my photometry.
"""

import matplotlib.pyplot as plt
import numpy as np

path_wfc3 = '/Users/brian.merino/Google Drive/My Drive/Clumpy_Research/Apperture_Test/HST_filt/UVIS/'
path_acs  = '/Users/brian.merino/Google Drive/My Drive/Clumpy_Research/Apperture_Test/HST_filt/ACS/'

#path_wfc3 = '/Users/brianmerino/Desktop/Apperture_Test/HST_filt/UVIS/'
#path_acs = '/Users/brianmerino/Desktop/Apperture_Test/HST_filt/ACS/'

# WFC3 filters
f105w   = 'F105W_IR_throughput.csv'
f110w   = 'F110W_IR_throughput.csv'
f125w   = 'F125W_IR_throughput.csv'
f140w   = 'F140W_IR_throughput.csv'
f160w   = 'F160W_IR_throughput.csv'
f225w   = 'F225W_UVIS_throughput.csv'
f275w   = 'F275W_UVIS_throughput.csv'
f336w   = 'F336W_UVIS_throughput.csv'
f390w   = 'F390W_UVIS_throughput.csv'

# ACS filters
f435w   = 'wfc_F435W.dat'
f475w   = 'wfc_F475W.dat'
f606w   = 'wfc_F606W.dat'
f625w   = 'wfc_F625W.dat'
f775w   = 'wfc_F775W.dat'
f814w   = 'wfc_F814W.dat'
f850lp  = 'wfc_F850LP.dat'


def HST_filt():
    x_f105w = np.loadtxt(path_wfc3+f105w, usecols=[1], skiprows=1, delimiter=',')
    y_f105w = np.loadtxt(path_wfc3+f105w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f110w = np.loadtxt(path_wfc3+f110w, usecols=[1], skiprows=1, delimiter=',')
    y_f110w = np.loadtxt(path_wfc3+f110w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f125w = np.loadtxt(path_wfc3+f125w, usecols=[1], skiprows=1, delimiter=',')
    y_f125w = np.loadtxt(path_wfc3+f125w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f140w = np.loadtxt(path_wfc3+f140w, usecols=[1], skiprows=1, delimiter=',')
    y_f140w = np.loadtxt(path_wfc3+f140w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f160w = np.loadtxt(path_wfc3+f160w, usecols=[1], skiprows=1, delimiter=',')
    y_f160w = np.loadtxt(path_wfc3+f160w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f225w = np.loadtxt(path_wfc3+f225w, usecols=[1], skiprows=1, delimiter=',')
    y_f225w = np.loadtxt(path_wfc3+f225w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f275w = np.loadtxt(path_wfc3+f275w, usecols=[1], skiprows=1, delimiter=',')
    y_f275w = np.loadtxt(path_wfc3+f275w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f336w = np.loadtxt(path_wfc3+f336w, usecols=[1], skiprows=1, delimiter=',')
    y_f336w = np.loadtxt(path_wfc3+f336w, usecols=[2], skiprows=1, delimiter=',')
    
    x_f390w = np.loadtxt(path_wfc3+f390w, usecols=[1], skiprows=1, delimiter=',')
    y_f390w = np.loadtxt(path_wfc3+f390w, usecols=[2], skiprows=1, delimiter=',')
    
    
    # The acs files don't need to skip rows because they don't have headers
    
    x_f435w = np.loadtxt(path_acs+f435w, usecols=[0], delimiter=' ')
    y_f435w = np.loadtxt(path_acs+f435w, usecols=[1], delimiter=' ')
    
    x_f475w = np.loadtxt(path_acs+f475w, usecols=[0], delimiter=' ')
    y_f475w = np.loadtxt(path_acs+f475w, usecols=[1], delimiter=' ')
    
    x_f606w = np.loadtxt(path_acs+f606w, usecols=[0], delimiter=' ')
    y_f606w = np.loadtxt(path_acs+f606w, usecols=[1], delimiter=' ')
    
    x_f625w = np.loadtxt(path_acs+f625w, usecols=[0], delimiter=' ')
    y_f625w = np.loadtxt(path_acs+f625w, usecols=[1], delimiter=' ')
    
    x_f775w = np.loadtxt(path_acs+f775w, usecols=[0], delimiter=' ')
    y_f775w = np.loadtxt(path_acs+f775w, usecols=[1], delimiter=' ')
    
    x_f814w = np.loadtxt(path_acs+f814w, usecols=[0], delimiter=' ')
    y_f814w = np.loadtxt(path_acs+f814w, usecols=[1], delimiter=' ')
    
    x_f850lp = np.loadtxt(path_acs+f850lp, usecols=[0], delimiter=' ')
    y_f850lp = np.loadtxt(path_acs+f850lp, usecols=[1], delimiter=' ')
    
    return x_f105w, y_f105w, x_f110w, y_f110w, x_f125w, y_f125w, x_f140w, y_f140w,\
           x_f160w, y_f160w, x_f225w, y_f225w, x_f275w, y_f275w, x_f336w, y_f336w,\
           x_f390w, y_f390w, x_f435w, y_f435w, x_f475w, y_f475w, x_f606w, y_f606w,\
           x_f625w, y_f625w, x_f775w, y_f775w, x_f814w, y_f814w, x_f850lp, y_f850lp


# To make sure that my code is working, I can uncomment this section to verify that the code is plotting correctly.
# Make sure to comment everything below this line when done checking code.
###################################################################################################
'''
x_f105w, y_f105w, x_f110w, y_f110w, x_f125w, y_f125w, x_f140w, y_f140w,\
x_f160w, y_f160w, x_f225w, y_f225w, x_f275w, y_f275w, x_f336w, y_f336w,\
x_f390w, y_f390w, x_f435w, y_f435w, x_f475w, y_f475w, x_f606w, y_f606w,\
x_f625w, y_f625w, x_f775w, y_f775w, x_f814w, y_f814w, x_f850lp, y_f850lp = HST_filt()

fig = plt.figure(figsize=(8,6))
p1 = fig.add_subplot(111)

p1.plot(x_f105w,  y_f105w,  label='F105W')
plt.plot(x_f110w,  y_f110w,  label='F110W')
plt.plot(x_f125w,  y_f125w,  label='F125W')
plt.plot(x_f140w,  y_f140w,  label='F140W')
plt.plot(x_f160w,  y_f160w,  label='F160W')
plt.plot(x_f225w,  y_f225w,  label='F225W')
plt.plot(x_f275w,  y_f275w,  label='F275W')
plt.plot(x_f336w,  y_f336w,  label='F336W')
plt.plot(x_f390w,  y_f390w,  label='F390W')
plt.plot(x_f475w,  y_f475w,  label='F475W')
plt.plot(x_f625w,  y_f625w,  label='F625W')
plt.plot(x_f775w,  y_f775w,  label='F775W')
plt.plot(x_f814w,  y_f814w,  label='F814W')
plt.plot(x_f850lp, y_f850lp, label='F850lp')

# This section will display approximately where the emission lines H-alpha
# and [OIII] should be after accounting for the source's redshift.
p3 = p1.twinx()
#z = 0.579

z = 0

H_alpha = (1+z)*6563  #Angstroms
H_beta  = (1+z)*4861.3  #Angstroms
OII     = (1+z)*3727  #Angstroms
OIII    = (1+z)*5007  #Angstroms

alpha1 = 0.5
p3.axvline(H_alpha, c='black', alpha=alpha1)
p3.axvline(OIII, c='black', alpha=alpha1)
p3.axvline(H_beta, c='black', alpha=alpha1)
p3.axvline(OII, c='black', alpha=alpha1)
plt.ylim(0,0.6)

p3.annotate('H-$\\alpha$', xy=(H_alpha, 0.56), xytext=(H_alpha, 0.57), rotation=90)
p3.annotate('H-$\\beta$', xy=(H_beta, 0.56), xytext=(H_beta-400, 0.57), rotation=90)
p3.annotate('[OII]', xy=(OII, 0.56), xytext=(OII, 0.565), rotation=90)
p3.annotate('[OIII]', xy=(OIII, 0.56), xytext=(OIII, 0.56), rotation=90)

p1.legend(bbox_to_anchor=(1, 0.09))
p1.set_xlim(1000,18000)
p1.set_ylim(0,0.6)
p3.set_ylim(0,0.6)
p3.axes.get_yaxis().set_visible(False)
p1.set_xlabel('Wavelength [Angstrom]')
p1.set_ylabel('Throughput')

plt.show()
'''