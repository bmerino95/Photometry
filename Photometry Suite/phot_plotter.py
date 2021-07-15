#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 11:54:53 2021

@author: brianmerino

NOTES
a611 does not have f625w
macs0744 does not have f625w and has 775w with ACS and wfc3uvis
macs1423 does not have f625w
"""

"""
2021/04/13 - This code plots the photometry of the galaxy and clumps on the same plot - BMM
"""

import configparser
import glob,os,sys
import time

import numpy as np
import matplotlib.pyplot as plt

#import astropy
from astropy.io import fits
from astropy.wcs import WCS
#from astropy.visualization import simple_norm
from astropy.table import Table
from astropy.table import QTable
#from photutils import CircularAperture
#from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
from photutils import SkyCircularAperture
import pyregion
import configparser

from HST_filters import HST_filt  # This will point to the HST throughput files


t0 = time.time()

# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()
config.read(cfg)

cluster = config.get('CONFIG','cluster')
ID = int(config.get('CONFIG','Id'))


#cluster = 'macs0329'
#ID = '1903'
y_max = 1.5

catalog_file = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/total_phot.cat'%(cluster,ID)
image_path = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/'%(cluster,ID)
plot_path = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/'%(cluster,ID)


catalog = Table.read(catalog_file,format = 'ascii.sextractor')

mag       = {}
magerr    = {}
pivot     = {}
bandwidth = {}

if cluster == 'a611' or cluster =='macs0717' or cluster =='macs0744':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

elif cluster == 'macs1423':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f850lp']

else:
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
                   'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
                   'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']

for Filter in filter_list:
    pf = fits.open(image_path+Filter+'.fits')
    head = pf[0].header
    photplam = head['photplam']
    photbw   = head['photbw']
    pivot[Filter]=float(photplam)
    bandwidth[Filter]=float(photbw)/2


fig = plt.figure(figsize=(8,6))
p1 = fig.add_subplot(111)

for i in range(0,len(catalog)):
    z = catalog['zb'][i]
    
    for item in filter_list:
        mag[item] = catalog[item+'_mag'][i]
        magerr[item] = catalog[item+'_magerr'][i]
        if np.isinf(mag[item]):
            mag[item] = -99
        
    mag0 = [mag[obj] for obj in mag]
    magerr0 = [magerr[obj] for obj in magerr]
    pivots = [pivot[obj] for obj in pivot]
    bandwidths = [bandwidth[obj] for obj in bandwidth]


    # determine range for plotting (filtering bad mags)
    mag1 = np.array(mag0)
    
    if i == 0:
        for j in range(0,len(mag1)):
            if np.isinf(mag1[j]):
                mag1[j]=-99
        
        #mask = mag2 > -99
        mask = mag1 > -50
        mag1_mask = np.compress(mask,mag1,0)
        y0 = np.min(mag1_mask)
        y1 = np.max(mag1_mask) + y_max
    
    obj = catalog['obj'][i]
    p1.errorbar(pivots, mag1, xerr=bandwidths, yerr = magerr0, marker='o', linestyle='', label="Merino_%s"%(obj), zorder=2)


# This section will direct the code to the HST filter throughputs files
x_f105w, y_f105w, x_f110w, y_f110w, x_f125w, y_f125w, x_f140w, y_f140w,\
x_f160w, y_f160w, x_f225w, y_f225w, x_f275w, y_f275w, x_f336w, y_f336w,\
x_f390w, y_f390w, x_f435w, y_f435w, x_f475w, y_f475w, x_f606w, y_f606w,\
x_f625w, y_f625w, x_f775w, y_f775w, x_f814w, y_f814w, x_f850lp, y_f850lp = HST_filt()


p1.set_title('CLASH id: ' + str(cluster)+  ' ' + str(ID) + ' Total,  z = ' + str(z))
p1.legend(loc=2)

p1.set_xlabel('Wavelength [Angstrom]')
p1.set_ylabel('Magnitude [AB]')
p1.set_ylim(y0 - 0.75,y1 + 0.75)
p1.invert_yaxis()

alpha = 0.4
p2 = p1.twinx()
p2.plot(x_f105w, y_f105w, zorder=1, label = 'F105W', c = 'green', alpha = alpha)
p2.plot(x_f110w, y_f110w, zorder=1, label = 'F110W', c = 'green', alpha = alpha)
p2.plot(x_f125w, y_f125w, zorder=1, label = 'F125W', c = 'green', alpha = alpha)
p2.plot(x_f140w, y_f140w, zorder=1, label = 'F140W', c = 'green', alpha = alpha)
p2.plot(x_f160w, y_f160w, zorder=1, label = 'F160W', c = 'green', alpha = alpha)
p2.plot(x_f225w, y_f225w, zorder=1, label = 'F225W', c = 'green', alpha = alpha)
p2.plot(x_f275w, y_f275w, zorder=1, label = 'F275W', c = 'green', alpha = alpha)
p2.plot(x_f336w, y_f336w, zorder=1, label = 'F336W', c = 'green', alpha = alpha)
p2.plot(x_f390w, y_f390w, zorder=1, label = 'F390W', c = 'blue', alpha = alpha)
p2.plot(x_f435w, y_f435w, zorder=1, label = 'F435W', c = 'blue', alpha = alpha)
p2.plot(x_f475w, y_f475w, zorder=1, label = 'F475W', c = 'blue', alpha = alpha)
p2.plot(x_f606w, y_f606w, zorder=1, label = 'F606W', c = 'blue', alpha = alpha)
p2.plot(x_f625w, y_f625w, zorder=1, label = 'F625W', c = 'blue', alpha = alpha)
p2.plot(x_f775w, y_f775w, zorder=1, label = 'F775W', c = 'blue', alpha = alpha)
p2.plot(x_f814w, y_f814w, zorder=1, label = 'F814W', c = 'blue', alpha = alpha)
p2.plot(x_f850lp, y_f850lp, zorder=1, label = 'F850lp', c = 'blue', alpha = alpha)
p2.set_xlim(1000, 17000)
p2.set_ylim(0,0.6)


fontsize = 8
p2.annotate('F105W', xy=(10500,0.02), xytext=(10000,0.06), fontsize=fontsize)
p2.annotate('F110W', xy=(11000,0.02), xytext=(10500,0.02), fontsize=fontsize)
p2.annotate('F125W', xy=(12500,0.02), xytext=(12000,0.04), fontsize=fontsize)
p2.annotate('F140W', xy=(14000,0.02), xytext=(13500,0.06), fontsize=fontsize)
p2.annotate('F160W', xy=(16000,0.02), xytext=(15500,0.02), fontsize=fontsize)
p2.annotate('F225W', xy=(2250,0.02), xytext=(1750,0.02), fontsize=fontsize)
p2.annotate('F275W', xy=(2750,0.02), xytext=(2250,0.04), fontsize=fontsize)
p2.annotate('F336W', xy=(3360,0.02), xytext=(2960,0.06), fontsize=fontsize)
p2.annotate('F390W', xy=(3900,0.02), xytext=(3400,0.02), fontsize=fontsize)
p2.annotate('F435W', xy=(4350,0.02), xytext=(3950,0.04), fontsize=fontsize)
p2.annotate('F475W', xy=(4750,0.02), xytext=(4450,0.06), fontsize=fontsize)
p2.annotate('F606W', xy=(6060,0.02), xytext=(5560,0.02), fontsize=fontsize)
p2.annotate('F625W', xy=(6250,0.02), xytext=(5750,0.04), fontsize=fontsize)
p2.annotate('F775W', xy=(7750,0.02), xytext=(7250,0.06), fontsize=fontsize)
p2.annotate('F814W', xy=(8140,0.02), xytext=(7640,0.02), fontsize=fontsize)
p2.annotate('F850lp', xy=(8500,0.02), xytext=(8000,0.04), fontsize=fontsize)



# This section will display approximately where the emission lines H-alpha
# and [OIII] should be after accounting for the source's redshift.
p3 = p1.twinx()
#z = 0.579

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


#plt.savefig(plot_dir + "/%s_%s.png" % (cluster,obj))
plt.savefig(plot_path + 'total.png')
plt.show()




t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')
