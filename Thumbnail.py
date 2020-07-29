#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 22:21:56 2020

@author: brianmerino
"""

import numpy as np
import matplotlib.pyplot as plt
#import sys
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord  
from photutils import SkyCircularAperture
import astropy.units as u
from astropy.table import Table
import time

t0 = time.time()

#Read in CLASH photometry for clumpy galaxies
iden,ra,dec,area = np.loadtxt('/Users/brianmerino/Desktop/CLUSTERS/master_catalog.txt',delimiter='\t',usecols=[3,0,1,4],skiprows=1,unpack=True)
cluster = np.loadtxt('/Users/brianmerino/Desktop/CLUSTERS/master_catalog.txt',dtype='str',delimiter='\t',usecols=[2],skiprows=1)

#Use circular area to calulcate radius
#Create circular aperatures in skycoordinates

coords,radius,sky_apers = [],[],[]
iden_int,rad = [],[]
obj,pixobs = [],[]


for l in range(0,len(area)):
    rad.append(round(np.sqrt(area[l]/np.pi),4))
    iden_int.append(int(iden[l]))
    coords.append(SkyCoord(ra[l], dec[l], unit="deg"))
    radius.append(rad[l] *u.pix)        # pixels
    sky_apers.append(SkyCircularAperture(coords[l], r=radius[l]))

#Makes sure that the correct number of apertures were created
print (len(sky_apers))

filter_list = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f625w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']
filter_list_a611 = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']
macs1423_filter_list = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f775w','f850lp','f105w','f110w','f125w','f140w','f160w']
Clusters = ['a209','a383','a611','macs0329','macs0416','macs0429','macs0717','macs0744','macs1115','macs1149','macs1206','macs1311','macs1423','macs1720','macs1931','macs2129','ms2137','rxj1347','rxj1532','rxj2129']

#Read in CLASH data so that I can compare it against mine
Master_catalog = Table.read('/Users/brianmerino/Desktop/CLUSTERS/master_catalog.txt',format = 'ascii')

filters = {}

#for n in range(0,len(ra)):
for n in range(0,1):
    test_string = cluster[n][1:-1]
    path = "/Users/brianmerino/Desktop/CLUSTERS/"+str(test_string)+"/data/hst/scale_65mas/"
    filter_path = str(path)+"hlsp_clash_hst_"
    
    if test_string == 'a611' or test_string == 'macs0744':
        for i in range(0,len(filter_list_a611)):
            if filter_list_a611[i] == 'f435w' or filter_list_a611[i] == 'f475w' or filter_list_a611[i] == 'f606w'\
            or filter_list_a611[i] == 'f775w' or filter_list_a611[i] == 'f814w' or filter_list_a611[i] == 'f850lp':
                filters['%s'%(filter_list_a611[i])] = str(filter_path)+'acs_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list_a611[i])
                
            elif filter_list_a611[i] == 'f105w' or filter_list_a611[i] == 'f110w' or filter_list_a611[i] == 'f125w'\
            or filter_list_a611[i] == 'f140w' or filter_list_a611[i] == 'f160w':
                filters['%s'%(filter_list_a611[i])] = str(filter_path)+'wfc3ir_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list_a611[i])
                
            elif filter_list_a611[i] == 'f225w' or filter_list_a611[i] == 'f275w' or filter_list_a611[i] == 'f336w'\
            or filter_list_a611[i] == 'f390w':
                filters['%s'%(filter_list_a611[i])] = str(filter_path)+'wfc3uvis_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list_a611[i])

            #obj.append([ra[n],dec[n]])
            obj = [[ra[n],dec[n]]]
            
            pf = astropy.io.fits.open(filters['%s'%(str(filter_list_a611[i]))])
            head = pf[0].header
            data = pf[0].data
            wcs = WCS(filters['%s'%(str(filter_list_a611[i]))])
            norm = simple_norm(data,'sqrt',percent=99.5)
            pw = WCS(header=head)
            pixobjs = pw.wcs_world2pix(obj,1)
            
            x1 = int(round(pixobjs[0][0]- (1.5*rad[n])))
            x2 = int(round(pixobjs[0][0]+ (1.5*rad[n])))
            y1 = int(round(pixobjs[0][1]- (1.5*rad[n])))
            y2 = int(round(pixobjs[0][1]+ (1.5*rad[n])))
            
            subimg = data[y1:y2,x1:x2]
            
            #This part will update the thumbnail's header
            head['crval1'],head['crval2'] = obj[0][0],obj[0][1]
            head['crpix1'] = int(pixobjs[0][0]) - x1
            head['crpix2'] = int(pixobjs[0][1]) - y1
            
            hdu = fits.PrimaryHDU(subimg,header=head)
            hdu.writeto('/Volumes/Home_Drive/Research/Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[n]),filter_list_a611[i]),overwrite=True)
            print('Finishing  %s %s' %(test_string,filter_list_a611[i]))
            
            # The plotting code below is to visually verify that the coordinates are correct.
            fig = plt.figure()
            p1 = fig.add_subplot(111, projection=wcs)
            p1.imshow(data,norm=norm,origin='lower')
            
            xi = pixobjs[0][0] - (1.5*rad[n])
            xf = pixobjs[0][0] + (1.5*rad[n])
            yi = pixobjs[0][1] - (1.5*rad[n])
            yf = pixobjs[0][1] + (1.5*rad[n])
            
            p1.set_xlim(xi,xf)
            p1.set_ylim(yi,yf)
            
            circle = plt.Circle((pixobjs[0][0], pixobjs[0][1]), rad[n], color='r', fill=False, lw=3)
            p1.add_artist(circle)
            
            plt.show()
            
    
    elif test_string == 'macs1423':
        for i in range(0,len(macs1423_filter_list)):
            if macs1423_filter_list[i] == 'f435w' or macs1423_filter_list[i] == 'f475w' or macs1423_filter_list[i] == 'f606w'\
            or macs1423_filter_list[i] == 'f775w' or macs1423_filter_list[i] == 'f850lp':
                filters['%s'%(macs1423_filter_list[i])] = str(filter_path)+'acs_'+str(test_string)+'_%s_v1_drz.fits'%(macs1423_filter_list[i])
                
            elif macs1423_filter_list[i] == 'f105w' or macs1423_filter_list[i] == 'f110w' or macs1423_filter_list[i] == 'f125w'\
            or macs1423_filter_list[i] == 'f140w' or macs1423_filter_list[i] == 'f160w':
                filters['%s'%(macs1423_filter_list[i])] = str(filter_path)+'wfc3ir_'+str(test_string)+'_%s_v1_drz.fits'%(macs1423_filter_list[i])
                
            elif macs1423_filter_list[i] == 'f225w' or macs1423_filter_list[i] == 'f275w' or macs1423_filter_list[i] == 'f336w'\
            or macs1423_filter_list[i] == 'f390w':
                filters['%s'%(macs1423_filter_list[i])] = str(filter_path)+'wfc3uvis_'+str(test_string)+'_%s_v1_drz.fits'%(macs1423_filter_list[i])

            #obj.append([ra[n],dec[n]])
            obj = [[ra[n],dec[n]]]
            
            pf = astropy.io.fits.open(filters['%s'%(str(macs1423_filter_list[i]))])
            head = pf[0].header
            data = pf[0].data
            wcs = WCS(filters['%s'%(str(macs1423_filter_list[i]))])
            norm = simple_norm(data,'sqrt',percent=99.5)
            pw = WCS(header=head)
            pixobjs = pw.wcs_world2pix(obj,1)
            
            x1 = int(round(pixobjs[0][0]- (1.5*rad[n])))
            x2 = int(round(pixobjs[0][0]+ (1.5*rad[n])))
            y1 = int(round(pixobjs[0][1]- (1.5*rad[n])))
            y2 = int(round(pixobjs[0][1]+ (1.5*rad[n])))
            
            subimg = data[y1:y2,x1:x2]
            
            #This part will update the thumbnail's header
            head['crval1'],head['crval2'] = obj[0][0],obj[0][1]
            head['crpix1'] = int(pixobjs[0][0]) - x1
            head['crpix2'] = int(pixobjs[0][1]) - y1
            
            hdu = fits.PrimaryHDU(subimg,header=head)
            hdu.writeto('/Volumes/Home_Drive/Research/Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[n]),macs1423_filter_list[i]),overwrite=True)
            print('Finishing  %s %s' %(test_string,macs1423_filter_list[i]))
            
            # The plotting code below is to visually verify that the coordinates are correct.
            fig = plt.figure()
            p1 = fig.add_subplot(111, projection=wcs)
            p1.imshow(data,norm=norm,origin='lower')
            
            xi = pixobjs[0][0] - (1.5*rad[n])
            xf = pixobjs[0][0] + (1.5*rad[n])
            yi = pixobjs[0][1] - (1.5*rad[n])
            yf = pixobjs[0][1] + (1.5*rad[n])
            
            p1.set_xlim(xi,xf)
            p1.set_ylim(yi,yf)
            
            circle = plt.Circle((pixobjs[0][0], pixobjs[0][1]), rad[n], color='r', fill=False, lw=3)
            p1.add_artist(circle)
            
            plt.show()
    
        
    else:
        for i in range(0,len(filter_list)):
            if filter_list[i] == 'f435w' or filter_list[i] == 'f475w' or filter_list[i] == 'f606w' or filter_list[i] == 'f625w'\
            or filter_list[i] == 'f775w' or filter_list[i] == 'f814w' or filter_list[i] == 'f850lp':
                filters['%s'%(filter_list[i])] = str(filter_path)+'acs_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list[i])
                
                
            elif filter_list[i] == 'f105w' or filter_list[i] == 'f110w' or filter_list[i] == 'f125w' or filter_list[i] == 'f140w'\
            or filter_list[i] == 'f160w':
                filters['%s'%(filter_list[i])] = str(filter_path)+'wfc3ir_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list[i])
                
                
            elif filter_list[i] == 'f225w' or filter_list[i] == 'f275w' or filter_list[i] == 'f336w' or filter_list[i] == 'f390w':
                filters['%s'%(filter_list[i])] = str(filter_path)+'wfc3uvis_'+str(test_string)+'_%s_v1_drz.fits'%(filter_list[i])
            
            
            obj = [[ra[n],dec[n]]]
            
            pf = fits.open(filters['%s'%(str(filter_list[i]))])
            head = pf[0].header
            data = pf[0].data
            wcs = WCS(filters['%s'%(str(filter_list[i]))])
            norm = simple_norm(data,'sqrt',percent=99.5)
            pw = WCS(header=head)
            pixobjs = pw.wcs_world2pix(obj,1)
            
            ###REMOVE THIS LATER
            pixobjs = [[2313.9906,2697.0374]]
            
            #Edit this numbrer to change the size of the thumbnail
            scalar = 7.0
            
            x1 = int(round(pixobjs[0][0]- (scalar*rad[n])))
            x2 = int(round(pixobjs[0][0]+ (scalar*rad[n])))
            y1 = int(round(pixobjs[0][1]- (scalar*rad[n])))
            y2 = int(round(pixobjs[0][1]+ (scalar*rad[n])))
            
            subimg = data[y1:y2,x1:x2]
            
            #This part will update the thumbnail's header
            head['crval1'],head['crval2'] = obj[0][0],obj[0][1]
            head['crpix1'] = int(pixobjs[0][0]) - x1
            head['crpix2'] = int(pixobjs[0][1]) - y1
            
            hdu = fits.PrimaryHDU(subimg,header=head)
            hdu.writeto('/Volumes/Home_Drive/Research/Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[n]),filter_list[i]),overwrite=True)
            print('Finishing %s %s %s' %(test_string,iden[n],filter_list[i]))
            
            
            # The plotting code below is to visually verify that the coordinates are correct.
            fig = plt.figure()
            p1 = fig.add_subplot(111, projection=wcs)
            p1.imshow(data,norm=norm,origin='lower')
            
            xi = pixobjs[0][0] - (scalar*rad[n])
            xf = pixobjs[0][0] + (scalar*rad[n])
            yi = pixobjs[0][1] - (scalar*rad[n])
            yf = pixobjs[0][1] + (scalar*rad[n])
            
            p1.set_xlim(xi,xf)
            p1.set_ylim(yi,yf)
            
            circle = plt.Circle((pixobjs[0][0], pixobjs[0][1]), rad[n], color='r', fill=False, lw=3)
            p1.add_artist(circle)
            
            plt.show()
            #print()
            #print(pixobjs)
            #print()
            #print(p1.get_xlim(),p1.get_ylim())

t1 = time.time()
total = t1-t0
total_min = total/60
print()
print(total, 'seconds')
print(total_min, 'minutes')