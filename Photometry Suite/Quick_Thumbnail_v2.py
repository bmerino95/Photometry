#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 16:09:55 2021

@author: brianmerino

This code will create a thumbnail of a source. This code will require the source's coordinates
and the size of the thumbnail.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord  
import configparser
import time

t0 = time.time()

def Quick(cluster,Filter,ra,dec,rad):
    path = '/Users/brianmerino/Desktop/CLUSTERS/%s/data/hst/scale_65mas/'%(cluster)
    
    if Filter == 'f105w' or Filter == 'f110w' or Filter == 'f125w' or Filter =='f140w' or Filter =='f160w':
        file = 'hlsp_clash_hst_wfc3ir_%s_%s_v1_drz.fits'%(cluster,Filter)
        file_wht = 'hlsp_clash_hst_wfc3ir_%s_%s_v1_wht.fits'%(cluster,Filter)
    
    elif Filter == 'f225w' or Filter == 'f275w' or Filter =='f336w' or Filter =='f390w':
        file = 'hlsp_clash_hst_wfc3uvis_%s_%s_v1_drz.fits'%(cluster,Filter)
        file_wht = 'hlsp_clash_hst_wfc3uvis_%s_%s_v1_wht.fits'%(cluster,Filter)
    
    else:
        file = 'hlsp_clash_hst_acs_%s_%s_v1_drz.fits'%(cluster,Filter)
        file_wht = 'hlsp_clash_hst_acs_%s_%s_v1_wht.fits'%(cluster,Filter)
    
    #name = Filter
    
    #Source's coordinates
    obj = [[ra, dec]] #Ra and Dec
    
    #Length is the amount of pixels to the center of the x or y axis
    #length = rad*4
    scale  = 1
    
    #Optional radius for plotting component at the end in pixels
    rad = rad
    
    #Load in image
    pf       = fits.open(path+file)
    pf_wht   = fits.open(path+file_wht)
    head     = pf[0].header
    head_wht = pf_wht[0].header
    data     = pf[0].data
    data_wht = pf_wht[0].data
    wcs      = WCS(path+file)
    norm     = simple_norm(data, 'sqrt', percent=99.5)
    pw       = WCS(header=head)
    pixobjs  = pw.wcs_world2pix(obj,1)
    
    
    x1 = int(round(pixobjs[0][0] - (scale*length)))
    x2 = int(round(pixobjs[0][0] + (scale*length)))
    y1 = int(round(pixobjs[0][1] - (scale*length)))
    y2 = int(round(pixobjs[0][1] + (scale*length)))
    
    
    subimg = data[y1:y2,x1:x2]
    subimg_wht = data_wht[y1:y2,x1:x2]
    
    
    #This part will update the thumbnail's header
    head['crval1'],head['crval2'] = obj[0][0],obj[0][1]
    head['crpix1'] = int(pixobjs[0][0]) - x1
    head['crpix2'] = int(pixobjs[0][1]) - y1
    
    head_wht['crval1'],head_wht['crval2'] = obj[0][0],obj[0][1]
    head_wht['crpix1'] = int(pixobjs[0][0]) - x1
    head_wht['crpix2'] = int(pixobjs[0][1]) - y1
    
    
    hdu = fits.PrimaryHDU(subimg,header=head)
    #save_to = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s.fits'%(cluster,Id,Filter)
    hdu.writeto(save_to,overwrite=True)
    
    hdu_wht = fits.PrimaryHDU(subimg_wht,header=head_wht)
    #save_to_wht = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s_wht.fits'%(cluster,Id,Filter)
    hdu_wht.writeto(save_to_wht,overwrite=True)
    
    print()
    print('Thumbnail saved as: '+save_to)
    
    
    # The plotting code below is to visually verify that the coordinates are correct.
    fig = plt.figure()
    p1 = fig.add_subplot(111, projection=wcs)
    p1.imshow(data,norm=norm,origin='lower',cmap='Greys_r')
    
    xi = int(round(pixobjs[0][0]- (scale*length)))
    xf = int(round(pixobjs[0][0]+ (scale*length)))
    yi = int(round(pixobjs[0][1]- (scale*length)))
    yf = int(round(pixobjs[0][1]+ (scale*length)))
    
    p1.set_title('Source')
    p1.set_xlim(xi,xf)
    p1.set_ylim(yi,yf)
    
    
    circle = plt.Circle((pixobjs[0][0], pixobjs[0][1]), rad, color='r', fill=False, lw=2)
    p1.add_artist(circle)
    
    #save_to2 = save_to[:-4]+'png'
    save_to2 = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/%s.png'%(cluster,Id,item)
    plt.savefig(save_to2,overwrite=True)
    plt.show()
    
# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()
config.read(cfg)

cluster = config.get('CONFIG','cluster')
Id = int(config.get('CONFIG','Id'))
ra  = float(config.get('CONFIG','ra'))
dec = float(config.get('CONFIG','dec'))
rad = float(config.get('CONFIG','rad'))


'''cluster = 'macs0329'
Id = '1903'
ra  = 52.41497083
dec = -2.214673056
rad = 35'''
length = rad*1

if cluster == 'a611' or cluster =='macs0717' or cluster == 'macs0744':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

elif cluster == 'macs1423':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f850lp']

else:
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
                   'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
                   'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']


for item in filter_list:
    save_to = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s.fits'%(cluster,Id,item)
    save_to_wht = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s_wht.fits'%(cluster,Id,item)
    Quick(cluster,item,ra,dec,rad)



t1 = time.time()
total = t1-t0
total_min = total/60
print()
print(total, 'seconds')
print(total_min, 'minutes')




