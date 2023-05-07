#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 14:08:44 2021

@author: brianmerino

This code is meant to be a quicker version of my Background_w_mask.py code. Instead of running over 
an entire catalog of galaxies, it runs over a few at a time.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy
from photutils import Background2D
from astropy.visualization import simple_norm
from photutils import make_source_mask
from astropy.stats import mad_std
from astropy.table import Table
import configparser
import time

t0 = time.time()


def Background(image, head, name, norm, sigma):
    
    # I am defining variables x and y to be the box size that Background2D uses.
    # Greg recommended using around 36 boxes in total, so I am taking 2*head['CRPIX1'] to be
    # the width of the image, then I divide by 3 so that 3x3 boxes will cover the image.
    x = int(2*head['CRPIX1']/6)
    y = int(2*head['CRPIX2']/6)
    
    sigma = sigma
    
    mask = make_source_mask(image, nsigma=sigma, npixels=3, dilate_size=10)
    bkg = Background2D(image, (x,y), sigma_clip = None, mask = mask,exclude_percentile=50)
    #bkg = Background2D(image, (x,y), mask=mask)
    
    
    # This adjusts the stretch for the background image
    bkg_min, bkg_max = np.min(bkg.background), np.max(bkg.background)
    
    # This changes the scale of the color bars. The default is too big
    scale = 0.85
    
    
    bs  = image - bkg.background
    
    fig, axs = plt.subplots(2, 4, figsize=(16,7))
    
    # This plots the original image with colorbar
    im1 = axs[0,0].imshow(image, origin='lower', cmap='Greys_r', norm=norm)
    axs[0,0].set_title('(a) %s'%(name))
    fig.colorbar(im1, ax=axs[0,0], shrink=scale)
    
    
    # This plots the background
    im3 = axs[0,1].imshow(bkg.background, origin='lower', cmap='Greys_r', norm=norm)
    axs[0,1].set_title('(b) Background')
    fig.colorbar(im3, ax=axs[0,1], shrink=scale)
    
    
    # This plots the background subtracted image and saves it as a new fits file
    im4 = axs[0,2].imshow(bs, origin='lower', cmap='Greys_r', norm=norm)
    axs[0,2].set_title('(c) Background Subtracted')
    fig.colorbar(im4, ax=axs[0,2], shrink=scale)
    
    hdu2 = astropy.io.fits.PrimaryHDU(bs,header=head)
    hdu2.writeto(save_bs,overwrite=True)
    
    
    # This plots the background
    im2 = axs[0,3].imshow(bkg.background, origin='lower', cmap='Greys_r' , vmin=bkg_min, vmax=bkg_max)
    axs[0,3].set_title('(d) Background')
    fig.colorbar(im2, ax=axs[0,3], shrink=scale)
    
    
    # This plots the mask
    axs[1,0].imshow(mask, origin='lower', cmap='Greys_r')
    axs[1,0].set_title('(e) %s sigma mask'%(sigma))
    #fig.colorbar(im2, ax=axs[0,2], shrink=scale)
    ratio = 0.8
    axs[1,0].set_aspect(1.0/axs[1,0].get_data_ratio()*ratio)
    
    
    # This plots a histogram of the pixels in the original image
    axs[1,1].hist(data.flatten(),35, label = 'Original')
    axs[1,1].set_title('(f) Histogram')
    axs[1,1].legend()
    #This histogram kept streching vertically. These two lines resolve that issue.
    ratio = 0.85
    axs[1,1].set_aspect(1.0/axs[1,1].get_data_ratio()*ratio)
    
    
    axs[1,2].hist(bs.flatten(),35, label = 'B.S.', color='orange')
    axs[1,2].set_title('(g) Histogram')
    axs[1,2].legend()
    #This histogram kept streching vertically. These two lines resolve that issue.
    ratio = 0.85
    axs[1,2].set_aspect(1.0/axs[1,2].get_data_ratio()*ratio)
    
    
    # This plots the original image with the grid used for background subtraction overlayed on top.
    # This plot needs to be last because the gridsize keeps applying itself to the last plot in the subplot
    im6 = axs[1,3].imshow(image, origin='lower', cmap='Greys_r', norm=norm)
    bkg.plot_meshes(outlines=True, color='red')
    axs[1,3].set_title('(h) Mesh Grid')
    fig.colorbar(im6, ax=axs[1,3], shrink=scale)
    
    plt.savefig(save_png)
    plt.show()
    
    print()
    print('Background = ',mad_std(image))
    
    hdu1 = astropy.io.fits.PrimaryHDU(bkg.background,header=head)
    hdu1.writeto(save_bkg,overwrite=True)
    
    


# List of all the filters that I am using
Filter_1 = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
          'f390w', 'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']

#Filter_1 = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
#            'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']

# Subset of filters for a611 and macs0744
Filter_2 = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
          'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

# Subset of filters for macs1423
Filter_3 = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
          'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f850lp']

# Read in catalog of clumpy galaxies
#catalog = Table.read('/Users/brianmerino/Desktop/Regions/Clumpy_ids/master_catalog_clash.txt',format='ascii.sextractor')

sigma = 1.5

# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()
config.read(cfg)

cluster = config.get('CONFIG','cluster')
ID = int(config.get('CONFIG','Id'))


#cluster = 'macs0329'
#ID = '1903'

if cluster == 'a611' or cluster == 'macs0717' or cluster == 'macs0744':
    FILTER = Filter_2
    
elif cluster =='macs1423':
    FILTER = Filter_3

else:
    FILTER = Filter_1

for item in FILTER:
    save_bs = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_files/bs/%s_%s_bs.fits'%(cluster,ID,cluster,item)
    save_bkg = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_files/bkg/%s_%s_bkg.fits'%(cluster,ID,cluster,item)
    save_png = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_files/%s_%s_bs.png'%(cluster,ID,cluster,item)
    
    #Read in the fits files
    image = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s.fits'%(cluster,ID,item)
    
    # Open the fits files
    pf = astropy.io.fits.open(image)
    
    # Read the headers
    head = pf[0].header
    
    # Access the data
    data = pf[0].data
    
    # This will adjust the strech for the plot
    norm = simple_norm(data, 'linear', percent=99.5)
    
    
    # Call the Background function
    # Background(image, head, name, norm, sigma)
    
    Background(data, head, '%s'%(item) , norm, sigma)


    print()
    print('Background subtracted image saved to '+save_bs)
    print('Background saved to ' +save_bkg)
    print()





print('DONE!')

t1 = time.time()
total = t1-t0
total_min = total/60
print()
print(round(total,3), 'seconds')
print(round(total_min,3), 'minutes')



