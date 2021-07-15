#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:43:36 2021

@author: brianmerino

This code will test how well photutils works at creating segmentation maps.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy
from photutils import Background2D
from astropy.visualization import simple_norm
from photutils import make_source_mask
from astropy.stats import mad_std
import configparser
import time

t0 = time.time()

# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()
config.read(cfg)

cluster = config.get('CONFIG','cluster')
ID = int(config.get('CONFIG','Id'))


#cluster = 'macs0329'
#ID = '1903'
Filter = 'f105w'
seg_num =2
sigma_0 = 1.5
npixels = 1

#Read in the fits file
image = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_Files/bs/%s_%s_bs.fits'%(cluster,ID,cluster,Filter)

# Open the fits files
pf = astropy.io.fits.open(image)

# Read the headers
head = pf[0].header

# Access the data
data = pf[0].data

# This will adjust the strech for the plot
norm = simple_norm(data, 'linear', percent=99.5)

#plt.imshow(data, origin='lower', cmap='Greys_r', norm=norm)



from photutils import detect_threshold
threshold = detect_threshold(data, nsigma=sigma_0)


from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
sigma = sigma_0 * gaussian_fwhm_to_sigma  # FWHM = 3.
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
segm = detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)
#segm = detect_sources(data, threshold, npixels=5)

#plt.imshow(segm, origin='lower', cmap='Greys_r')

from photutils import deblend_sources
segm_deblend = deblend_sources(data, segm, npixels=npixels,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.001)

print(segm_deblend)
#np.where(segm_deblend == 1, segm_deblend.labels, 0*segm_deblend.labels)

#np.equal(segm_debled=1, segm_deblend)


test = np.where(np.not_equal(segm, seg_num), 0, 1)

'''
#When the galaxy is broken up into parts, I can conbine them like this
test1 = np.where(np.not_equal(segm, seg_num), 0, 1)
test2 = np.where(np.not_equal(segm, 3), 0, 1)
#test3 = np.where(np.not_equal(segm, 5), 0, 1)
#test4 = np.where(np.not_equal(segm, 8), 0, 1)
#test5 = np.where(np.not_equal(segm, 6), 0, 1)

test = test1+test2#+test3+test4+test5
'''

masked = test*data

hdu = astropy.io.fits.PrimaryHDU(masked,header=head)
hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/segmentation_%s.fits'%(cluster,ID,Filter),overwrite=True)


#test = segm_deblend

#test[ind] = 0

from photutils import source_properties
#cat = source_properties(data, segm_deblend)
cat = source_properties(data, segm)

tbl = cat.to_table()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 8), sharey=True)
ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
ax1.set_title('(a) %s'%(Filter))
#cmap = segm.make_cmap(seed=123)
ax2.imshow(segm, origin='lower', interpolation='nearest')

copy = segm

for i in range(0,len(tbl['cxx'])):
    ax2.scatter(tbl['xcentroid'],tbl['ycentroid'], marker='o', c='r')

ax2.set_title('(b) Segmentation Image')

#ax3.imshow(segm_deblend, origin='lower', interpolation='nearest')
#ax3.set_title('Deblended Segmentation Image')

ax3.imshow(masked, origin='lower', interpolation='nearest', cmap='Greys_r')
ax3.set_title('(c) Masked Image')

plt.savefig('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/segmentation_%s.png'%(cluster,ID,Filter),overwrite=True)

plt.show()


print()
print(tbl)
print()

for i in range(0,len(tbl['cxx'])):
    print(i+1,(tbl['xcentroid'][i], tbl['ycentroid'][i]))

########
#This part of the code should apply the mask to the other thumbnails.
########
if cluster == 'a611' or cluster =='macs0717' or cluster =='macs0744':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
               'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
               'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

elif cluster == 'macs1423':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f850lp']

else:
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
               'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
               'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']



#This loop will apply the mask to each filter.
for item in filter_list:
    #Read in the fits file
    image2 = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s.fits'%(cluster,ID,item)
    
    # Open the fits files
    pf2 = astropy.io.fits.open(image2)
    
    # Read the headers
    head2 = pf2[0].header
    
    # Access the data
    data2 = pf2[0].data
    
    wht_image = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s_wht.fits'%(cluster,ID,Filter)
    pf_wht = astropy.io.fits.open(wht_image)
    head_wht = pf_wht[0].header
    data_wht = pf_wht[0].data
    
    final_mask = test*data2
    final_mask_wht = test*data_wht
    
    hdu = astropy.io.fits.PrimaryHDU(final_mask,header=head2)
    hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/masked/galaxy_%s.fits'%(cluster,ID,item),overwrite=True)
    
    hdu_wht = astropy.io.fits.PrimaryHDU(final_mask_wht,header=head_wht)
    hdu_wht.writeto('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/masked/galaxy_%s_wht.fits'%(cluster,ID,item),overwrite=True)
    
    
    plt.imshow(final_mask,origin='lower',cmap='Greys_r')
    plt.title(item)
    plt.show()



plt.imshow(copy, origin='lower', interpolation='nearest')
plt.show()

print('X: '+str(head['CRPIX1'])+'\t Y: '+str(head['CRPIX2']))


t1 = time.time()
total = t1-t0
total_min = total/60
print()
print(round(total,3), 'seconds')
print(round(total_min,3), 'minutes')



