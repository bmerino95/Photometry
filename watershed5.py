#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 11:43:12 2021

@author: brianmerino

This is a modified version of my watershed4.py code that I created during my first visit 
to Carnegie in 2019. This version will eventually be integrated into my working photometry code.
"""

import sys
import pyfits
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage.morphology import watershed
from skimage.feature import peak_local_max
import astropy

from math_utils import RMS,MAD


#np.set_printoptions(threshold=sys.maxsize)

# http://en.wikipedia.org/wiki/Watershed_%28image_processing%29
# http://cmm.ensmp.fr/~beucher/wtshed.html

# watershed.py hlsp_clash_hst_acs_rxj2248_f606w_v1_drz_thumb.fits
# watershed.py hlsp_clash_hst_wfc3ir_rxj2248_f105w_v1_drz_thumb.fits

# watershed3.py hlsp_clash_hst_acs_a383_f625w_v1_drz_thumb.fits hlsp_clash_hst_acs_a383_f625w_v1_wht_thumb.fits


Filter = 'F850lp'

sf = '/Users/brianmerino/Desktop/Apperture_Test/%s.fits'%(Filter)
#nf = sys.argv[2]

ext = 0


pf = astropy.io.fits.open(sf)
image = pf[ext].data
head = pf[ext].header
xl,yl = image.shape

#pf = pyfits.open(nf)
#inv_var = pf[ext].data

#nsigma = 3.0
nsigma = 1.5
#nsigma = 10.0

min_img = np.min(image)
max_img = np.max(image)
print( "min value in image =",min_img)
print( "max value in image =",max_img)


#var = 1./inv_var
#sigma = np.sqrt(var)

#sigma0 = RMS(sigma)

sigma0 = RMS(image)
#sigma0 = 2000
#sigma0 = 0.03       # may want to fix sigma instead of calulate it from the image

print()
print( "sigma for image =",sigma0)




#levels = np.arange(nsigma*sigma,max_img,sigma)
#print levels


threshold = nsigma*sigma0

mask = image > threshold

fig = plt.figure(figsize=(10,10))
p1 = fig.add_subplot(2,2,1)
p1.imshow(image, cmap=cm.jet, vmin=-sigma0,vmax=threshold,interpolation='nearest',origin="lower")
plt.title('Original Image')

# Using http://scipy-lectures.github.io/advanced/image_processing/index.html

# Another guide:
# http://wiki.scipy.org/Cookbook/Watershed
print()
print()
local_maxi = peak_local_max(image,  min_distance=1,
                      threshold_abs=threshold,
                      indices=False, exclude_border=False)

markers_im, nb_markers = ndimage.label(local_maxi)
print( "Number of regions found by peak_local_max =",nb_markers) # how many regions
labels = watershed(-image, markers_im ,mask=mask)
print()
print( "Labels variable =")
print( labels)
print()
print( "Unique values in labels array =")
print( np.unique(labels,return_counts=True))
print()
print( "Max number of labels =",np.max(labels))
print()
print()

###################################################################
# sizes
print()
print()
print( range(nb_markers + 1))
print( "Sizes (pixels) =")
sizes = ndimage.sum(mask, labels, range(nb_markers + 1))
print( sizes) # size of each region (number of pixels)
#print( len(sizes)
print()
print()
print( range(1, nb_markers + 1))
print( "Mean values (DN) =")
mean_vals = ndimage.sum(image, labels, range(1, nb_markers + 1))
print( mean_vals) # mean value of each region
#print( len(mean_vals)
print()
print()
###################################################################

masked = mask*image


p2 = fig.add_subplot(2,2,2)
p2.imshow(markers_im, cmap=cm.jet,interpolation='nearest',origin="lower")
plt.title('Position of Centroids')

p3 = fig.add_subplot(2,2,3)
p3.imshow(labels, cmap=cm.jet,interpolation='nearest',origin="lower")
plt.title('Segmentation Map')

p4 = fig.add_subplot(2,2,4)
p4.imshow(masked, cmap=cm.jet,interpolation='nearest',origin="lower")
plt.title('Segmentation Mask Over Original Image')

hdu = astropy.io.fits.PrimaryHDU(masked,header=head)
hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/Masked_files/Watershed/%s.fits'%(Filter),overwrite=True)

plt.savefig('/Users/brianmerino/Desktop/Apperture_Test/Plots/Watershed/%s.png'%(Filter),overwrite=True)
plt.show()
