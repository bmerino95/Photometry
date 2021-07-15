#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 12:36:48 2021

@author: brianmerino


This is a modified version of my watershed4.py code that I created during my first visit 
to Carnegie in 2019. This version will eventually be integrated into my working photometry code.


This version of my watershed code will read in a region file that contains an estimate of where
the center of the clumps falls and then reports which centroid falls closeset to these regions.


USEAGE:
The code should estimate the size of the clumps and produce a clump mask.
The mask will then be applied to the rest of the thumbnails to conserve flux
"""


import sys
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from astropy.visualization import simple_norm
import astropy
from astropy.wcs import WCS
import configparser
import time

from math_utils import RMS#,MAD
from filter_chooser import chooser  # This will choose filter that best captures halpha

t0 = time.time()

#np.set_printoptions(threshold=sys.maxsize)

# http://en.wikipedia.org/wiki/Watershed_%28image_processing%29
# http://cmm.ensmp.fr/~beucher/wtshed.html

# watershed.py hlsp_clash_hst_acs_rxj2248_f606w_v1_drz_thumb.fits
# watershed.py hlsp_clash_hst_wfc3ir_rxj2248_f105w_v1_drz_thumb.fits

# watershed3.py hlsp_clash_hst_acs_a383_f625w_v1_drz_thumb.fits hlsp_clash_hst_acs_a383_f625w_v1_wht_thumb.fits

def Watershed(filt,x_coord,y_coord,counter):
    
    Filter = filt
    
    #sf = '/Users/brianmerino/Desktop/Apperture_Test/BS_files/bs/%s_%s_%s_bs.fits'%(cluster,ID,Filter)
    sf = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_files/bs/%s_%s_bs.fits'%(cluster,ID,cluster,Filter)
    #nf = sys.argv[2]
    
    ext = 0
    
    
    pf = astropy.io.fits.open(sf)
    image = pf[ext].data
    head = pf[ext].header
    norm = simple_norm(image, 'sqrt', percent=99.5)
    xl,yl = image.shape
    
    # Convert region coords to pixel coords
    w = WCS(head)
    
    x1, y1 = w.wcs_world2pix(x_coord, y_coord,1)
    
    #wcs_world2pix swaps x and y values
    x2 = round(y1.item())
    y2 = round(x1.item())


    # This is used to verify the correct number of coordinate sets have been detected.
    print()
    print('LOOK HERE')
    print(x2,y2)
    print()
    print()


    #pf = pyfits.open(nf)
    #inv_var = pf[ext].data
    
    #nsigma = 3.0
    nsigma = sigma
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
    #p1.imshow(image, cmap=cm.jet, vmin=-sigma0,vmax=threshold,interpolation='nearest',origin="lower")
    p1.imshow(image, cmap='Greys_r', interpolation='nearest',origin="lower")
    plt.title('(a) ' +str(Filter))
    
    
    ###################################################################
    # Using http://scipy-lectures.github.io/advanced/image_processing/index.html
    
    # Another guide:
    # http://wiki.scipy.org/Cookbook/Watershed

    local_maxi = peak_local_max(image,  min_distance=min_distance,
                          threshold_abs=threshold,
                          indices=False, exclude_border=False)
    
    markers_im, nb_markers = ndimage.label(local_maxi)
    
    print()
    print()
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
    # Find which centroids correspond to clumps
    # Find out which centroid is closest to the center of each regions by using the distance formula
    print('Labels that match centroid of visually inspected clumps:')
    
    clump_label = []
    clump_x     = []
    clump_y     = []

    x_top_range    = x2+5
    y_top_range    = y2+5
    x_bottom_range = x2-5
    y_bottom_range = y2-5
        
    dist = 5
    x_coord,y_coord=0,0
    
    for a in range(x_bottom_range,x_top_range):
        for b in range(y_bottom_range,y_top_range):
            if markers_im[a][b]>0:
                
                dist2 = np.sqrt((x2-a)**2+(y2-b)**2)
                
                if dist > dist2:
                    dist=dist2
                    x_coord=a
                    y_coord=b
                
                
                clump_label.append(markers_im[x_coord][y_coord])
                clump_x.append(x_coord)
                clump_y.append(y_coord)
                print('Label:',markers_im[x_coord][y_coord])
                print('X , Y:',x_coord,',',y_coord,'\n')
    
    ###################################################################
    # Drop everything that is not a clump from centroid and segmentation map
    
    for d in range(0,xl):
        for e in range(0,yl):
            if markers_im[d][e] not in clump_label:
                markers_im[d][e] = 0
            
            if labels[d][e] not in clump_label:
                labels[d][e] = 0


    p2 = fig.add_subplot(2,2,2)
    p2.imshow(markers_im, cmap=cm.jet,interpolation='nearest',origin="lower")
    plt.title('(b) Position of Centroids')
    
    p3 = fig.add_subplot(2,2,3)
    p3.imshow(labels, cmap=cm.jet,interpolation='nearest',origin="lower")
    plt.title('(c) Segmentation Map')

    ###################################################################
    # Sizes
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
    # Create new mask to isolate clumps
    
    for f in range(0,xl):
        for g in range(0,yl):
            if labels[f][g] > 0:
                labels[f][g] = 1
    
    masked = labels*image
    
    ###################################################################

    #Old mask that displays the extent of the galaxy
    #masked = mask*image 
    
    
    p4 = fig.add_subplot(2,2,4)
    p4.imshow(masked, cmap='Greys_r',interpolation='nearest',origin="lower")
    plt.title('(d) Segmentation Mask Over Original Image')
    
    #hdu = astropy.io.fits.PrimaryHDU(masked,header=head)
    #hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/Clumps/Watershed/Plots/%s.fits'%(Filter),overwrite=True)
    
    #plt.savefig('/Users/brianmerino/Desktop/Apperture_Test/Clumps/Plots/%s_%s_clump_%s_%s.png'%(cluster,ID,counter,Filter),overwrite=True)
    plt.savefig('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/masked/clump_%s_%s.png'%(cluster,ID,counter,Filter),overwrite=True)
    plt.show()
    
    
    return labels

###################################################################
###################################################################
#User inputs

# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()    
config.read(cfg)

cluster = config.get('CONFIG','cluster')
ID = int(config.get('CONFIG','Id'))
z = float(config.get('CONFIG','z'))

#cluster = 'macs0329'
#ID = 1903
min_distance = 1
sigma = 1.5
wht = 1
more_than_one = 1 #Set this to 0 if a region file has 1 target.
###################################################################
###################################################################
if cluster == 'a611' or cluster =='macs0717' or cluster =='macs0744':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

elif cluster == 'macs1423':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w',  'f850lp']

else:
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
                   'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
                   'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']


###################################################################
# Read in region file
import pyregion
region_path = '/Users/brianmerino/Desktop/Individual_Clumps/%s/'%(cluster)
region_name = region_path+'%s.reg'%(ID)

#region_name = '/Users/brianmerino/Desktop/Apperture_Test/a383-1686/bulge.reg'

r = pyregion.open(region_name)
print(r)


###################################################################
# Find out which filter to mask a mask with
chosen = chooser(z)
print(chosen)
chosen = 'f606w'

###################################################################
# Save coordinates of clumps
x ,y = [],[]
for i in range(0,len(r)):
    x.append(r[i].coord_list[0])
    y.append(r[i].coord_list[1])

print(x)
print(y)



###################################################################
# Run watershed over filter where halpha falls

# This counter will assign a number to each clump
counter = 0

if more_than_one:
    end = len(x)-1
else:
    end = len(x)-0

# The first for loop isolate each clump individually and return a mask
for i in range(0,end):
    mask = Watershed(chosen,x[i],y[i], counter)
    
    
    # The second for loop will apply the mask to each thumbnail and save the isolated clump
    for item in filter_list:
        sf = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/BS_files/bs/%s_%s_bs.fits'%(cluster,ID,cluster,item)
        #sf = '/Users/brianmerino/Desktop/Apperture_Test/a383-1686/BS_files/bs/%s_bulge_%s_bs.fits'%(cluster,item)
        ext = 0
        pf = astropy.io.fits.open(sf)
        image = pf[ext].data
        head  = pf[ext].header
        
        clump = mask*image
        
        hdu = astropy.io.fits.PrimaryHDU(clump,header=head)
        hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/masked/clump_%s_%s.fits'%(cluster,ID,counter,item),overwrite=True)
        #hdu.writeto('/Users/brianmerino/Desktop/Apperture_Test/a383-1686/masked/bulge_%s.fits'%(item),overwrite=True)
        
        
        
        if wht:
            sf_wht = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/Thumbnails/%s_wht.fits'%(cluster,ID,item)
            pf_wht=astropy.io.fits.open(sf_wht)
            image_wht = pf_wht[ext].data
            head_wht = pf_wht[ext].header
            
            clump_wht = mask*image_wht
            
            hdu_wht = astropy.io.fits.PrimaryHDU(clump_wht,header=head_wht)
            hdu_wht.writeto('/Users/brianmerino/Desktop/Apperture_Test/%s-%s/masked/clump_%s_%s_wht.fits'%(cluster,ID,counter,item),overwrite=True)
        
        
    
    counter += 1

print("\n")
print(str(counter)+" clump(s)")




t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')