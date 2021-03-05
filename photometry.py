#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:29:17 2020

@author: brianmerino

This code is meant to automate my previous photometry code. It should plot the photometry for each filter and then
read in the master catalog that I created for my clumpy galaxies and compare the two results. This version of the 
code does NOT include background subtraction. Instead, I created a separate code to subtract the background since 
that process was very time intensive. 


NEEDS:
Master catalog (for CLASH information)
Interesting Sources region file

Eventually:
Molino and Conner information will be added to final plots

NOTES
a611 does not have f625w
macs0744 does not have f625w and has 775w with ACS and wfc3uvis
macs1423 does not have f625w
"""

"""
2020/08/22 - Cleaned up code - GLW
"""

import configparser   # Python 2.7?
import glob,os,sys
import time

import numpy as np
import matplotlib.pyplot as plt

#import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.table import Table
#from photutils import CircularAperture
#from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
from photutils import SkyCircularAperture


t0 = time.time()

####################################
# constants
####################################
c = 2.9979E10       # cm/s
Ang = 1E-8          # cm
Jy = 1.0E-23        # erg/s/cm^2/Hz
mJy = 1e-3          # Jy
uJy = 1e-6          # Jy
Mpc = 3.086e24      # cm

## PHOTFLAM: inverse sensitivity (erg * cm**(-2) * s**(-1) Angstrom**(-1)
## PHOTPLAM: Pivot wavelength

####################################
# variables and flags
####################################
cfg = "photometry.cfg"
bkgd_subtract = 1   # use background subtracted frames for photometry
ext = 0
####################################


# read-in cfg file with paramters 
config = configparser.ConfigParser()
config.read(cfg)
#path to background subtracted iamges
img_dir = config.get('CONFIG','img_dir')
sub_img_dir = config.get('CONFIG','sub_img_dir')
plot_dir = config.get('CONFIG','plot_dir')
catalog_file = config.get('CONFIG','catalog_file')


def hst_phot(photflam,photplam,dn,dn_err=None,unit="Jy"):
    ## http://www.stsci.edu/hst/acs/analysis/zeropoints
    ABMAG_ZP=-2.5*np.log10(photflam)-5*np.log10(photplam)-2.408 
    #print("%.1f" % (photplam),)
    #print("ABMAG_ZP =", ABMAG_ZP)
    m = -2.5*np.log10(dn) + ABMAG_ZP
    #print(m)
    #C = photflam                              # erg/s/cm^2/Ang
    #C = photflam/Ang*(photplam*Ang)**2/c      # erg/s/cm^2/Hz
    C = photflam/Ang*(photplam*Ang)**2/c/Jy   # Jy
    if unit == "Jy": Fnu = dn*C
    elif unit == "mJy": Fnu = dn*C/mJy
    elif unit == "uJy": Fnu = dn*C/uJy
    ABmag = -2.5*np.log10(Fnu*Jy) - 48.60
    #print("%.3f" % (ABmag),)
    if dn_err:
        if unit == "Jy": Fnu_err = dn_err*C
        elif unit == "mJy": Fnu_err = dn_err*C/mJy
        elif unit == "uJy": Fnu_err = dn_err*C/uJy

        ABmag_err = 2.5*(Fnu_err/(Fnu*np.log(10)))
        #print("%.3f" % (ABmag_err))


        #print("%.4f %.6e %.6e" % (photplam/1E4,Fnu,Fnu_err))
        #print("%.6e %.6e" % (Fnu,Fnu_err))
        #return Fnu,Fnu_err
        return ABmag,ABmag_err
    else:
        #print("%.6e" % Fnu)
        #return Fnu
        return ABmag


# Read in CLASH data so that I can compare it against mine
mcatalog = Table.read(catalog_file,format = 'ascii.sextractor')

# unique set of clusters

clusters = list(set(mcatalog['Cluster']))
clusters.sort()
print()
print(clusters)
print()
print(len(clusters))
print()

#print(mcatalog)
#print(mcatalog['RA'])
#print(mcatalog['area'])


#print(mcatalog['f606w_mag'])
##

# Master dictionary, everything will live here
mdict = {}


'''
NOTE: For some reason, one of the filters (f814w) in macs1423 is giving me an 
issue that none of the other images are giving me. I might need to exclude it from this code
if I cant figure out what is going on. 
'''

print(mcatalog.colnames)


filts = []
mag_cols = []
magerr_cols = []

# generate lists of filts and magnitude and magnitude error columns from master catalog
for col in mcatalog.colnames:
    if col.endswith('mag'):
        mag_cols.append(col)
        magerr_cols.append(col+"err")
        filts.append(col.replace('_mag',''))

#print(filts)
#print(len(filts))
#print(mag_cols)
#print(magerr_cols)


# loop over cluster (with the goal of opening each image only once)
#for i,cluster in enumerate(clusters):

#    for j,filt in enumerate(filts):



print(mcatalog)


# test with the first 4 catalog entries
#mcatalog = mcatalog[:4]
#print(mcatalog)


# loop over catalog object
for i,mcat_row in enumerate(mcatalog):
    #print(i,mcat_row)
    #print(mcat_row[mag_cols])

    cluster = mcat_row['Cluster']
    obj = mcat_row['Molino_ID']
    #print(cluster)

    mdict[obj] = {}
    mdict[obj]["keys"] = {"mag":[], "magerr":[], "pivot":[], "bandwidth":[]}

    for j,filt in enumerate(filts):
        print(j,filt)

        # without f625w
        if cluster == 'a611' or cluster == 'macs0744':
            if filt == "f625w": continue

        # without f625w and f814w
        elif cluster == 'macs1423':
            if filt == "f625w" or filt == 'f814w': continue
        

        if bkgd_subtract:
            # background subtracted frames                                      
            drz_file = '/Users/brianmerino/Desktop/grid_test/Elliptical_Backgrounds/' + mcat_row['Molino_ID'] +'_'+ filt + '.fits'
            wht_file = '/Users/brianmerino/Desktop/Galaxies_w_appertures/' + cluster +'_' + mcat_row['Molino_ID'] +'_'+ filt + '_wht.fits'
            print("nothing here now...")

        else:
            # original frames
            #drz_str = img_dir + "/" + cluster + '/data/hst/scale_65mas/*_%s_v1_drz.fits' %(filt)
            #wht_str = img_dir + "/" + cluster + '/data/hst/scale_65mas/*_%s_v1_wht.fits' %(filt)
            
            drz_str = '/Users/brianmerino/Desktop/CLUSTERS/' + cluster + '/data/hst/scale_65mas/*_%s_v1_drz.fits' %(filt)
            wht_str = '/Users/brianmerino/Desktop/CLUSTERS/' + cluster + '/data/hst/scale_65mas/*_%s_v1_wht.fits' %(filt)
            
            
            drz_file = glob.glob(drz_str)[0]
            wht_file = glob.glob(wht_str)[0]
            
            print(drz_file)
            print(wht_file)
        
        #print(os.path.exists(drz_file))
        #print(os.path.exists(wht_file))
        
        pf = fits.open(drz_file)
        image = pf[ext].data
        head = pf[ext].header
        photflam = head['photflam']
        photplam = head['photplam']
        photbw   = head['photbw']
        
        pf = fits.open(wht_file)
        #  inverse variance weight images
        wht = pf[ext].data
        
        noise = 1./np.sqrt(wht)
        
        wcs = WCS(head)
        wcs.sip = None
        
        ra = mcat_row['RA']
        dec = mcat_row['DEC']
        area = mcat_row['area']
        
        coords = SkyCoord(ra, dec, unit="deg")
        radius = round(np.sqrt(area/np.pi),4) * u.pix # pixels
        sky_aper = SkyCircularAperture(coords, r=radius)
        ap_sum = aperture_photometry(image,sky_aper,error=noise,wcs=wcs) # flux + flux error in counts or electrons
        flux = ap_sum["aperture_sum"]
        flux_err = ap_sum["aperture_sum_err"]
        mag, magerr = hst_phot(photflam, photplam, flux, dn_err=flux_err) # AB magnitude
        
        
        #print(mag,magerr)
        #print(photplam)
        
        # store values in a dict --> eventually write to a table...
        mdict[obj][filt+'_mag'] = float(mag)
        mdict[obj][filt+'_magerr'] = float(magerr)
        mdict[obj][filt+'_pivot'] = float(photplam)
        mdict[obj][filt+'_bandwidth'] = float(photbw)/2.
        
        mdict[obj]["keys"]["mag"].append(filt+'_mag')
        mdict[obj]["keys"]["magerr"].append(filt+'_magerr')
        mdict[obj]["keys"]["pivot"].append(filt+'_pivot')
        mdict[obj]["keys"]["bandwidth"].append(filt+'_bandwidth')
        
        # subtract Molino phot from Merino's phot
        
        
        
        
    fig = plt.figure()
    p1 = fig.add_subplot(111)

    # dictionary and table keys
    mag_keys = mdict[obj]["keys"]["mag"]
    magerr_keys = mdict[obj]["keys"]["magerr"]
    pivot_keys = mdict[obj]["keys"]["pivot"]
    bandwidth_keys = mdict[obj]["keys"]["bandwidth"]

    pivot = [mdict[obj][key] for key in pivot_keys]
    bandwidth = [mdict[obj][key] for key in bandwidth_keys]

    # Merino photometry (from dictionary)
    mag1 = [mdict[obj][key] for key in mag_keys]
    magerr1 = [mdict[obj][key] for key in magerr_keys]

    # Molino photometry (from table)
    mag2 = [mcat_row[key] for key in mag_keys]
    magerr2 = [mcat_row[key] for key in magerr_keys]

    # Difference between photometry -> Merino - Molino
    delta = [mdict[obj][key] - mcat_row[key] for key in mag_keys]

    # determine range for plotting (filtering bad mags)
    mag2 = np.array(mag2)
    mask = mag2 != 99
    mag2_mask = np.compress(mask,mag2,0)
    y0 = np.min(mag2_mask)
    y1 = np.max(mag2_mask)
    #print mag2
    #print mag2_mask
    #print y0
    #print y1


    p1.errorbar(pivot, mag1, yerr=magerr1, xerr=bandwidth, marker='o', linestyle='', label="B. Merino")
    p1.errorbar(pivot, mag2, yerr=magerr2, xerr=bandwidth, marker='o', linestyle='', label="Molino")
    p1.legend(loc=4)

    p1.set_xlabel('Wavelength [Angstrom]')
    p1.set_ylabel('Magnitude [AB]')
    p1.set_title('Molino id: ' + str(obj))
    p1.set_ylim(y0 - 0.75,y1 + 0.75)
    p1.invert_yaxis()

    plt.savefig(plot_dir + "/%s_%s.png" % (cluster,obj) )
    #plt.show()

    print(cluster, mcat_row['Molino_ID'], ra, dec, area, mcat_row['redshift'], delta, \
          file = open('/Users/brianmerino/Desktop/grid_test/elliptical_phot_2.cat','a'))

        
t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')