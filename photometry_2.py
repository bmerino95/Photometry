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

from HST_filters import HST_filt  # This will point to the HST throughput files


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
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/photometry.cfg"
bkgd_subtract = 1   # use background subtracted frames for photometry
ext = 0
####################################


# read-in cfg file with paramters 
config = configparser.ConfigParser()
config.read(cfg)
#path to background subtracted iamges
img_dir = config.get('CONFIG','img_dir')
sub_img_dir = config.get('CONFIG','sub_img_dir')
seg_img_dir = config.get('CONFIG','seg_img_dir')
plot_dir = config.get('CONFIG','plot_dir')
catalog_file = config.get('CONFIG','catalog_file')
hst_filt = config.get('CONFIG','hst_filt')


def hst_phot(photflam,photplam,dn,dn_err=None,unit="Jy"):
    ## http://www.stsci.edu/hst/acs/analysis/zeropoints
    #ABMAG_ZP=-2.5*np.log10(photflam)-5*np.log10(photplam)-2.408 
    #print("%.1f" % (photplam),)
    #print("ABMAG_ZP =", ABMAG_ZP)
    #m = -2.5*np.log10(dn) + ABMAG_ZP
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
print(clusters)
print(len(clusters))

#print(mcatalog)
#print(mcatalog['RA'])
#print(mcatalog['area'])


#print(mcatalog['f606w_mag'])
##

# Master dictionary, everything will live here
mdict = {}
mdict_seg = {}


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
    obj = mcat_row['id']
    z = mcat_row['redshift']
    #print(cluster)

    mdict[obj] = {}
    mdict[obj]["keys"] = {"mag":[], "magerr":[], "pivot":[], "bandwidth":[]}
    
    mdict_seg[obj] = {}
    mdict_seg[obj]["keys"] = {"mag":[], "pivot":[], "bandwidth":[]}

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
            #drz_file = sub_img_dir + cluster + '_%s.fits' % filt
            #wht_file = sub_img_dir + cluster + '_%s.fits' % filt   ????
            
            # I am running a test that focuses on a383, so I slightly modified the two lines above
            drz_file = sub_img_dir + 'drz/%s.fits' % filt
            wht_file = sub_img_dir + 'wht/%s.fits' % filt
            
            
            drz_seg_str = seg_img_dir + 'drz/%s.fits' % filt
            wht_seg_str = seg_img_dir + 'drz/%s.fits' % filt


        else:
            # original frames
            drz_str = img_dir + "/" + cluster + '/data/hst/scale_65mas/*_%s_v1_drz.fits' % filt
            wht_str = img_dir + "/" + cluster + '/data/hst/scale_65mas/*_%s_v1_wht.fits' % filt
            
            drz_file = glob.glob(drz_str)[0]
            wht_file = glob.glob(wht_str)[0]
            
            print(drz_file)
            print(wht_file)

        #print(os.path.exists(drz_file))
        #print(os.path.exists(wht_file))



        #This part of the code will perform circular apperture photometry on the source
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


        # store values in a dict --> eventually write to a table...
        mdict[obj][filt+'_mag'] = float(mag)
        mdict[obj][filt+'_magerr'] = float(magerr)
        mdict[obj][filt+'_pivot'] = float(photplam)
        mdict[obj][filt+'_bandwidth'] = float(photbw)/2.

        mdict[obj]["keys"]["mag"].append(filt+'_mag')
        mdict[obj]["keys"]["magerr"].append(filt+'_magerr')
        mdict[obj]["keys"]["pivot"].append(filt+'_pivot')
        mdict[obj]["keys"]["bandwidth"].append(filt+'_bandwidth')



        #This part of the code will use segementation maps to measure the photometry
        pf = fits.open(drz_seg_str)
        image = pf[ext].data
        head = pf[ext].header
        photflam = head['photflam']
        photplam = head['photplam']
        photbw   = head['photbw']
        
        pf = fits.open(wht_seg_str)
        #  inverse variance weight images
        wht = pf[ext].data

        noise = 1./np.sqrt(wht)

        flux_seg = np.sum(image)
        mag_seg = hst_phot(photflam, photplam, flux_seg, dn_err=None)


        # store values in a dict --> eventually write to a table...
        mdict_seg[obj][filt+'_mag'] = float(mag_seg)
        mdict_seg[obj][filt+'_pivot'] = float(photplam)
        mdict_seg[obj][filt+'_bandwidth'] = float(photbw)/2.

        mdict_seg[obj]["keys"]["mag"].append(filt+'_mag')
        mdict_seg[obj]["keys"]["pivot"].append(filt+'_pivot')
        mdict_seg[obj]["keys"]["bandwidth"].append(filt+'_bandwidth')



    fig = plt.figure(figsize=(8,6))
    p1 = fig.add_subplot(111)

    # dictionary and table keys
    mag_keys = mdict[obj]["keys"]["mag"]
    magerr_keys = mdict[obj]["keys"]["magerr"]
    pivot_keys = mdict[obj]["keys"]["pivot"]
    bandwidth_keys = mdict[obj]["keys"]["bandwidth"]
    
    
    mag_seg_keys = mdict_seg[obj]["keys"]["mag"]
    pivot_keys_seg = mdict_seg[obj]["keys"]["pivot"]
    bandwidth_keys_seg = mdict_seg[obj]["keys"]["bandwidth"]



    pivot = [mdict[obj][key] for key in pivot_keys]
    bandwidth = [mdict[obj][key] for key in bandwidth_keys]

    # Merino photometry (from dictionary)
    mag1 = [mdict[obj][key] for key in mag_keys]
    magerr1 = [mdict[obj][key] for key in magerr_keys]

    # Molino photometry (from table)
    mag2 = [mcat_row[key] for key in mag_keys]
    magerr2 = [mcat_row[key] for key in magerr_keys]

    # Watershed photometry
    mag3 = [mdict_seg[obj][key] for key in mag_seg_keys]


    # determine range for plotting (filtering bad mags)
    mag2 = np.array(mag2)
    mask = mag2 > -99
    mag2_mask = np.compress(mask,mag2,0)
    y0 = np.min(mag2_mask)
    y1 = np.max(mag2_mask)
    #print mag2
    #print mag2_mask
    #print y0
    #print y1


    # This section will direct the code to the HST filter throughputs files
    x_f105w, y_f105w, x_f110w, y_f110w, x_f125w, y_f125w, x_f140w, y_f140w,\
    x_f160w, y_f160w, x_f225w, y_f225w, x_f275w, y_f275w, x_f336w, y_f336w,\
    x_f390w, y_f390w, x_f435w, y_f435w, x_f475w, y_f475w, x_f606w, y_f606w,\
    x_f625w, y_f625w, x_f775w, y_f775w, x_f814w, y_f814w, x_f850lp, y_f850lp = HST_filt()



    p1.errorbar(pivot, mag1, yerr=magerr1, xerr=bandwidth, marker='o', linestyle='', label="Circular App.", zorder=2)
    p1.errorbar(pivot, mag2, yerr=magerr2, xerr=bandwidth, marker='*', linestyle='', label="Molino", zorder=2)
    p1.scatter(pivot, mag3, marker='o', label='Watershed', c='r', zorder=2)
    p1.set_title('CLASH id: ' + str(obj) + ', z = ' + str(z))
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


    plt.savefig(plot_dir + "/%s_%s.png" % (cluster,obj))
    #plt.show()


    # This section establishes the name of each column for the catalog
    Cluster, ID, RA, DEC, redshift = [],[],[],[],[]

    f105w_mag,f110w_mag,f125w_mag,f140w_mag,f160w_mag = [],[],[],[],[]
    f225w_mag,f275w_mag,f336w_mag,f390w_mag,f435w_mag = [],[],[],[],[]
    f475w_mag,f606w_mag,f625w_mag,f775w_mag,f814w_mag,f850lp_mag = [],[],[],[],[],[]

    f105w_err,f110w_err,f125w_err,f140w_err,f160w_err = [],[],[],[],[]
    f225w_err,f275w_err,f336w_err,f390w_err,f435w_err = [],[],[],[],[]
    f475w_err,f606w_err,f625w_err,f775w_err,f814w_err,f850lp_err = [],[],[],[],[],[]

    Cluster.append(cluster)
    ID.append(obj)
    RA.append(ra)
    DEC.append(dec)
    redshift.append(z)
    f105w_mag.append(round(mdict[obj]['f105w_mag'],3))
    f110w_mag.append(round(mdict[obj]['f110w_mag'],3))
    f125w_mag.append(round(mdict[obj]['f125w_mag'],3))
    f140w_mag.append(round(mdict[obj]['f140w_mag'],3))
    f160w_mag.append(round(mdict[obj]['f160w_mag'],3))
    f225w_mag.append(round(mdict[obj]['f225w_mag'],3))
    f275w_mag.append(round(mdict[obj]['f275w_mag'],3))
    f336w_mag.append(round(mdict[obj]['f336w_mag'],3))
    f390w_mag.append(round(mdict[obj]['f390w_mag'],3))
    f435w_mag.append(round(mdict[obj]['f435w_mag'],3))
    f475w_mag.append(round(mdict[obj]['f475w_mag'],3))
    f606w_mag.append(round(mdict[obj]['f606w_mag'],3))
    f625w_mag.append(round(mdict[obj]['f625w_mag'],3))
    f775w_mag.append(round(mdict[obj]['f775w_mag'],3))
    f814w_mag.append(round(mdict[obj]['f814w_mag'],3))
    f850lp_mag.append(round(mdict[obj]['f850lp_mag'],3))

    f105w_err.append(round(mdict[obj]['f105w_magerr'],3))
    f110w_err.append(round(mdict[obj]['f110w_magerr'],3))
    f125w_err.append(round(mdict[obj]['f125w_magerr'],3))
    f140w_err.append(round(mdict[obj]['f140w_magerr'],3))
    f160w_err.append(round(mdict[obj]['f160w_magerr'],3))
    f225w_err.append(round(mdict[obj]['f225w_magerr'],3))
    f275w_err.append(round(mdict[obj]['f275w_magerr'],3))
    f336w_err.append(round(mdict[obj]['f336w_magerr'],3))
    f390w_err.append(round(mdict[obj]['f390w_magerr'],3))
    f435w_err.append(round(mdict[obj]['f435w_magerr'],3))
    f475w_err.append(round(mdict[obj]['f475w_magerr'],3))
    f606w_err.append(round(mdict[obj]['f606w_magerr'],3))
    f625w_err.append(round(mdict[obj]['f625w_magerr'],3))
    f775w_err.append(round(mdict[obj]['f775w_magerr'],3))
    f814w_err.append(round(mdict[obj]['f814w_magerr'],3))
    f850lp_err.append(round(mdict[obj]['f850lp_magerr'],3))


    t = QTable([Cluster,ID,RA,DEC,redshift,f105w_mag,f105w_err,f110w_mag,f110w_err,f125w_mag,f125w_err,\
                f140w_mag,f140w_err,f160w_mag,f160w_err,f225w_mag,f225w_err,f275w_mag,f275w_err,\
                f336w_mag,f336w_err,f390w_mag,f390w_err,f435w_mag,f435w_err,f475w_mag,f475w_err,f606w_mag,f606w_err,\
                f625w_mag,f625w_err,f775w_mag,f775w_err,f814w_mag,f814w_err,f850lp_mag,f850lp_err],\
                names=('Cluster','ID','RA','DEC','redshift','f105w_mag','f105w_err','f110w_mag','f110w_err','f125w_mag','f125w_err',\
                'f140w_mag','f140w_err','f160w_mag','f160w_err','f225w_mag','f225w_err','f275w_mag','f275w_err',\
                'f336w_mag','f336w_err','f390w_mag','f390w_err','f435w_mag','f435w_err','f475w_mag','f475w_err','f606w_mag','f606w_err',\
                'f625w_mag','f625w_err','f775w_mag','f775w_err','f814w_mag','f814w_err','f850lp_mag','f850lp_err'))


print(t)
path = '/Users/brianmerino/Desktop/Apperture_Test/'
t.write(path+'a383_phot.csv', overwrite=True)



t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')
