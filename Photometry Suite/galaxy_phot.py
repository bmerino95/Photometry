#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 11:28:29 2021

NOTES
a611 does not have f625w
macs0744 does not have f625w and has 775w with ACS and wfc3uvis
macs1423 does not have f625w
"""

"""
2021/03/18 - This version of the code performs photometry over segmented images of the entire galaxy - BMM
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
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/clump.cfg"
bkgd_subtract = 1   # use background subtracted frames for photometry
ext = 0
####################################


# read-in cfg file with paramters 
config = configparser.ConfigParser()
config.read(cfg)
img_dir = config.get('CONFIG','img_dir')
sub_img_dir_drz = config.get('CONFIG','sub_img_dir_drz')
#sub_img_dir_wht = config.get('CONFIG','sub_img_dir_wht')
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

# This cfg file will contain the variables that each code needs. 
cfg = "/Users/brianmerino/Desktop/Apperture_Test/Code/general.cfg"

config = configparser.ConfigParser()
config.read(cfg)

cluster = config.get('CONFIG','cluster')
ID = int(config.get('CONFIG','Id'))


#cluster = 'macs0329'
#ID = '1903'
z = float(config.get('CONFIG','z'))
counter = -1 #In the region files, the galactic coords are listed last
drop_uvis = 1
Area = -99 # don't have these values yet
mag,magerr = {},{}
pivot = {}
bandwidth = {}



#path = '/Users/brianmerino/Desktop/Apperture_Test/a383_clump%s_phot.csv'%(counter)
path = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/galaxy_phot.csv'%(cluster,ID)

if os.path.exists(path):
    os.remove(path)
    
    # Set up header for catalog
    print('#1 \t cluster \n#2 \t id \n#3 \t obj \n#4\t RA \n#5 \t DEC \n#6 \t area \n\
#7 \t zb \n#8 \t f225w_mag \n#9 \t f225w_magerr \n#10 \t f275w_mag \n\
#11 \t f275w_magerr \n#12 \t f336w_mag \n#13 \t f336w_magerr \n#14 \t f390w_mag \n\
#15 \t f390w_magerr \n#16 \t f435w_mag \n#17 \t f435w_magerr \n#18 \t f475w_mag \n\
#19 \t f475w_magerr \n#20 \t f606w_mag \n#21 \t f606w_magerr \n#22 \t f625w_mag \n\
#23 \t f625w_magerr \n#24 \t f775w_mag \n#25 \t f775w_magerr \n#26 \t f814w_mag \n\
#27 \t f814w_magerr \n#28 \t f850lp_mag \n#29 \t f850lp_magerr \n#30 \t f105w_mag \n\
#31 \t f105w_magerr \n#32 \t f110w_mag \n#33 \t f110w_magerr \n#34 \t f125w_mag \n\
#35 \t f125w_magerr \n#36 \t f140w_mag \n#37 \t f140w_magerr \n#38 \t f160w_mag \n\
#39 \t f160w_magerr', file=open(path, "a"))

else:
    # Set up header for catalog
    print('#1 \t cluster \n#2 \t id \n#3 \t obj \n#4\t RA \n#5 \t DEC \n#6 \t area \n\
#7 \t zb \n#8 \t f225w_mag \n#9 \t f225w_magerr \n#10 \t f275w_mag \n\
#11 \t f275w_magerr \n#12 \t f336w_mag \n#13 \t f336w_magerr \n#14 \t f390w_mag \n\
#15 \t f390w_magerr \n#16 \t f435w_mag \n#17 \t f435w_magerr \n#18 \t f475w_mag \n\
#19 \t f475w_magerr \n#20 \t f606w_mag \n#21 \t f606w_magerr \n#22 \t f625w_mag \n\
#23 \t f625w_magerr \n#24 \t f775w_mag \n#25 \t f775w_magerr \n#26 \t f814w_mag \n\
#27 \t f814w_magerr \n#28 \t f850lp_mag \n#29 \t f850lp_magerr \n#30 \t f105w_mag \n\
#31 \t f105w_magerr \n#32 \t f110w_mag \n#33 \t f110w_magerr \n#34 \t f125w_mag \n\
#35 \t f125w_magerr \n#36 \t f140w_mag \n#37 \t f140w_magerr \n#38 \t f160w_mag \n\
#39 \t f160w_magerr', file=open(path, "a"))



def phot(image_drz,noise,photflam,photplam):
    flux = np.sum(image_drz)
    
    flux_err = noise
    print(flux_err)
    
    if noise == 0:
        SN = 0
    else:
        SN = flux/noise
    
    if item == 'f225w' or item == 'f275w' or item == 'f336w' or item == 'f390w':
        if drop_uvis:
            mag = -99
            magerr = -99
        
        elif SN < 8:
            mag = -99
            magerr= -99
        
        else:
            mag, magerr = hst_phot(photflam, photplam, flux, dn_err=flux_err) # AB magnitude
    
    else:
        if SN < 3:
            mag = -99
            magerr = -99
    
        else:
            mag, magerr = hst_phot(photflam, photplam, flux, dn_err=flux_err) # AB magnitude
    
    return round(mag,4),round(magerr,4)


if cluster == 'a611' or cluster =='macs0717' or cluster =='macs0744':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
               'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
               'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

elif cluster == 'macs1423':
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w',\
                   'f390w', 'f435w', 'f475w', 'f606w', 'f775w', 'f814w', 'f850lp']

else:
    filter_list = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w',\
               'f225w', 'f275w', 'f336w', 'f390w', 'f435w',\
               'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp']



region_path = '/Users/brianmerino/Desktop/Individual_Clumps/%s/'%(cluster)
region_name = region_path+'%s.reg'%(ID)
#region_name = '/Users/brianmerino/Desktop/Apperture_Test/%s-%s/clump.reg'%(cluster,ID)
r = pyregion.open(region_name)


for item in filter_list:
    #pf = fits.open(sub_img_dir_drz+'clump_'+str(counter)+'_'+item+'.fits')
    pf = fits.open(sub_img_dir_drz+'%s-%s/masked/galaxy_'%(cluster,ID)+item+'.fits')
    image = pf[ext].data
    head = pf[ext].header
    photflam = head['photflam']
    photplam = head['photplam']
    photbw   = head['photbw']
    
    pf_wht = fits.open(sub_img_dir_drz+'%s-%s/masked/galaxy_'%(cluster,ID)+item+'_wht.fits')
    #  inverse variance weight images
    wht = pf_wht[ext].data
    
    quotient = 1/wht
    
    for i in range(0,quotient.shape[0]):
        for j in range(0,quotient.shape[1]):
            if np.isinf(quotient[i][j]):
                quotient[i][j]=0
    
    SUM = np.sum(quotient)
    
    noise = np.sqrt(SUM)
    
    wcs = WCS(head)
    wcs.sip = None
    
    mag[item+'_mag'],magerr[item+'_magerr']=phot(image,noise,photflam,photplam)
    pivot[item]=float(photplam)
    bandwidth[item]=float(photbw)/2
    
    #counter +=1


fig = plt.figure(figsize=(8,6))
p1 = fig.add_subplot(111)


# Merino photometry (from dictionary)
mag1 = [mag[obj] for obj in mag]
magerr1 = [magerr[obj] for obj in magerr]
pivots = [pivot[obj] for obj in pivot]
bandwidths = [bandwidth[obj] for obj in bandwidth]



# determine range for plotting (filtering bad mags)
mag1 = np.array(mag1)

for i in range(0,len(mag1)):
    if np.isinf(mag1[i]):
        mag1[i]=-99

#mask = mag2 > -99
mask = mag1 > -50
mag1_mask = np.compress(mask,mag1,0)
y0 = np.min(mag1_mask)
y1 = np.max(mag1_mask)
#print mag2
#print mag2_mask
#print y0
#print y1


# This section will direct the code to the HST filter throughputs files
x_f105w, y_f105w, x_f110w, y_f110w, x_f125w, y_f125w, x_f140w, y_f140w,\
x_f160w, y_f160w, x_f225w, y_f225w, x_f275w, y_f275w, x_f336w, y_f336w,\
x_f390w, y_f390w, x_f435w, y_f435w, x_f475w, y_f475w, x_f606w, y_f606w,\
x_f625w, y_f625w, x_f775w, y_f775w, x_f814w, y_f814w, x_f850lp, y_f850lp = HST_filt()



p1.errorbar(pivots, mag1, xerr=bandwidths, yerr=magerr1, marker='o', linestyle='', label="Merino", zorder=2)
#p1.set_title('CLASH id: ' + str(cluster)+  ' ' + str(ID) + ' Clump '+str(counter)+', z = ' + str(z))
p1.set_title('CLASH id: ' + str(cluster)+  ' ' + str(ID) + ' Galaxy,  z = ' + str(z))
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
plt.savefig(plot_dir + '/%s-%s/galaxy.png'%(cluster,ID))
plt.show()


# This section establishes the name of each column for the catalog
Cluster, CLASHID, RA, DEC, area, redshift = [],[],[],[],[],[]

f105w_mag,f110w_mag,f125w_mag,f140w_mag,f160w_mag = [],[],[],[],[]
f225w_mag,f275w_mag,f336w_mag,f390w_mag,f435w_mag = [],[],[],[],[]
f475w_mag,f606w_mag,f625w_mag,f775w_mag,f814w_mag,f850lp_mag = [],[],[],[],[],[]

f105w_err,f110w_err,f125w_err,f140w_err,f160w_err = [],[],[],[],[]
f225w_err,f275w_err,f336w_err,f390w_err,f435w_err = [],[],[],[],[]
f475w_err,f606w_err,f625w_err,f775w_err,f814w_err,f850lp_err = [],[],[],[],[],[]

Cluster.append(cluster)
CLASHID.append(ID)
RA.append(r[counter].coord_list[0])
DEC.append(r[counter].coord_list[1])
area.append(Area)
redshift.append(z)
f105w_mag.append(round(mag['f105w_mag'],3))
f110w_mag.append(round(mag['f110w_mag'],3))
f125w_mag.append(round(mag['f125w_mag'],3))
f140w_mag.append(round(mag['f140w_mag'],3))
f160w_mag.append(round(mag['f160w_mag'],3))
f225w_mag.append(round(mag['f225w_mag'],3))
f275w_mag.append(round(mag['f275w_mag'],3))
f336w_mag.append(round(mag['f336w_mag'],3))
f390w_mag.append(round(mag['f390w_mag'],3))
f435w_mag.append(round(mag['f435w_mag'],3))
f475w_mag.append(round(mag['f475w_mag'],3))
f606w_mag.append(round(mag['f606w_mag'],3))

if cluster=='a611' or cluster=='macs0717' or cluster=='macs1423' or cluster =='macs0744':
    f625w_mag.append(-99)
else:
    f625w_mag.append(round(mag['f625w_mag'],3))

f775w_mag.append(round(mag['f775w_mag'],3))

if cluster=='macs1423':
    f814w_mag.append(-99)
else:
    f814w_mag.append(round(mag['f814w_mag'],3))

f850lp_mag.append(round(mag['f850lp_mag'],3))


f105w_err.append(round(magerr['f105w_magerr'],3))
f110w_err.append(round(magerr['f110w_magerr'],3))
f125w_err.append(round(magerr['f125w_magerr'],3))
f140w_err.append(round(magerr['f140w_magerr'],3))
f160w_err.append(round(magerr['f160w_magerr'],3))
f225w_err.append(round(magerr['f225w_magerr'],3))
f275w_err.append(round(magerr['f275w_magerr'],3))
f336w_err.append(round(magerr['f336w_magerr'],3))
f390w_err.append(round(magerr['f390w_magerr'],3))
f435w_err.append(round(magerr['f435w_magerr'],3))
f475w_err.append(round(magerr['f475w_magerr'],3))
f606w_err.append(round(magerr['f606w_magerr'],3))

if cluster=='a611' or cluster=='macs0717' or cluster=='macs1423' or cluster =='macs0744':
    f625w_err.append(-99)
else:
    f625w_err.append(round(magerr['f625w_magerr'],3))

f775w_err.append(round(magerr['f775w_magerr'],3))

if cluster=='macs1423':
    f814w_err.append(-99)
else:
    f814w_err.append(round(magerr['f814w_magerr'],3))

f850lp_err.append(round(magerr['f850lp_magerr'],3))



print(Cluster[0],' ',CLASHID[0],' ',counter,' ',RA[0],' ',DEC[0],' ',area[0],' ',redshift[0],' ',\
      f225w_mag[0],' ',f225w_err[0],' ',f275w_mag[0],' ',f275w_err[0],' ',f336w_mag[0],' ',f336w_err[0],' ',f390w_mag[0],' ',f390w_err[0],' ',\
      f435w_mag[0],' ',f435w_err[0],' ',f475w_mag[0],' ',f475w_err[0],' ',f606w_mag[0],' ',f606w_err[0],' ',\
      f625w_mag[0],' ',f625w_err[0],' ',f775w_mag[0],' ',f775w_err[0],' ',f814w_mag[0],' ',f814w_err[0],' ',\
      f850lp_mag[0],' ',f850lp_err[0],' ',f105w_mag[0],' ',f105w_err[0],' ',f110w_mag[0],' ',f110w_err[0],' ',\
      f125w_mag[0],' ',f125w_err[0],' ',f140w_mag[0],' ',f140w_err[0],' ',f160w_mag[0],' ',f160w_err[0], file=open(path, "a"))



t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')
