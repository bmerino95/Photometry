#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 18:25:01 2020

@author: brianmerino
"""
import numpy as np
import sys
import astropy
from photutils import Background2D
import time

t0 = time.time()

#NOTES
#a611 does not have f625w
#macs0744 does not have f625w and has 775w with ACS and wfc3uvis
#macs1423 does not have f625w and f814w doesn't work

#Filters and Cluster that are being used
filter_list = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f625w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']
filter_list_a611 = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']
macs1423_filter_list = ['f225w','f275w','f336w','f390w','f435w','f475w','f606w','f775w','f850lp','f105w','f110w','f125w','f140w','f160w']

cluster,iden,area = np.loadtxt('/Users/brianmerino/Desktop/CLUSTERS/master_catalog.txt',delimiter='\t',usecols=[2,3,4],skiprows=1,unpack=True)
path = "/Volumes/data7/Research/Thumbnails/"

#Create dictionaries
filters = {}
bkg     = {}
bs      = {}

for j in range(0,len(cluster)):
    test_string = cluster[j][1:-1]
    
    if cluster[j] == ' a611 ' or cluster[j] == ' macs0744 ':
        for i in range(0,len(filter_list_a611)):
            filters['%s'%(filter_list_a611[i])] = str(path)+'%s_%s_%s.fits'%(test_string,int(iden[j]),filter_list_a611[i])
            
            pf = astropy.io.fits.open(filters['%s'%(str(filter_list_a611[i]))])
            head = pf[0].header
            data = pf[0].data
            bkg['bkg%s'%str(filter_list_a611[i])] = Background2D(data, (int(head['CRPIX1']/2),int(head['CRPIX2']/2)))
            bs['%s_bs'%str(filter_list_a611[i])] =data - bkg['bkg%s'%str(filter_list_a611[i])].background
            hdu1 = astropy.io.fits.PrimaryHDU(bkg['bkg%s'%str(filter_list_a611[i])].background,header=head)
            hdu1.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s_bkg.fits'%(test_string,int(iden[j]),filter_list_a611[i]),overwrite=True)
            hdu2 = astropy.io.fits.PrimaryHDU(bs['%s_bs'%str(filter_list_a611[i])],header=head)
            hdu2.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[j]),filter_list_a611[i]),overwrite=True)
            print('Finished  %s %s %s' %(cluster[j],int(iden[j]),filter_list_a611[i]))

    elif cluster[j] == ' macs1423 ':
        for i in range(0,len(macs1423_filter_list)):
            filters['%s'%(macs1423_filter_list[i])] = str(path)+'%s_%s_%s.fits'%(test_string,int(iden[j]),macs1423_filter_list[i])
            
            pf = astropy.io.fits.open(filters['%s'%(str(macs1423_filter_list[i]))])
            head = pf[0].header
            data = pf[0].data
            bkg['bkg%s'%str(macs1423_filter_list[i])] = Background2D(data, (int(head['CRPIX1']/2),int(head['CRPIX2']/2)))
            bs['%s_bs'%str(macs1423_filter_list[i])] =data - bkg['bkg%s'%str(macs1423_filter_list[i])].background
            hdu1 = astropy.io.fits.PrimaryHDU(bkg['bkg%s'%str(macs1423_filter_list[i])].background,header=head)
            hdu1.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s_bkg.fits'%(test_string,int(iden[j]),macs1423_filter_list[i]),overwrite=True)
            hdu2 = astropy.io.fits.PrimaryHDU(bs['%s_bs'%str(macs1423_filter_list[i])],header=head)
            hdu2.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[j]),macs1423_filter_list[i]),overwrite=True)
            print('Finished  %s %s %s' %(cluster[j],int(iden[j]),macs1423_filter_list[i]))

    else:
        for i in range(0,len(filter_list)):
            filters['%s'%(filter_list[i])] = str(path)+'%s_%s_%s.fits'%(test_string,int(iden[j]),filter_list[i])

            pf = astropy.io.fits.open(filters['%s'%(str(filter_list[i]))])
            head = pf[0].header
            data = pf[0].data
            bkg['bkg%s'%str(filter_list[i])] = Background2D(data, (int(head['CRPIX1']/2),int(head['CRPIX2']/2)))
            bs['%s_bs'%str(filter_list[i])] =data - bkg['bkg%s'%str(filter_list[i])].background
            hdu1 = astropy.io.fits.PrimaryHDU(bkg['bkg%s'%str(filter_list[i])].background,header=head)
            hdu1.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s_bkg.fits'%(test_string,int(iden[j]),filter_list[i]),overwrite=True)
            hdu2 = astropy.io.fits.PrimaryHDU(bs['%s_bs'%str(filter_list[i])],header=head)
            hdu2.writeto('/Volumes/Home_Drive/Research/Subtracted_Thumbnails/%s_%s_%s.fits'%(test_string,int(iden[j]),filter_list[i]),overwrite=True)
            print('Finished  %s %s %s' %(cluster[j],int(iden[j]),filter_list[i]))

t1 = time.time()
total = t1-t0
total_min = total/60
print()
print(total, 'seconds')
print(total_min, 'minutes')

