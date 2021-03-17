#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:17:22 2020

@author: brianmerino

This is an updated version of my catalog reader. It reads in all of CLASH's WFC3-IR catalogs
and 19 region files that contain the CLASH ids of the clumpy galaxies that I found by eye.
"""

"""
3/15/2021 - I have updated this code so that it creates a master catalog of CLASh photometry, a
Molino photometry (if the target is included in Molino's dataset'), and an extra catalog that 
states which objects are in both datasets. The extra catalog includes a flag and both ids so
that I can identify it quickly. -BMM
"""

from astropy.table import Table
import numpy as np
import time

t0 = time.time()

members = ['a209','a383','a611','macs0329','macs0416','macs0429','macs0717','macs0744','macs1115','macs1149','macs1206','macs1311','macs1423','macs1720','macs1931','macs2129','ms2137','rxj1347','rxj1532','rxj2129']
Clusters = ['a209','a383','a611','macs0329','macs0416','macs0429','macs0717','macs0744','macs1115','macs1149','macs1206','macs1311','macs1423','macs1720','macs1931','macs2129','ms2137','rxj1347','rxj1532','rxj2129']

for i in range(0,len(members)):
    Clusters[i] = Table.read('/Users/brianmerino/Desktop/CLUSTERS/'+str(members[i])+'/catalogs/hst/hlsp_clash_hst_ir_'+str(members[i])+'_cat.txt',format = 'ascii')

# Read in master region file
#test = Table.read('/Users/brianmerino/Desktop/CLUSTERS/ID_and_cluster.txt',format = 'ascii')
#print(test)

save_path = "/Users/brianmerino/Desktop/Regions/Clumpy_ids/"
final_path_clash  = save_path + 'master_catalog_clash.txt'
final_path_molino = save_path + 'master_catalog_molino.txt'
final_path_both   = save_path + 'master_catalog_clash_and_molino.txt'

count_clash = 0
count_molino = 0

print('#1 \t RA \n#2 \t Dec \n#3 \t clusterName \n#4 \t CLASHID \n\
#5 \t area \n#6 \t zb_1 \n#7 \t F225W_WFC3_PHOTOZ \n#8 \t dF225W_WFC3_PHOTOZ \n\
#9 \t F275W_WFC3_PHOTOZ \n#10 \t dF275W_WFC3_PHOTOZ \n#11 \t F336W_WFC3_PHOTOZ \n#12 \t dF336W_WFC3_PHOTOZ \n\
#13 \t F390W_WFC3_PHOTOZ \n#14 \t dF390W_WFC3_PHOTOZ \n#15 \t F435W_ACS_PHOTOZ \n#16 \t dF435W_ACS_PHOTOZ \n\
#17 \t F475W_ACS_PHOTOZ \n#18 \t dF475W_ACS_PHOTOZ \n#19 \t F606W_ACS_PHOTOZ \n#20 \t dF606W_ACS_PHOTOZ \n\
#21 \t F625W_ACS_PHOTOZ \n#22 \t dF625W_ACS_PHOTOZ \n#23 \t F775W_ACS_PHOTOZ \n#24 \t dF775W_ACS_PHOTOZ \n\
#25 \t F814W_ACS_PHOTOZ \n#26 \t dF814W_ACS_PHOTOZ \n#27 \t F850LP_ACS_PHOTOZ \n#28 \t dF850LP_ACS_PHOTOZ \n\
#29 \t F105W_WFC3_PHOTOZ \n#30 \t dF105W_WFC3_PHOTOZ \n#31 \t F110W_WFC3_PHOTOZ \n#32 \t dF110W_WFC3_PHOTOZ \n\
#33 \t F125W_WFC3_PHOTOZ \n#34 \t dF125W_WFC3_PHOTOZ \n#35 \t F140W_WFC3_PHOTOZ \n#36 \t dF140W_WFC3_PHOTOZ \n\
#37 \t F160W_WFC3_PHOTOZ \n#38 \t dF160W_WFC3_PHOTOZ', file=open(final_path_molino, "a"))

print('#1 \t RA \n#2 \t Dec \n#3 \t cluster \n#4 \t id \n\
#5 \t area \n#6 \t zb \n#7 \t f225w_mag \n#8 \t f225w_magerr \n\
#9 \t f275w_mag \n#10 \t f275w_magerr \n#11 \t f336w_mag \n#12 \t f336w_magerr \n\
#13 \t f390w_mag \n#14 \t f390w_magerr \n#15 \t f435w_mag \n#16 \t f435w_magerr \n\
#17 \t f475w_mag \n#18 \t f475w_magerr \n#19 \t f606w_mag \n#20 \t f606w_magerr \n\
#21 \t f625w_mag \n#22 \t f625w_magerr \n#23 \t f775w_mag \n#24 \t f775w_magerr \n\
#25 \t f814w_mag \n#26 \t f814w_magerr \n#27 \t f850lp_mag \n#28 \t f850lp_magerr \n\
#29 \t f105w_mag \n#30 \t f105w_magerr \n#31 \t f110w_mag \n#32 \t f110w_magerr \n\
#33 \t f125w_mag \n#34 \t f125w_magerr \n#35 \t f140w_mag \n#36 \t f140w_magerr \n\
#37 \t f160w_mag \n#38 \t f160w_magerr \n#39 \t CLASH&Molino \n#40 \t MolinoID', file=open(final_path_both, "a"))

print('#1 \t RA \n#2 \t Dec \n#3 \t cluster \n#4 \t id \n\
#5 \t area \n#6 \t zb \n#7 \t f225w_mag \n#8 \t f225w_magerr \n\
#9 \t f275w_mag \n#10 \t f275w_magerr \n#11 \t f336w_mag \n#12 \t f336w_magerr \n\
#13 \t f390w_mag \n#14 \t f390w_magerr \n#15 \t f435w_mag \n#16 \t f435w_magerr \n\
#17 \t f475w_mag \n#18 \t f475w_magerr \n#19 \t f606w_mag \n#20 \t f606w_magerr \n\
#21 \t f625w_mag \n#22 \t f625w_magerr \n#23 \t f775w_mag \n#24 \t f775w_magerr \n\
#25 \t f814w_mag \n#26 \t f814w_magerr\n#27 \t f850lp_mag \n#28 \t f850lp_magerr \n\
#29 \t f105w_mag \n#30 \t f105w_magerr \n#31 \t f110w_mag \n#32 \t f110w_magerr \n\
#33 \t f125w_mag \n#34 \t f125w_magerr \n#35 \t f140w_mag \n#36 \t f140w_magerr \n\
#37 \t f160w_mag \n#38 \t f160w_magerr', file=open(final_path_clash, "a"))

# Create for loop to extract catalog information for each object
for i in range(0,len(members)):
    # Read in text files that contain the CLASH ids of the clumpy galaxies
    new_test = np.loadtxt(save_path+'%s_clumpy.reg'%(members[i]))
    test = []
    
    for j in range(0,len(new_test)):
        test.append(int(new_test[j]))


    # This section creates the CLASH sub catalog
    for m in range(0,len(test)):
        k = i
        
        # This section creates the Molino sub catalog
        path_start = '/Users/brianmerino/Desktop/CLUSTERS/'
        path_middle = '/catalogs/molino/hlsp_clash_hst_ir_'
        path_end = '_cat-molino.txt'
        
        group = members[k]
        
        cat_path = path_start + group + path_middle + group + path_end
        molino = Table.read(cat_path,format='ascii.sextractor')
        
        ra,dec = round(Clusters[k]['RA'][test[m]-1],4),round(Clusters[k]['Dec'][test[m]-1],4)
        
        for j in range(0,len(molino['RA'])):
            if ra == molino['RA'][j] and dec == molino['Dec'][j]:
                
                print(molino['RA'][j],' ',molino['Dec'][j],' ',molino['clusterName'][j],' ',molino['CLASHID'][j],' ',\
                      molino['area'][j],' ',molino['zb_1'][j],' ',molino['F225W_WFC3_PHOTOZ'][j],' ',molino['dF225W_WFC3_PHOTOZ'][j],' ',\
                      molino['F275W_WFC3_PHOTOZ'][j],' ',molino['dF275W_WFC3_PHOTOZ'][j],' ',molino['F336W_WFC3_PHOTOZ'][j],' ',molino['dF336W_WFC3_PHOTOZ'][j],' ',\
                      molino['F390W_WFC3_PHOTOZ'][j],' ',molino['dF390W_WFC3_PHOTOZ'][j],' ',molino['F435W_ACS_PHOTOZ'][j],' ',molino['dF435W_ACS_PHOTOZ'][j],' ',\
                      molino['F475W_ACS_PHOTOZ'][j],' ',molino['dF475W_ACS_PHOTOZ'][j],' ',molino['F606W_ACS_PHOTOZ'][j],' ',molino['dF606W_ACS_PHOTOZ'][j],' ',\
                      molino['F625W_ACS_PHOTOZ'][j],' ',molino['dF625W_ACS_PHOTOZ'][j],' ',molino['F775W_ACS_PHOTOZ'][j],' ',molino['dF775W_ACS_PHOTOZ'][j],' ',\
                      molino['F814W_ACS_PHOTOZ'][j],' ',molino['dF814W_ACS_PHOTOZ'][j],' ',molino['F850LP_ACS_PHOTOZ'][j],' ',molino['dF850LP_ACS_PHOTOZ'][j],' ',\
                      molino['F105W_WFC3_PHOTOZ'][j],' ',molino['dF105W_WFC3_PHOTOZ'][j],' ',molino['F110W_WFC3_PHOTOZ'][j],' ',molino['dF110W_WFC3_PHOTOZ'][j],' ',\
                      molino['F125W_WFC3_PHOTOZ'][j],' ',molino['dF125W_WFC3_PHOTOZ'][j],' ',molino['F140W_WFC3_PHOTOZ'][j],' ',molino['dF140W_WFC3_PHOTOZ'][j],' ',\
                      molino['F160W_WFC3_PHOTOZ'][j],' ',molino['dF160W_WFC3_PHOTOZ'][j], file=open(final_path_molino, "a"))
                
                print(Clusters[k]['RA'][test[m]-1],' ', Clusters[k]['Dec'][test[m]-1],' ', members[k],' ',\
                      Clusters[k]['id'][test[m]-1],' ',Clusters[k]['area'][test[m]-1],' ',Clusters[k]['zb'][test[m]-1],' ',\
                      Clusters[k]['f225w_mag'][test[m]-1],' ', Clusters[k]['f225w_magerr'][test[m]-1],' ', Clusters[k]['f275w_mag'][test[m]-1],' ',\
                      Clusters[k]['f275w_magerr'][test[m]-1],' ', Clusters[k]['f336w_mag'][test[m]-1],' ', Clusters[k]['f336w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f390w_mag'][test[m]-1],' ', Clusters[k]['f390w_magerr'][test[m]-1],' ', Clusters[k]['f435w_mag'][test[m]-1],' ',\
                      Clusters[k]['f435w_magerr'][test[m]-1],' ', Clusters[k]['f475w_mag'][test[m]-1],' ', Clusters[k]['f475w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f606w_mag'][test[m]-1],' ',\
                      Clusters[k]['f606w_magerr'][test[m]-1],' ', Clusters[k]['f625w_mag'][test[m]-1],' ', Clusters[k]['f625w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f775w_mag'][test[m]-1],' ', Clusters[k]['f775w_magerr'][test[m]-1],' ', Clusters[k]['f814w_mag'][test[m]-1],' ',\
                      Clusters[k]['f814w_magerr'][test[m]-1],' ',  Clusters[k]['f850lp_mag'][test[m]-1],' ', Clusters[k]['f850lp_magerr'][test[m]-1],' ',\
                      Clusters[k]['f105w_mag'][test[m]-1],' ', Clusters[k]['f105w_magerr'][test[m]-1],' ', Clusters[k]['f110w_mag'][test[m]-1],' ',\
                      Clusters[k]['f110w_magerr'][test[m]-1],' ', Clusters[k]['f125w_mag'][test[m]-1],' ', Clusters[k]['f125w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f140w_mag'][test[m]-1],' ', Clusters[k]['f140w_magerr'][test[m]-1],' ', Clusters[k]['f160w_mag'][test[m]-1],' ',\
                      Clusters[k]['f160w_magerr'][test[m]-1],' ', 1, ' ', molino['CLASHID'][j],file=open(final_path_both, "a"))
                    
                count_molino+=1
                
                
        print(Clusters[k]['RA'][test[m]-1],' ', Clusters[k]['Dec'][test[m]-1],' ', members[k],' ',\
                      Clusters[k]['id'][test[m]-1],' ',Clusters[k]['area'][test[m]-1],' ',Clusters[k]['zb'][test[m]-1],' ',\
                      Clusters[k]['f225w_mag'][test[m]-1],' ', Clusters[k]['f225w_magerr'][test[m]-1],' ', Clusters[k]['f275w_mag'][test[m]-1],' ',\
                      Clusters[k]['f275w_magerr'][test[m]-1],' ', Clusters[k]['f336w_mag'][test[m]-1],' ', Clusters[k]['f336w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f390w_mag'][test[m]-1],' ', Clusters[k]['f390w_magerr'][test[m]-1],' ', Clusters[k]['f435w_mag'][test[m]-1],' ',\
                      Clusters[k]['f435w_magerr'][test[m]-1],' ', Clusters[k]['f475w_mag'][test[m]-1],' ', Clusters[k]['f475w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f606w_mag'][test[m]-1],' ',\
                      Clusters[k]['f606w_magerr'][test[m]-1],' ', Clusters[k]['f625w_mag'][test[m]-1],' ', Clusters[k]['f625w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f775w_mag'][test[m]-1],' ', Clusters[k]['f775w_magerr'][test[m]-1],' ', Clusters[k]['f814w_mag'][test[m]-1],' ',\
                      Clusters[k]['f814w_magerr'][test[m]-1],' ',  Clusters[k]['f850lp_mag'][test[m]-1],' ', Clusters[k]['f850lp_magerr'][test[m]-1],' ',\
                      Clusters[k]['f105w_mag'][test[m]-1],' ', Clusters[k]['f105w_magerr'][test[m]-1],' ', Clusters[k]['f110w_mag'][test[m]-1],' ',\
                      Clusters[k]['f110w_magerr'][test[m]-1],' ', Clusters[k]['f125w_mag'][test[m]-1],' ', Clusters[k]['f125w_magerr'][test[m]-1],' ',\
                      Clusters[k]['f140w_mag'][test[m]-1],' ', Clusters[k]['f140w_magerr'][test[m]-1],' ', Clusters[k]['f160w_mag'][test[m]-1],' ',\
                      Clusters[k]['f160w_magerr'][test[m]-1],file=open(final_path_clash, "a"))
                    
        count_clash+=1

    print('Finishing %s'%(members[k]))
    print('Number of Molino Objects:', count_molino)
    print('Number of CLASH objects :', count_clash)
    print()

print('Done')

t1 = time.time()
total = t1-t0
total_min = total/60.
print()
print(total, 'seconds')
print(total_min, 'minutes')