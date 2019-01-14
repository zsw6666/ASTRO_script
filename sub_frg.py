import numpy as np
from astropy.io import fits
from sys import argv
import os
'''
This script used to subtract fringe 

objlist: science image list
'''

# path,fitslistname,frgname=argv[1:]
path='/Users/shiwuzhang/ASTRO/script'
os.chdir(path)
fitslistname='list'
frgname='fringe.fits'
os.chdir(path)
fitslist=open(fitslistname).readlines()
frgfits=fits.open(frgname)

for fitsname in fitslist:
    scifits=fits.open(fitsname[:-1])
    sub_frgfits=fits.HDUList()
    sub_frgfits.append(scifits[0])
    for i in range(1,len(scifits)):
        scifits[i].data-=frgfits[i-1].data
        sub_frgfits.append(scifits[i])

sub_frgfits.writeto('subfrg.fits',overwrite=True)

