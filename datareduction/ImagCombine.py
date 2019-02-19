import numpy as np
from astropy.io import fits
from sys import argv
import os
'''
This script used to combine lbc image 

fitslist: list contain the images' name
path: path of the images and list
'''

path,fitslistname=argv[1:]
os.chdir(path)
print(path)
print(fitslistname)
fitslist=open(fitslistname).readlines()

objmasklist=[]
backgroundlist=[]

#mask the most bright source to make the fringe image better
for i in range(len(fitslist)):
      print('sextractor fits '+str(i)+'.....')
      sex1='sex '+fitslist[i][:-1]+' -c bok001.sex -SATUR_LEVEL 60000 -DETECT_MINAREA 5 -DETECT_MAXAREA 5000 -DETECT_THRESH 2.2 -ANALYSIS_THRESH 0.2 -FILTER Y -FILTER_NAME gauss_2.5_5x5.conv -CHECKIMAGE_TYPE -OBJECTS -CHECKIMAGE_NAME '+'obj_'+fitslist[i]
      objmasklist.append('obj_'+fitslist[i][:-1])
      os.system(sex1)
      sex2='sex '+fitslist[i][:-1]+' -c bok001.sex -SATUR_LEVEL 60000 -DETECT_MINAREA 5 -DETECT_MAXAREA 5000 -DETECT_THRESH 2.2 -ANALYSIS_THRESH 0.2 -FILTER Y -FILTER_NAME gauss_2.5_5x5.conv -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+'bk_'+fitslist[i]
      backgroundlist.append('bk_'+fitslist[i][:-1])
      os.system(sex2)
     
    





#make fringe image by median,because it adopt dither mode
outimg=fits.HDUList()
backimg=fits.HDUList()
for i in range(1,5):
    print('the '+str(i)+'ext.....')
    ext=[]
    back=[]
    for fitsi in range(len(objmasklist)):
	    print('the '+str(fitsi)+' fits......')
	    fitsfile=fits.open(objmasklist[fitsi])
	    backfile=fits.open(backgroundlist[fitsi])
	    ext.append(fitsfile[i].data)
	    back.append(backfile[i].data)
	    backfile.close()
	    fitsfile.close()
    ext=np.array(ext)
    back=np.array(back)
    new_ext=np.median(ext,axis=0)
    new_back=np.median(back,axis=0)
    outimg.append(fits.ImageHDU(new_ext))
    backimg.append(fits.ImageHDU(new_back))

print(len(outimg))
backimg.writeto('background.fits',overwrite=True)
outimg.writeto('fringe.fits',overwrite=True)
print('done')
os.system('rm obj*fits')




