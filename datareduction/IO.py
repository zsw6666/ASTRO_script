from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


def txt_io(path,txtname):
    with open(path+'/'+txtname) as txtfile:
        containlist=[]
        for i in txtfile.readlines():
            containlist.append(i[:-1])
    return containlist

def Accessfits(fits_name,ext=0):
    '''
    get image of the fits you can also appoint which extension
    :param fits_name: name of fits you should cd to the dictionary of the fits first
    :return: 2D numpy array which is the image
    '''
    fits_file=fits.open(fits_name,'readonly')
    data=fits_file[ext].data
    header=fits_file[ext].header
    wcs=WCS(header)
    fits_file.close()
    return data,header,wcs

def UpdateImg(fits_name,new_img,ext=0):
    '''
    update the image
    :param fits_name: same with GetImg
    :param new_img: same with GetImg
    :param ext: same with UpdateImg
    :return: no return
    '''
    fits_file = fits.open(fits_name, 'update')
    fits_file[ext].data=new_img
    fits_file.close()
    return None