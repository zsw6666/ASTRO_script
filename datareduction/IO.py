from astropy.io import fits
import numpy as np
import os


def txt_io(path,txtname,mark=None,celltype=None):
    '''
    read data from txt and return numpy array
    :param path: path of the txt
    :param celltype: data type for each column,for exaple:['i8','f8','S5']
    i8: int 8, f8: float 8, S5: string 5
    :param txtname: its name
    :param mark: how data in txt seperated(spaceline or ,)
    :return: numpy array of data in txt
    '''
    os.chdir(path)
    txt_file=np.genfromtxt(txtname,dtype=celltype,delimiter=mark)
    return txt_file


def GetImg(fits_name,ext=0):
    '''
    get image of the fits you can also appoint which extension
    :param fits_name: name of fits you should cd to the dictionary of the fits first
    :return: 2D numpy array which is the image
    '''
    fits_file=fits.open(fits_name,'readonly')
    data=fits_file[ext].data
    fits_file.close()
    return data

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