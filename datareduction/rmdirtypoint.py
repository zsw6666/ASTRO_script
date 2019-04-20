from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
from . import IO

'''
this script is used to remove the dirty points in iamges,
dirty points are some structrues in image whose size is not large and value is 
low then the background. it's typical shape is circle
written by SW Zhang
2019.04.20
'''





def ArrayStat(array):
    '''
    calculate the statistical value of the array
    :param array:
    :return: statistical value of this array
    '''
    average=np.mean(array)
    std=np.std(array)
    return average,std


def RmlargeValue(array):
    '''
    remove pixels whose value is big enough to interfere your result
    :param array:
    :return: new array the value of some large-value pixels is reset to the average
    '''
    array_mean,array_std=ArrayStat(array)
    array[np.where(array>=array_mean+3.*array_std)]=array_mean
    return array

def RepeatRmLV(array,n):
    '''
    repeat RmlargeValue n times to ensure that the large-value pixels is cleaned
    :param array:
    :param n: times you want to repeat RmlargeValue
    :return: new array without large-value pixels
    '''
    for i in range(n):
        array=RmlargeValue(array)
    return array

def LocateDP(array,sigma):
    '''
    locate the pysical coordinates of dirty points
    :param array: array which might contain dirty points
    :param sigma: value you use to determine what kind of pixel is dirty point
    :return: poision in this array of dirty point
    '''
    array_mean,array_std=ArrayStat(array)
    position=np.where(array<=array_mean-sigma*array_std)
    return position

def CorrectDirtyPoint(img):
    '''
    this function use the functions above to remove the dirty points in the input image
    the idea is simple, size of dirty point is ~10*10 pixel, select 10 raws in a team,
    calculate the median along the direction of raw, if there is a valley and the width of
    this valley is larger than 10 pixels, identifying this as dirty point and substitude it's
    value with surrounding pixels.
    :param img: original image with dirty point
    :return: image dirty points of which have been masked
    '''
    i=9
    while i<=len(img)-1:
        line=np.median(img[i-9:i,:],axis=0)
        # line_min=np.where(line==np.min(line))
        line=RepeatRmLV(line,3)
        dp_position=LocateDP(line,3)
        if len(dp_position[0])>1:#if there's only one pixel in this line is identified as dp, it can't be
                                 #dp, as a result, move on to the next line, instead, find the general location of it.
            start = np.min(dp_position[0])
            end = np.max(dp_position[0])
        else:
            i+=9
            continue
        if (end-start)>=10 and (end-start)<=30 and len(dp_position[0])>=10:# size of dp cannot be to large and too small
            img[i-12:i+5,start-2:end+3]=np.median(line)
            print('remove dirty point')
        else:
            pass
        i+=9
    return img

def Run(path,fitsname):
    '''
    run the process to remove dirty points in images
    :param path: dictionary which contain images
    :param fitsname: name of fits
    :return: none
    '''
    os.chdir(path)
    img=CorrectDirtyPoint(IO.GetImg(fitsname))
    IO.UpdateImg(fitsname,img)
    return None


if __name__=='__main__':
    Run('/Users/shiwuzhang/ASTRO','lbcb180515_LBT-2018A-C1657-1_02_03.fits')

