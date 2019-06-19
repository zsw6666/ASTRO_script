'''
this script is used to produce the combined bias image from all bias.fits
created by Shiwu Zhang
2019/05/23
'''
import os
import numpy as np
from datareduction import IO
from astropy.io import fits



def Biasfitsread(path,biastxt,ext):
    '''
    this function is used to read all the bias images
    :param path: dictionary of all bias.fits
    :param biastxt: .txt file comtain all the name of bias.fits
    :return: image dictionaory, each cell represent an extension, it value is a 3D numpy array, each slice of this array is a bias image of this extension
    '''
    os.chdir(path)
    biaslist=IO.txt_io(path,biastxt)
    biasimgdic={x:[] for x in range(1,ext+1)}
    for biasname in biaslist:
        biasfits = fits.open(biasname)
        for extnum in range(1,ext+1):
            biasimg=biasfits[extnum].data
            biasimgdic[extnum].append(biasimg)
        biasfits.close()
    biasimgdic={x:np.array(biasimgdic[x]) for x in biasimgdic}
    return biasimgdic

Biasfitsread('/Users/shiwuzhang/work&study/ASTRO','bias.txt',4)