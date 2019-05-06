import numpy as np
from astropy.io import fits
from . import IO


def img_sub(img1,img2):

    img_sub=img1-img2#img1 and img2 must be np.array and have the same size
    return img_sub

