import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from pylab import pcolor,colorbar
from mpl_toolkits.mplot3d import Axes3D

__all__=['PlotMap','Colormap']

def PlotMap(cubedata,emissionline=None):
    '''
    plot the 2D image of the cube by summing along the wavelength axis to get the total flux for each pixel
    :param cubedata: 3D cube data
    :return: 2D image array
    '''

    if emissionline=='lyman':
        map=np.sum(cubedata[950:1010,:,:],axis=0)
    elif emissionline=='CV':
        map=np.sum(cubedata[1050:1330,:,:],axis=0)
    else:
        map = np.sum(cubedata[:, :, :], axis=0)


    return map

def Colormap(map,x,y):
    plt.figure()
    img=pcolor(x,y,map)
    colorbar()
    return img

def ThreeDmap(twodmap,x,y):
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, twodmap, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    plt.show()