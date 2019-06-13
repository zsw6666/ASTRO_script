import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.visualization import ZScaleInterval,LinearStretch,ImageNormalize
from pylab import pcolor,colorbar
from mpl_toolkits.mplot3d import Axes3D
from cube import CubeData
from . import ImgInterSmo


'''
This script contains functions used to plot image
'''


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
        map = np.median(cubedata[:, :, :], axis=0)


    return map

def Twodplotimg(map,x,y,vmin=None,vmax=None):
    plt.figure()
    ax=pcolor(y,x,map.T,cmap=cm.jet,vmax=vmax,vmin=vmin)
    cbar=colorbar()
    return ax,cbar

def Mapprepro(twodmap,x,y,internum):
    twodmap = ImgInterSmo.InSmimg(twodmap, internum, [3., 3.])
    x, y = ImgInterSmo.Onedinterpolation(x, internum)[1:-1], ImgInterSmo.Onedinterpolation(y, internum)[1:-1]
    return twodmap,x,y

def Threedplotimg(twodmap,x,y):
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = Axes3D(fig)
    img=ax.plot_surface(X, Y, twodmap, rstride=1, cstride=1, cmap=cm.jet)
    return img

def Cutmap(datacube,cutaxis,waverange,xaxis,yaxis):
    cutflux, cutwavelength = CubeData.CubeCut(datacube, cutaxis, 'manual', waverange)
    twodimg = PlotMap(cutflux)
    # twodimg ,xaxis,yaxis=Mapprepro(twodimg,xaxis,yaxis,internum=0)
    return twodimg,xaxis,yaxis

def Scaleimgconverter(img):
    norm=ImageNormalize(img,interval=ZScaleInterval(),stretch=LinearStretch())
    return norm

def Contourgenerator(handle=None,contourimg=None,wcs=None,levels=None,color=None):
    contour=handle.contour(contourimg,transform=handle.get_transform(wcs),levels=levels,linewidths=.5,colors=color)

    return handle,contour