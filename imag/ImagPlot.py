import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval,LinearStretch,ImageNormalize
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
    map = np.sum(cubedata[:, :, :], axis=0)


    return map

def Twodplotimg(map,x,y,subrow=1,subclo=1,vmin=None,vmax=None,markpoint='mark',xlabel=r'x',ylabel=r'y',cbarlabel='cbar',subtitle=None):
    fig,AX=plt.subplots(subrow,subclo,num=len(map))
    fig.text(0.5, 0.04, xlabel, ha='center')
    fig.text(0.04, 0.5, ylabel, va='center', rotation='vertical')
    l_map=0
    if (subclo is 1) and (subrow is 1):
        AX=[AX]
        map=[map]
    else:
        AX = AX.flatten()
    for i in range(len(AX)):
        AX[i].set(adjustable='box-forced', aspect='auto')
        if l_map<=len(map)-1:
            img=AX[i].pcolor(x,y,map[l_map],cmap='jet',vmax=vmax,vmin=vmin)
            AX[i].set_title(subtitle[l_map])
            cbar = fig.colorbar(img, ax=AX[i])
            cbar.set_label(cbarlabel)
        else:
            fig.delaxes(AX[i])
        if markpoint is not None:
            AX[i].scatter(0.,0.,marker='*',color='magenta',s=100)

        l_map += 1
    # cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    # cbar=fig.colorbar(img,cax)
    # cbar.set_label(cbarlabel)
    fig.tight_layout(rect=[0,0,.92,1])
    plt.show()


    return None

def Mapprepro(twodmap,x,y,internum):
    twodmap = ImgInterSmo.InSmimg(twodmap, internum, [3., 3.])
    x, y = ImgInterSmo.Onedinterpolation(x, internum)[1:-1], ImgInterSmo.Onedinterpolation(y, internum)[1:-1]
    return twodmap,x,y

def Threedplotimg(twodmap,x,y):
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = Axes3D(fig)
    img=ax.plot_surface(X, Y, twodmap, rstride=1, cstride=1, cmap='jet')
    return img

def Cutmap(datacube,cutaxis,waverange,xaxis,yaxis):
    cutflux, cutwavelength = CubeData.CubeCut(datacube, cutaxis, 'manual', waverange)
    twodimg = PlotMap(cutflux)
    return twodimg,xaxis,yaxis

def Scaleimgconverter(img):
    norm=ImageNormalize(img,interval=ZScaleInterval(),stretch=LinearStretch())
    return norm

def Contourgenerator(handle=None,contourimg=None,wcs=None,levels=None,color=None,labels=None):
    contour=handle.contour(contourimg,transform=handle.get_transform(wcs),levels=levels,linewidths=.5,colors=color)

    return contour,handle