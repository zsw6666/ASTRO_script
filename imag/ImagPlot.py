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
    map = np.mean(cubedata[:, :, :], axis=0)


    return map

def ImgCut(image,cutlist):

    image_shape=np.shape(image)
    image[:, :int(cutlist[0] * image_shape[0])] = np.nan
    image[:, -1 - int(cutlist[1] * image_shape[0]):] = np.nan
    image[:int(cutlist[2] * image_shape[0]), :] = np.nan
    image[-1 - int(cutlist[3]* image_shape[0]):, :] = np.nan

    return image

def Twodplotimg(map,x=None,y=None,contourmap=None,subrow=1,subclo=1,vmin=None,vmax=None,contourlevel=None,
                xlabel=r'x',ylabel=r'y',cbarlabel='cbar',subtitle=None,title=None):
    fig,AX=plt.subplots(subrow,subclo,num=len(map),sharey=True,sharex=True)
    fig.text(0.5, 0.3, xlabel, ha='center',fontsize=25.)
    fig.text(0.08, 0.6, ylabel, va='center', rotation='vertical',fontsize=25.)
    fig.suptitle(title,fontsize=20.)

    l_map=0
    if (subclo is 1) and (subrow is 1):
        AX=[AX]
        map=[map]
    else:
        AX = AX.flatten()
    for i in range(len(AX)):
        AX[i].set(adjustable='box-forced', aspect='auto')
        if l_map<=len(map)-1:
            img=AX[i].pcolor(x,y,map[l_map],cmap='gist_ncar',vmax=vmax,vmin=vmin)
            AX[i].set_title(subtitle[l_map],fontsize=20.)
            AX[i].tick_params(labelsize=15.)
            cbar = fig.colorbar(img, ax=AX[i],orientation="horizontal", pad=0.2)#
            cbar.set_label(cbarlabel[l_map],fontsize=25.)
            cbar.ax.tick_params(labelsize=15.)
        else:
            fig.delaxes(AX[i])
        if contourmap is not None:
            con=AX[i].contour(x,y,contourmap,levels=contourlevel,colors='black')
            AX[i].clabel(con, inline=1, fontsize=15.)

        l_map += 1
    plt.subplots_adjust(wspace=0, hspace=0)
    return fig,AX

def Mapprepro(twodmap,x,y,internum):
    twodmap = ImgInterSmo.InSmimg(twodmap, internum, [3., 3.])
    x, y = ImgInterSmo.Arrayinterpolation(x, internum)[1:-1], ImgInterSmo.Arrayinterpolation(y, internum)[1:-1]
    return twodmap,x,y

def Threedplotimg(twodmap,x,y):
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = Axes3D(fig)
    img=ax.plot_surface(X, Y, twodmap, rstride=1, cstride=1, cmap='jet')
    return img

def Cutmap(datacube,cutaxis,waverange,xaxis,yaxis,threshold=1e-4):
    cutflux, cutwavelength = CubeData.CubeCut(datacube, cutaxis, 'manual', waverange,threshold)
    twodimg = PlotMap(cutflux)
    return twodimg,xaxis,yaxis

def Scaleimgconverter(img):
    norm=ImageNormalize(img,interval=ZScaleInterval(),stretch=LinearStretch())
    return norm

def Contourgenerator(handle=None,contourimg=None,wcs=None,levels=None,color=None,labels=None):
    contour=handle.contour(contourimg,transform=handle.get_transform(wcs),levels=levels,linewidths=1.5,colors=color)

    return contour,handle


def Gimgplot(fig=None,AX=None,imglist=None,x=None,y=None,
             contourmap=None,contourlevel=None,xlabel=None,
             ylabel=None,cbarlabel=None,imgname=None,contourmark=None):
    fig.text(0.48, 0.27, xlabel, ha='center', fontsize=25.)
    fig.text(0.05, 0.6, ylabel, va='center', rotation='vertical', fontsize=25.)
    for i in range(len(AX)):
        img = AX[i].pcolor(x, y, imglist[i], cmap='gist_ncar')
        if contourmark[i] is not None:
            con=AX[i].contour(x,y,contourmap,levels=contourlevel,colors='black')
            AX[i].clabel(con, inline=1, fontsize=15.)
        if imgname is not None:
            AX[i].text(-18, 7, imgname[i], fontsize=25.,color='red')
        cbar = fig.colorbar(img, ax=AX[i], orientation="horizontal",aspect=20)  #
        cbar.set_label(cbarlabel[i], fontsize=25.)
        cbar.ax.tick_params(labelsize=15.)
        AX[i].tick_params(labelsize=15.)
        # AX[i].set_aspect(1)