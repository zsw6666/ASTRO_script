import os
import numpy as np
from spectral import SpectralData
from spectral_cube import SpectralCube
from imag import ImagPlot,ImgInterSmo
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import Angle
from astropy.io import fits

'''
This script is used to access the datacube
'''

def AccessCube(path,cube_name):
    os.chdir(path)
    cube = SpectralCube.read(cube_name)
    return cube

def ReadCube(path, cube_name):
    '''
    load the fits and print the basic information of this cube
    :param path: location of the fits
    :param cube_name: name of the fits
    :return: 3D data and wavelength axis
    '''

    cube=AccessCube(path,cube_name)
    wcs = cube.world[:]
    print(cube)
    return cube._data, cube.spectral_axis, wcs, cube

def Editheader(path,name,keyword,newvalue):
    os.chdir(path)
    fitsfile=fits.open(name,mode='update')
    fitsfile[0].header[keyword]=newvalue
    fitsfile.close()
    return None

def Readheader(path,name):
    os.chdir(path)
    fitsfile = fits.open(name)
    header=fitsfile[0].header
    fitsfile.close()
    return header


def CubeCut(cube=None, cube_wavelength=None, emissionline=None,waverange=None):
    '''
    cut the cube and return the part we need
    :param cube: cube data
    :param cube_wavelength: wavelength
    :param emissionline: select different part according to this parameter
    :param rang: if set emission line to manual,then this parameter must be given
    :return: part of the initial cube and corresponding wavelength
    '''
    cutrange = SpectralData.Findwavelength(cube_wavelength, waverange)
    if emissionline == 'manual':
        cube_cut, wavelength_cut = cube[cutrange[0]:cutrange[1],:,:],cube_wavelength[cutrange[0]:cutrange[1]]
    else:
        # cube_cut, wavelength_cut = cube[:, 2:67, :], cube_wavelength
        cube_cut, wavelength_cut = cube[:, :, :], cube_wavelength

    return cube_cut, wavelength_cut

def FLux(cube=None, wavelength_axis=None,gain=None, exposuretime=None,emissionline=None):
    """
            convert number of electrons per second of pixel to flux
            :param cube: cube data
            :param emissionline: this parameter determine which part of the given cube we need
            :param wavelength_axis: wavelength used to calculate flux
            :param gain: ccd gain number of electrons / number of photons
            :return: flux and corresponding wavelength
    """
    # cube_all = CubeCut(cube=cube, cube_wavelength=wavelength_axis,
    #                    emissionline=emissionline)

    cubedata, wavelength_cut = cube[:,2:67,:],wavelength_axis
    photons = cubedata/gain
    flux=[]
    for i in range(len(wavelength_cut)):
        flux.append(((photons[i]*const.h*(const.c/wavelength_cut[i]))/(exposuretime*u.s)).to(u.erg/u.s))
    flux=np.array(flux)*(u.erg/u.s)

    return flux, wavelength_cut
def ContinuumEst(flux, ran=[1300, 1500]):
    '''
    estimate the continuum component of flux
    :param flux: total flux
    :param ran: wavelength range used to calculate continuum
    :return: continuum component
    '''

    continuum = np.mean(flux[ran[0]:ran[1], :, :], axis=0)
    continuum_std = np.std(flux[ran[0]:ran[1], :, :], axis=0)

    return continuum, continuum_std
def ContinuumSub(flux, flux_all):
    '''
    subtract continuum component
    :param flux: total flux
    :param flux_all: used to estimate continuum
    :return: flux and its standard deviation
    '''

    flux_sub = np.zeros(np.shape(flux))
    continuum, continuum_std = ContinuumEst(flux_all)
    for i in range(np.shape(flux)[0]):
        flux_sub[i, :, :] = flux[i, :, :]-continuum
    flux_std = continuum_std
    return flux, flux_std
def BackgroundEstimation(image, region):
    '''
    estimate background value of the given image
    :param image: wait to be analysis
    :param region: which part we use to estimate background
    :return: background and its standard deviation
    '''

    background = np.mean(
        image[region[0, 0]:region[0, 1], region[1, 0]:region[1, 1]])
    std_background = np.std(
        image[region[0, 0]:region[0, 1], region[1, 0]:region[1, 1]])

    return background, std_background
def ImageBackgroundSubtraction(image, background):
    '''
    subtract background from the given image
    :param image: wait to be subtracting
    :param background: background value
    :return: background-subrtacted image
    '''

    image = image-background
    return image
def CubeBackgroundSubtraction(cube):
    '''
    subtract background from the given cube
    :param cube: cube waited to be subtract
    :return: background-subtracted cube
    '''
    cube_bk = np.zeros(np.shape(cube))
    for i in range(len(cube)):
        background, background_std = BackgroundEstimation(
            cube[i], np.array([[11, 20], [5, 9]]))
        cube_bk[i] = ImageBackgroundSubtraction(cube[i], background)

    return cube_bk

def WCS(wcsmap):
    '''
    calculate ra and dec of this cube
    :param wcsmap: 3D cube for wavelength ra and dec
    :return: ra and dec all 1D array
    '''
    ramap=np.mean(wcsmap[2],axis=0)
    ra=np.mean(ramap[:,:],axis=1)
    decmap=np.mean(wcsmap[1],axis=0)
    dec=np.mean(decmap,axis=0)

    return ra,dec

def Findposition(ramap,decmap,postion):
    '''
    find the physical coordinate corresponding to ra and dec
    :param ramap: 2D array of ra
    :param decmap: 2D array of dec
    :param postion: position selected
    :return: physical coordinate of the position selected
    '''

    ra_physical=np.where(abs(postion[0]-ramap)<=5.5e-5)
    dec_physical=np.where(abs(postion[1]-decmap)<=2e-4)
    return ra_physical[0][0],dec_physical[0][0]

def MarkSource(twodmap,ramap,decmap,position):
    '''
    mark the position selected
    :param twodmap: image
    :param ramap: 2D array of ra
    :param decmap: 2D array of dec
    :param position: postion selected
    :return: marked 2D image
    '''

    ra_physical,dec_physical=Findposition(ramap,decmap,position)
    twodmap[ra_physical-round(5/2):ra_physical+round(5/2),dec_physical-round(3/2):dec_physical+round(3/2)]=0.5*np.min(twodmap)
    return twodmap

def Findabsorption(fluxcube,cubewavelength,abrange,upperrange,lowerrange,xaxis,yaxis):
    medianimg,axis1,axis2=ImagPlot.Cutmap(fluxcube,cubewavelength,abrange,xaxis,yaxis)
    upperimg,axis1,axis2=ImagPlot.Cutmap(fluxcube,cubewavelength,upperrange,xaxis,yaxis)
    lowerimg,xaxis,yaxis=ImagPlot.Cutmap(fluxcube,cubewavelength,lowerrange,xaxis,yaxis)
    umimg=upperimg-medianimg
    lmimg=lowerimg-medianimg
    img=np.array([umimg,lmimg])
    img=np.median(img,axis=0)
    return img,umimg,lmimg,xaxis,yaxis

def Coordinateconvert(position):
    ra=Angle(str(position[0])+'d')
    dec=Angle(str(position[1])+'d')

    return ra.hms,dec.dms

def Imgreadout(img,header,name):
    '''
    convert the image to fits file
    :param img: image of the fits file
    :param header: header fo the fits file
    :param name: name of the fits file
    :return: None
    '''

    hdulist=fits.HDUList()
    imghdu=fits.ImageHDU(data=img,header=header)
    hdulist.append(imghdu)
    hdulist.writeto(name)
    return None

def Cubeweightedmean(cube,weight):
    '''
    use to calculate the flux-weighted velocity map
    :param cube: datacube
    :param weight: velocity correspond to the wavelength
    :return: flux-weighted velocity map
    '''

    #for each wavelength, smooth the image firstly
    cube=ImgInterSmo.CubeSmooth(cube,[1.5,0.428])

    #sum the intensity for each pixel as the denominator
    #find the pxiel with negative flux and repalce them with a very high value,
    #for these kind of pixel, there's no emiision line, with this step when calculate
    #the velocity, it will be a very small value(very close to zero) which will not influence
    #the velocity map
    totalmap = np.sum(cube, axis=0)
    totalmap[np.where(totalmap<=0.)]=1e8

    #for image of each wavelength, multiply it with the weight(image is the flux array and weight is the velocity corresponding
    #to the wavelength )
    for i in range(len(cube)):
        cube[i,:,:]=cube[i,:,:]*weight[i]

    #divide the multipiled cube by the totalmap to generate the flux-weighted velocity map
    meanmap=np.sum(cube,axis=0)/totalmap
    # meanmap[np.where(meanmap==np.min(meanmap))]=0.
    return meanmap

def Angle2distance(ra,dec):
    '''
    convert ra,dec to angle distance to the center of the image
    :param ra: ra
    :param dec: dec
    :return: angle distance from the center of the image
    '''

    ra,dec=ra.to(u.rad),dec.to(u.rad)
    ra_dis,dec_dis=(ra-np.median(ra)).to(u.arcsec),(dec-np.mean(dec)).to(u.arcsec)

    return ra_dis,dec_dis

def CubeNoiseFilter(cube,N=5,wn=.3):
    '''
    reduce the noise for the spectra of each pixel
    :param cube: datacube
    :param N: control the filter
    :param wn: also used to control the filter
    :return: noise-filtered datacube
    '''

    shape=np.shape(cube)

    #for each pixel reduce the noise
    for i in range(shape[1]):
        for j in range(shape[2]):
            cube[:,i,j]=ImgInterSmo.NoiseFilter(cube[:,i,j],N,wn)
    return cube

def Cubeseeinglimit(cube,size=[3.,1.]):
    '''
    because the length and width of each pixel stand for different angle scale, for ra the angle scale
    for each pixel is 0.38 arcsec, so use three pixels along ra axis which corresponds to 1.14 arcsec
    for dec, the angle scale for single pxiel is 1.35 arcsec. So use 3*1 Rectangle as a cell.
    :param cube: datacube
    :param size: size of the cell
    :return: datacube after considering the seeing
    '''

    cube_shape=np.shape(cube)

    #for each wavelength
    for i in range(cube_shape[0]):
        cube[i,:,:]=ImgInterSmo.Imgseeinglimit(cube[i,:,:],size)
    return cube