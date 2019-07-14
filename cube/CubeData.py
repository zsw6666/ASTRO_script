import os
import numpy as np
import matplotlib.pyplot as plt
from spectral import SpectralData
from spectral_cube import SpectralCube
from imag import ImgInterSmo
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

def WCSextractor(path,fielname):
    _, _, wcs, _=ReadCube(path,fielname)
    ra,dec=WCS(wcs)
    ra_dis,dec_dis=Angle2distance(ra, dec)
    return ra_dis.value,dec_dis.value

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

def Continumsubtractor(cubecontinum,cubesub):

    img_median=np.median(cubecontinum,axis=0)
    cubesub=cubesub-img_median
    return cubesub

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

def Cubeweightedmean(cube,weight,mask_cube):
    '''
    use to calculate the flux-weighted velocity map
    :param cube: datacube
    :param weight: velocity correspond to the wavelength
    :return: flux-weighted velocity map
    '''


    raw_cube=cube.copy()

    # we use optimal extraction method to generate the flux-weighted map,
    # to generate the optimal flux map, we first multiply the cube with the mask
    # sum the cube up along the axis of wavelength and get the flux map,
    # then by dividing the flux map with the sumed-up mask we ge the intensity map
    cube=cube*mask_cube
    totalmap = np.sum(cube, axis=0)
    mask=np.sum(mask_cube,axis=0)
    mask[np.where(mask==0.)]=1.
    avermap=totalmap/mask




    #for image of each wavelength, multiply it with the weight(image is the flux array and weight is the velocity corresponding
    #to the wavelength )
    cube_vel=cube.copy()
    # cube_vel=CubeNoiseFilter(cube_vel,3,.2)
    for i in range(len(cube)):
        cube_vel[i,:,:]=cube_vel[i,:,:]*weight[i]

    #divide the multipiled cube by the totalmap to generate the flux-weighted velocity map
    meanmap=np.sum(cube_vel,axis=0)/totalmap
    # meanmap[np.where(meanmap==np.min(meanmap))]=0.
    return meanmap,avermap

def Angle2distance(ra,dec,refpoint):
    '''
    convert ra,dec to angle distance to the center of the image
    :param ra: ra
    :param dec: dec
    :return: angle distance from the center of the image
    '''

    ra,dec=ra.to(u.rad),dec.to(u.rad)
    ra_dis,dec_dis=(ra-refpoint[0]*u.deg).to(u.arcsec),(dec-refpoint[1]*u.deg).to(u.arcsec)

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

def Cubeseeinglimit(cube,size=[6.,4.]):
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

def Cubebadpixelremovor(cube,mean=0.,sigma=1.):
    '''
    some pixels' value is negative, some should be zero, this function remove these kinds of
    pixels
    :param cube: data cube
    :param mean: average value of this cube
    :param sigma: standard deviation
    :return: bad-pixel-removed cube
    '''

    cube[np.where(cube<(mean+sigma))]=0.
    return cube

def Regionfilter(img1,img2,sigma_num=0.3):
    '''
    select the region in image 2 whose value is beyond the threshold in image 1
    :param img1: reference image
    :param img2: image wait for filter
    :param sigma_num: threshold
    :return: filtered image
    '''
    img_return=img2.copy()
    img_return[np.where(img1<sigma_num)]=np.nan
    return img_return

def BadvalueEstimator(cube):
    '''
    estimate the threshold which is used to remove the bad pixel in image
    :param cube: used to estimate the threshold
    :return: threshold
    '''


    img_sigma=np.std(cube)
    return img_sigma

def Maskgenerator(data_cube,threshold):
    '''
    generathe the mask cube for the data cube, it check each slice in cube
    and select pixels whose value is beyond the threshold of its own slice
    :param data_cube: wait for mask
    :param threshold: threshold
    :return: mask cube for the data cube
    '''

    #generate the mask cube
    mask_cube=np.zeros(np.shape(data_cube))

    # for each slice, generate a mask with the given threshold
    for data, mask in zip(data_cube,mask_cube):
        index_set=np.where(data>=threshold)
        mask[index_set]=1.

    return mask_cube

def Maskgenerator2(data_cube,n_sigma):
    '''
    generate the mask cube, it check the spectra of each point, only select
    pixels whose value is beyond the threshold of this spectra
    :param data_cube: wait for mask
    :param n_sigma: threshold
    :return: mask cube
    '''

    #generate the mask cube
    mask_shape = np.shape(data_cube)
    mask_cube = np.zeros(mask_shape)

    #check the spectra of each point and select the pixel satisfied the condition
    for i in range(mask_shape[1]):
        for j in range(mask_shape[2]):
            spec_aver = np.mean(data_cube[:, i, j])
            spec_std = np.std(data_cube[:, i, j])
            index_set = np.where(data_cube[:, i, j] >= spec_aver + n_sigma * spec_std)
            mask_cube[index_set, i, j] += 1

    return mask_cube