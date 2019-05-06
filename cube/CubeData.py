import os
import numpy as np
from datareduction import WCS
from spectral_cube import SpectralCube
from astropy import constants as const
from astropy import units as u

__all__=['ReadCube','CubeCut','FLux','ContinuumSub','ContinuumEst','BackgroundEstimation','ImageBackgroundSubtraction','CubeBackgroundSubtraction','WCS']

def ReadCube(path, cube_name):
    '''
    load the fits and print the basic information of this cube
    :param path: location of the fits
    :param cube_name: name of the fits
    :return: 3D data and wavelength axis
    '''

    os.chdir(path)
    cube = SpectralCube.read(cube_name)
    wcs = cube.world[:]
    print(cube)
    # ra1, dec1 = (wcs[:, :, 0].T)[::-1, :], (wcs[:, :, 1].T)[:, ::-1]
    # ra=cube.world[0,:,:][2].value
    # dec=cube.world[0,:,:][1].value
    return cube._data, cube.spectral_axis, wcs, cube


def CubeCut(cube=None, cube_wavelength=None, emissionline=None, rang=None):
    '''
    cut the cube and return the part we need
    :param cube: cube data
    :param cube_wavelength: wavelength
    :param emissionline: select different part according to this parameter
    :param rang: if set emission line to manual,then this parameter must be given
    :return: part of the initial cube and corresponding wavelength
    '''
    if emissionline == 'lyman':
        # cube_cut,wavelength_cut=cube[950:1010,2:67,:],cube_wavelength[950:1010]
        cube_cut, wavelength_cut = cube[940:990,
                                        2:67, :], cube_wavelength[940:990]
    elif emissionline == 'CV':
        # print(cube_wavelength[1030:1200])
        # cube_cut,wavelength_cut=cube[1030:1200,2:67,:],cube_wavelength[1030:1200]
        cube_cut, wavelength_cut = cube[1030:1500,
                                        2:67, :], cube_wavelength[1030:1500]
    elif emissionline == 'manual':
        cube_cut, wavelength_cut = cube[rang[0]:rang[1],2:67,:],cube_wavelength[rang[0]:rang[1]]
    else:
        cube_cut, wavelength_cut = cube[:, 2:67, :], cube_wavelength

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
    ra=np.mean(ramap[2:67,:],axis=1)
    decmap=np.mean(wcsmap[1],axis=0)
    dec=np.mean(decmap,axis=0)

    return ra,dec