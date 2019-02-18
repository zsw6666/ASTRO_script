import os
import numpy as np
import WCS
from spectral_cube import SpectralCube
from astropy import constants as const
from astropy import units as u


def ReadCube(path, cube_name):
    '''
    load the fits and print the basic information of this cube
    :param path: location of the fits
    :param cube_name: name of the fits
    :return: 3D data and wavelength axis
    '''

    os.chdir(path)
    cube = SpectralCube.read(cube_name)
    wcs = WCS.WCS(cube.header, cube.shape)
    print(cube)
    ra1, dec1 = (wcs[:, :, 0].T)[::-1, :], (wcs[:, :, 1].T)[:, ::-1]
    # ra=cube.world[0,:,:][2].value
    # dec=cube.world[0,:,:][1].value
    return cube._data, cube.spectral_axis, dec1, ra1, cube


def CubeCut(cube=None, cube_wavelength=None, ra=None, dec=None, emissionline=None, rang=None):
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
        cube_cut, wavelength_cut = cube[rang[0, 0]:,
                                        rang[1, 0]:rang[1, 1], rang[2, 0]:rang[2, 1]]
    else:
        cube_cut, wavelength_cut = cube[:, 2:67, :], cube_wavelength
    if ra is not None and dec is not None:
        ra_cut = ra[2:67, :]
        dec_cut = dec[2:67, :]

    else:
        ra_cut = None
        dec_cut = None
    return cube_cut, wavelength_cut, ra_cut, dec_cut

def FLux(cube=None, wavelength_axis=None, ra=None, dec=None, gain=None, emissionline=None):
    """
            convert number of electrons per second of pixel to flux
            :param cube: cube data
            :param emissionline: this parameter determine which part of the given cube we need
            :param wavelength_axis: wavelength used to calculate flux
            :param gain: ccd gain number of electrons / number of photons
            :return: flux and corresponding wavelength
    """
    cube_all = CubeCut(cube=cube, cube_wavelength=wavelength_axis,
                       ra=ra, dec=dec, emissionline=emissionline)
    cubedata, wavelength_cut = cube_all[0], cube_all[1]
    photons = cubedata/gain
    wavelength_axis = np.array((wavelength_cut).to(u.meter).data)
    for i in range(np.shape(wavelength_axis)[0]):
        photons[i] = photons[i]*const.h*(const.c.data/wavelength_axis[i])

    flux = photons

    return flux, wavelength_cut, cube_all[2], cube_all[3]

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