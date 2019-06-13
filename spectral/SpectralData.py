import numpy as np
import matplotlib.pyplot as plt
from cube import CubeData
from astropy import constants as const

'''
This script used to access with the spectra
'''


def Spectral(cube, source_coordinate, wavelength):
    '''
    :param cube: cube data
    :param source_coordinate: coordinate of the sources
    :param wavelength: wavelength
    :return: None
    '''
    spectral_set = SourceSpectral(cube, source_coordinate)
    SpectralPlot(spectral_set, wavelength, source_coordinate)
    return None

def SourceSpectral(cube,ssize,position,ramap=None,decmap=None):
    '''
    output the 1D spectral of selected sources
    :param cube: cube data which is a 3D numpy array
    :param position: where you want the spectral(physical position)
    :return: spectrals
    '''
    if ramap is not None and decmap is not None:
        position=CubeData.Findposition(ramap,decmap,position)

    #for each wavelength using the median of these pxiels' value as the intensity at this wavelength
    selected_region=cube[:,position[0]-round(ssize[0]/2.):position[0]+round(ssize[0]/2.),position[1]-round(ssize[1]/2.):position[1]+round(ssize[1]/2.)]
    onedspectral=np.median(np.median(selected_region,axis=1),axis=1)

    return onedspectral

def Mapspectral(fluxcube,ssize):
    '''
    generate t spectras for the whole datacube
    :param fluxcube: datacube used to produce the spectra
    :param ssize: size whin which to calculate the spectra
    :return: spectra_list contains the spectra for each position
    '''

    width=int(np.shape(fluxcube)[1]/ssize[0])
    length=int(np.shape(fluxcube)[2]/ssize[1])
    spectral_list=[]
    position_list=[]

    #for each position, generate a spectra and put it into the spectra)list
    for i in range(1,width):
        for j in range(1,length):
            position=[i*ssize[0],j*ssize[1]]
            ondspectral=SourceSpectral(fluxcube,position=position,ssize=ssize)
            spectral_list.append(ondspectral)
            position_list.append(position)
    return spectral_list,position_list

def SpectralPlot(onedspectral, wavelength,wavecut=None,title=None):
    '''
    accept the spectral from SourceSpectral and plot them
    :param spectral_set: spectral of each sources
    :param wavelength: wavelength used to plot spectral
    :param source_coordinate: mark the the spectral to tell to which source it belongs
    :return: None
    '''
    plt.title(title)
    plt.xlabel(r'$wavelength \ (\AA)$')
    plt.ylabel(r'$flux \ (erg/s/m^{2}/ \AA)$')
    lowlimit,uplimit=Findwavelength(wavelength,wavecut)
    # img=plt.plot(wavelength[800:],onedspectral[800:])
    img = plt.plot(wavelength[lowlimit:uplimit], onedspectral[lowlimit:uplimit])
    return img

def wavelength2velocity(wavelength,z,intrinsicwavelength):
    '''
    convert wavelength to velocity
    :param wavelength: wavelength interval
    :param z: redshift use to for this calculation
    :param intrinsicwavelength: the intrinsic wavelength for lyman-lapha it's 1216 \AA
    :return: velocity array
    '''

    #calculate the observed wavelength at that redshift
    lamda0 = intrinsicwavelength * (z + 1)
    velocity = const.c.value * (wavelength- lamda0) / lamda0
    return velocity

def WavelengthSelectPlot(spectral,waverange):
    '''
    plot the boundary of which range you select
    :param spectral: a spectral
    :param waverange: interval of wavelength which you select
    :return: two plot boundaries on image
    '''
    lineupper=np.max(spectral)
    lowerwave=np.repeat(waverange[0],100)
    upperwave=np.repeat(waverange[1],100)
    lowerline=plt.plot(lowerwave,np.linspace(np.min(spectral),lineupper,100),c='r')
    upperline=plt.plot(upperwave,np.linspace(np.min(spectral),lineupper,100),c='r')
    return lowerline,upperline

def Findwavelength(wavelength,waverange):
    '''
    find the physical position of the selected wavelength in the wavelength interval
    :param wavelength: wavelength axis read from the IFU fits
    :param waverange: interval which you select
    :return: the physical position of the interval's boundaries
    '''
    lowerlimit, upperlimit = np.where(abs(wavelength - waverange[0]) <= 1e-4), np.where(
        abs(wavelength - waverange[1]) <= 1e-4)
    return lowerlimit[0][0],upperlimit[0][0]