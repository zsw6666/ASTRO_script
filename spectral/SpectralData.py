import numpy as np
import matplotlib.pyplot as plt
import pyspeckit
from imag import ImgInterSmo
from cube import CubeData
from astropy import constants as const
from astropy import units as u
from scipy import stats

'''
This script used to access with the spectra
'''

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
    if wavecut is not None:
        lowlimit,uplimit=Findwavelength(wavelength,wavecut)
    # img=plt.plot(wavelength[800:],onedspectral[800:])
        img = plt.plot(wavelength[lowlimit:uplimit], onedspectral[lowlimit:uplimit])
    else:
        img=plt.plot(wavelength,onedspectral)
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

def Findwavelength(wavelength,waverange,threshold=1e-4):
    '''
    find the physical position of the selected wavelength in the wavelength interval
    :param wavelength: wavelength axis read from the IFU fits
    :param waverange: interval which you select
    :return: the physical position of the interval's boundaries
    '''
    lowerlimit, upperlimit = np.where(abs(wavelength - waverange[0]) <= threshold), np.where(
        abs(wavelength - waverange[1]) <= threshold)
    return lowerlimit[0][0],upperlimit[0][0]

def Specplotfit(filename,wavecut,fitmodel,guesses=None,rewave=None,position=None,axis=None,title='',
                color='black',xlabel='velocity km/s',ylabel='flux $erg/s/cm^{2}/ \AA$',fontsize=20.,annotate=''):
    '''
    new function for spectrum inspection for a specific position and aperture
    :param filename: data cube fits name
    :param wavecut: wavelength range
    :param fitmodel: model to fit the emission and absorption line
    :param guesses: initial parameter for the fit model
    :param rewave: referece wavelength for the convertion from wavelength to velocity
    :param position: position we want to extract wavelength from
    :param axis: figure.axes
    :param title: title of the plot
    :param color: color of line
    :param xlabel: x label
    :param ylabel: y label
    :param fontsize: size of words
    :param annotate:
    :return: None
    '''
    #read cube, select specific wavelength slices
    #mark the position, convert wavelength to velocity
    #normalize flux and subtract baseline
    cube=pyspeckit.Cube(filename)
    cube.cube=cube.cube*1e-16
    subcube=cube.slice(wavecut[0],wavecut[1],unit='Angstrom')
    spectra=subcube.get_apspec(position)
    spectra.xarr.refX=rewave*u.AA
    spectra.xarr.velocity_convention='optical'
    spectra.xarr.convert_to_unit(u.km/u.s)
    spectra.flux=spectra.flux/np.max(spectra.flux)
    spectra.baseline(order=0,subtract=False)


    # spectra.smooth(2)
    spectra.flux = spectra.flux / np.max(spectra.flux)

    #plot spectra
    if axis is None:
        spectra.plotter(plt.figure(1),xlabel='',ylabel='',title=title,color=color,linewidht=1.5)
    else:
        spectra.plotter(axis=axis, xlabel='', ylabel='', title=title,
                        color=color,clear=False,linewidth=1.5)


    #fit spectra with gaussian function usually
    spectra.specfit(fittype=fitmodel,lw=1.5)
    spectra.specfit.plot_components(add_baseline=True,component_yoffset=0.,lw=1.5)
    # spectra.specfit.plotresiduals(axis=axis,clear=False,yoffset=0,label=False,linewidth=1.5,color='gray')
    axis.set_xlabel('',fontsize=fontsize)
    axis.set_ylabel('',fontsize=fontsize)
    axis.annotate(annotate,(-2050,1.1), fontsize=fontsize)
    axis.tick_params(labelsize=15.)
    return None

def Slitspectrum(datacube,maskcube,x,y,wavelengthrange,z,intrinsicwavelength,horizontal=True):
    '''
    extrat spectra from a long slit
    :param datacube: IFU data
    :param maskcube: maskcube
    :param x: spatial axis of the long slit
    :param y: spatial axis of the long slit
    :param wavelengthrange: wavelength range of the spectra
    :param z:redshift of the target source
    :param intrinsicwavelength: intrinsic wavelength of the emission line
    :param horizontal: direction of the long slit
    if horizontal is True, then the long slit is horizontal
    else vertical
    :return: 2d long slit spectra
    '''
    if horizontal:
        twodspec=CubeData.OptimalextractImg(datacube,maskcube,1)
        spatial_axis=y
    else:
        twodspec=CubeData.OptimalextractImg(datacube,maskcube,2)
        spatial_axis=x
    velocityrange=wavelength2velocity(wavelengthrange,z,intrinsicwavelength)
    return twodspec[:,::-1],velocityrange,spatial_axis

def RegionSpectrum(datacube,maskcube=None):
    '''
    extract spectra from a small region
    :param datacube: wait for spectra extraction
    :param maskcube:
    :return: extracted spectra
    '''
    if maskcube is None:
        maskcube=np.ones(np.shape(datacube))
    spec1=np.sum(np.sum(datacube,axis=1),axis=1)
    mask1=np.sum(np.sum(maskcube,axis=1),axis=1)
    mask1[np.where(mask1==0)]=1
    return spec1/mask1

def Specnoiseestor(spectra):
    '''
    estimate noise level of an array
    :param spectra: input spectra
    :return: noise standard deviation
    '''
    #use highpass filter to remove the signal and only
    #noise left
    noisespec = ImgInterSmo.NoiseFilter(spectra, 5, 0.1, 'highpass')
    return np.std(noisespec),noisespec
