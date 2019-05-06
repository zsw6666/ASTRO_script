import numpy as np
import matplotlib.pyplot as plt
from cube import CubeData

__all__=['Spectral','SourceSpectral','SpectralPlot','TwoDSpectral','Mapspectral']

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


def SourceSpectral(cube, position,size):
    '''
    output the 1D spectral of selected sources
    :param cube: cube data which is a 3D numpy array
    :param position: where you want the spectral(physical position)
    :return: spectrals
    '''
    selected_region=cube[:,position[0]-round(size[0]/2.):position[0]+round(size[0]/2.),position[1]-round(size[1]/2.):position[1]+round(size[1]/2.)]
    onedspectral=np.sum(np.sum(selected_region,axis=1),axis=1)

    return onedspectral

def Mapspectral(fluxcube,size):
    width=int(np.shape(fluxcube)[1]/size[0])
    length=int(np.shape(fluxcube)[2]/size[1])
    spectral_list=[]
    position_list=[]
    k=0
    for i in range(1,width):
        for j in range(1,length):
            k+=1
            position=[i*size[0],j*size[1]]
            ondspectral=SourceSpectral(fluxcube,position,size)
            print(k)
            print(position)
            spectral_list.append(ondspectral)
            position_list.append(position)
    return spectral_list,position_list
def SpectralPlot(onedspectral, wavelength,title=None):
    '''
    accept the spectral from SourceSpectral and plot them
    :param spectral_set: spectral of each sources
    :param wavelength: wavelength used to plot spectral
    :param source_coordinate: mark the the spectral to tell to which source it belongs
    :return: None
    '''
    wavelength = wavelength.value
    plt.title(title)
    plt.xlabel(r'$wavelength \ A$')
    plt.ylabel(r'$flux \ erg/s/m^{2}/ \AA$')
    img=plt.plot(wavelength,onedspectral)
    return img

def TwoDSpectral(cube, rang=None):
    '''
    extract the 2D spectra from the cube
    :param cube:
    :return: 2D spectral array
    '''
    cube = CubeData.CubeCut(cube, emissionline='manual', rang=rang)
    twod_spectral = np.sum(cube, axis=2)
    return twod_spectral