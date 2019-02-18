import numpy as np
import matplotlib.pyplot as plt
import CubeData


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


def SourceSpectral(cube, source_coordinate):
    '''
    output the 1D spectral of selected sources
    :param cube: cube data which is a 3D numpy array
    :param source_coordinate: coordinate of selected sources which is a 2D numpy array
    :return: spectrals
    '''
    spectral_set = []
    for individual in source_coordinate:
        # add neighboor pixels' value up to imporve SN
        cut = cube[:, individual[0]-1:individual[0] +
                   2, individual[1]-1:individual[1]+2]
        spectral_individual = np.sum(np.sum(cut, axis=1), axis=1)
        spectral_set.append(spectral_individual)
    spectral_set = np.array(spectral_set)

    return spectral_set


def SpectralPlot(spectral_set, wavelength, source_coordinate):
    '''
    accept the spectral from SourceSpectral and plot them
    :param spectral_set: spectral of each sources
    :param wavelength: wavelength used to plot spectral
    :param source_coordinate: mark the the spectral to tell to which source it belongs
    :return: None
    '''
    n_source = len(spectral_set)
    wavelength = wavelength.value
    for i in range(1, n_source+1):
        plt.figure(i)
        plt.plot(wavelength, spectral_set[i-1])
        plt.text(x=0.99*np.median(wavelength), y=np.mean(spectral_set[i-1]), s=str(
            source_coordinate[i-1, 0])+','+str(source_coordinate[i-1, 1]))
        plt.xlabel(r'wavelength')
        plt.ylabel(r'flux')
        plt.title(r'spectral')
    plt.show()
    return None

def TwoDSpectral(cube, rang=None):
    '''
    extract the 2D spectra from the cube
    :param cube:
    :return: 2D spectral array
    '''
    cube = CubeData.CubeCut(cube, emissionline='manual', rang=rang)
    twod_spectral = np.sum(cube, axis=2)
    return twod_spectral