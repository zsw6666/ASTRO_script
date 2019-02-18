import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr


def PlotMap(cubedata,cube_wavelength,waveprint=True,emissionline=None):
    '''
    plot the 2D image of the cube by summing along the wavelength axis to get the total flux for each pixel
    :param cubedata: 3D cube data
    :return: 2D image array
    '''

    if emissionline=='lyman':
        map=np.sum(cubedata[950:1010,2:67,:],axis=0)
        wavelengthrange=cube_wavelength[950:1010]
    elif emissionline=='CV':
        map=np.sum(cubedata[1050:1330,2:67,:],axis=0)
        wavelengthrange=cube_wavelength[1050:1330]
    else:
        map = np.sum(cubedata[:, 2:67, :], axis=0)
        wavelengthrange=cube_wavelength

    if waveprint:
        print(wavelengthrange)


    return map

def Colormap(map, **krg):
    '''
    plot the map
    :param map: array of the map
    :return: none
    '''

    # set the axis and color
    ra = np.linspace(krg['ra'][1], krg['ra'][0], np.shape(map)[0])
    dec = np.linspace(krg['dec'][0], krg['dec'][1], np.shape(map)[1])
    X, Y = np.meshgrid(dec, ra)
    Y = Y[::-1, :]
    bounds = np.linspace(np.min(map), np.max(map), 25)
    norm = clr.BoundaryNorm(boundaries=bounds, ncolors=310)

    # plot it
    fig = plt.figure(1)
    ax = plt.gca()
    ax.set_aspect(1)
    ax.invert_xaxis()

    pcm = ax.pcolormesh(X, Y, map[::-1, :], norm=norm, cmap='RdBu_r')
    fig.colorbar(pcm, ax=ax, extend='both', orientation='vertical')
    plt.grid(axis='x')
    plt.grid(axis='y')
    # plt.xlabel(r'$wavelength$')
    # plt.ylabel(r'$distance$')
    # plt.axis('off')
    plt.show()
    return None