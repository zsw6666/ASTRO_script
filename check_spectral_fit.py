import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from spectral_cube import SpectralCube
from astropy import constants as const
from astropy import units as u
from photutils import DAOStarFinder as sourcefinder
from astropy.io import fits
import os

path='/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI'
cube_name='1441+4003_comb_ss_icubes.fits'

def ReadCube(path,cube_name):
    '''
    load the fits and print the basic information of this cube
    :param path: location of the fits
    :param cube_name: name of the fits
    :return: 3D data and wavelength axis
    '''

    os.chdir(path)
    cube = SpectralCube.read(cube_name)
    print(cube)
    return cube._data,cube.spectral_axis,cube


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


    # plt.imshow(map)
    # plt.show()
    return map


def FLux(cubedata,wavelength_axis,gain):
    """
            convert number of electrons per second of pixel to flux
            :param cubedata: cube data
            :param wavelength_axis: wavelength used to calculate flux
            :param gain: ccd gain number of electrons / number of photons
            :return: map of flux
    """

    photons=cubedata[:,2:67,:]/gain
    wavelength_axis=np.array((wavelength_axis).to(u.meter).data)
    for i in range(np.shape(wavelength_axis)[0]):
        photons[i]=photons[i]*const.h*(const.c.data/wavelength_axis[i])

    flux=photons

    return flux


def ContinuumEst(flux,ran=[1300,1500]):
    '''
    estimate the continuum component of flux
    :param flux: total flux
    :param ran: wavelength range used to calculate continuum
    :return: continuum component
    '''

    continuum=np.mean(flux[ran[0]:ran[1],:,:],axis=0)

    return continuum

def ContinuumSub(flux,continuum):
    '''
    subtract continuum component
    :param flux: total flux
    :param continuum: comtinuum component
    :return: continuum subtracted flux
    '''

    flux_sub=np.zeros(np.shape(flux))
    for i in range(np.shape(flux)[0]):
        flux_sub[i,:,:]=flux[i,:,:]-continuum
    return flux


def FindSource(map,sigma=3.):
    '''
    this function is used to find the point-like source
    :param map: image's 2-dimensional array
    :param sigma: limit beyond which we consider as the source
    :return: source's coordinate in image
    '''
    # convert 2D array to 1D for convenience
    photon = map.flatten()

    # calculate the standard deviation and mean value of pixels
    photon_std = np.std(photon)
    photon_mean=np.mean(photon)
    boundary=photon_mean+sigma*photon_std

    #use the module already exist
    daofind=sourcefinder(fwhm=3.,threshold=sigma*photon_std)
    source=daofind(map-photon_mean)

    y,x=[int(i) for i in np.round(source['xcentroid'].data)],[int(i) for i in np.round(source['ycentroid'].data)]
    coordinate_map=np.vstack((x,y)).T

    return coordinate_map,boundary


def Spectral(cubedata,cube_wavelength,coordinate_map,image,emissionline=False):
    '''
    plot the spectral for each selected regions
    :param cubedata:
    :param cube_wavelength: wavelength axis
    :param coordinate_map: coordinates of selected pixels
    :param emissionline: if lyman==true then only search pixels which are bright at wavelength between 4000 to 4100
                  else search for pixels at the whole wavelength region
    :return:
    '''

    #for the whole wavelength region
    if emissionline==False:
        Spectral_set=IndividualPixelSpectral(cubedata,coordinate_map)
        for i in range(len(Spectral_set)):
            if i%3==0 and i>0:
                plt.figure('map')
                plt.imshow(image)
                plt.figure(i)
                plt.subplot(311)
                plt.plot(cube_wavelength[800:], Spectral_set[i - 1][800:])
                plt.title(str(coordinate_map[i-1][0]) + ',' + str(coordinate_map[i-1][1]))
                plt.subplot(312)
                plt.plot(cube_wavelength[800:], Spectral_set[i - 2][800:])
                plt.title(str(coordinate_map[i-2][0]) + ',' + str(coordinate_map[i-2][1]))
                plt.subplot(313)
                plt.plot(cube_wavelength[800:], Spectral_set[i - 3][800:])
                plt.title(str(coordinate_map[i-3][0]) + ',' + str(coordinate_map[i-3][1]))
                plt.show()
    #only for lyman-alpha wavelength region
    else:
        Spectral_set = EmissionLinePixelSepctral(cubedata, coordinate_map)
        print(coordinate_map)
        plt.figure(2)
        plt.plot(cube_wavelength[:], Spectral_set[:])
        plt.show()
    return None

def EmissionLinePixelSepctral(cubedata,coordinate_map):
    '''
    sum the flux for the lyman-alpha emitter at each wavelength
    :param cubedata: cube data
    :param coordinate_map: contain pixels which represent lyman-alpha region
    :return: spectral for lyman-alpha region
    '''

    flux=cubedata[:,coordinate_map[:,0],coordinate_map[:,1]]
    # plt.plot(range(len(flux[:, 10])), flux[:, 10])
    # plt.show()
    spectral=np.sum(flux,axis=1)
    return spectral

def IndividualPixelSpectral(cubedata,coordinate_map):
    '''
    sum the flux at each wavelength for selected pixels
    :param cubedata: cube data
    :param coordinate_map: pixels selected by FindSource
    :return: spectral for individual pixels
    '''

    Spectral_set=np.zeros(1600)
    for co in coordinate_map:
        flux=cubedata[:, co[0] - 2:co[0] + 2, co[1] - 2:co[1] + 2]
        # flux = cubedata[:, co[0] - 1:co[0] + 1, co[1] - 1:co[1] + 1]
        # flux = cubedata[:, co[0], co[1]]
        spectral = np.sum(np.sum(flux, axis=1), axis=1)
        Spectral_set = np.vstack((Spectral_set, spectral))
    Spectral_set = Spectral_set[1:]
    return Spectral_set


def LocateTarget(cube,target_coordinate):

    map=PlotMap(cube._data)
    wavelength,dec,ra=cube.world[:]
    deta_dec,deta_ra=np.abs(dec.value[0,:,:]-target_coordinate[1]),np.abs(ra.value[0,:,:]-target_coordinate[0])

    x,y=np.where(deta_ra==np.min(deta_ra))[0],np.where(deta_dec==np.min(deta_dec))[1]

    map[x,y]=10

    plt.imshow(map)
    plt.show()


def AxiesInterpolation(map,axies=0):
    '''
    interpolate value to the map
    :param map: original image
    :param axies: axies along which interpolated
    :return: image interpolated along chosen axis
    '''

    #choose which axis to interpolate along
    if axies==1:
        map=map.T

    #create new array waited values
    interpolated_map=np.zeros((2*np.shape(map)[0]-1,np.shape(map)[1]))
    #apply the mean value of the neighbor pixel
    inter_axies=(map[:-1,:]+map[1:,:])/2.

    #put the original image and interpolated image to the new image
    for axiesi in range(np.shape(interpolated_map)[0]):
        if axiesi%2==0:
            interpolated_map[axiesi,:]=map[int(axiesi/2),:]
        else:
            interpolated_map[axiesi,:]=inter_axies[int((axiesi-1)/2),:]

    #transverse the image
    if axies==1:
        final_map=interpolated_map.T
    else:
        final_map=interpolated_map

    return final_map

def MapInterpolation(map):
    '''
    interpolate value to the map
    :param map: original map
    :return: interpolated map
    '''

    #interploate original map along x axis and y axis
    for i in range(3):
        map=AxiesInterpolation(map)
        map=AxiesInterpolation(map,axies=1)
    map=map[::-1]

    return map


def GussianKernel(sigma):
    '''
    This is guassian kernel function (standard 2D gussian distribution)
    first generate N random numbers which obey gussian distribution
    and counts the number of points(Xcounts) in each x bins
    then define counts in each x bins as Ny and counts the number of points(Ycounts) in each y bins
    finally save this number in image

    :param sigma: standard deviation
    :return: 2D standard guassian kernel function (5*5)
    '''

    mux,muy=0,0
    scalex,scaley=sigma
    kernel=np.zeros((5,5)) # initilize the image
    Px=np.random.normal(loc=mux,scale=scalex,size=10000)  # N photons obey standard distribution
    Xcounts,Xbins=np.histogram(Px,np.linspace(-5,5,6))  # statistic number of photons in each X bins

    #for n in each x bins,generate n points obey normal distribution and save the value in image varaible
    Ybins=np.linspace(-5,5,6)
    for i in range(len(Xcounts)):
        Py=np.random.normal(loc=muy,scale=scaley,size=Xcounts[i])
        Ycounts,Ybins=np.histogram(Py,Ybins)
        kernel[:,i]=Ycounts # rewrite the image
    kernel=kernel/np.sum(kernel)


    return kernel

def GussianFilter(map,kernel):
    '''
    This function use gussian kernel to smooth the image

    :param map: original image
    :param kernel: gussian kernel
    :return: smoothed image
    '''

    #to calculate some pixels on the boundary of image,it first expands the original image
    kernel_map=np.zeros(np.array(np.shape(map))+np.array((4,4)))
    kernel_map[2:np.shape(kernel_map)[0]-2,2:np.shape(kernel_map)[1]-2]=map

    #create new map
    map_shape=np.shape(map)
    new_map=np.zeros(map_shape)

    #calculate the value for each pixel
    for x in range(2,map_shape[0]):
        for y in range(2,map_shape[1]):
            i_map=kernel_map[x-2:x+3,y-2:y+3]
            pixel_value=kernel*kernel_map[x-2:x+3,y-2:y+3]
            s=np.sum(kernel*kernel_map[x-2:x+3,y-2:y+3])
            new_map[x-2,y-2]=np.sum(kernel*kernel_map[x-2:x+3,y-2:y+3])

    return new_map[:-2,:-2]


def Colormap(map,boundary):
    '''
    plot the map
    :param map: array of the map
    :return: none
    '''

    #set the axis and color
    Y,X=np.mgrid[0:np.shape(map)[0],0:np.shape(map)[1]]
    bounds=np.linspace(np.min(map),np.max(map),10)
    norm = clr.BoundaryNorm(boundaries=bounds, ncolors=256)

    #plot it
    fig=plt.figure(1)
    plt.xlim(xmax=175.)
    ax=plt.gca()
    ax.set_aspect(1)
    pcm=ax.pcolormesh(X,Y,map,norm=norm,cmap='RdBu_r')
    plt.contour(X,Y,map,[boundary],alpha=.75,color='b')
    fig.colorbar(pcm, ax=ax, extend='both', orientation='vertical')
    plt.axis('off')
    plt.show()

gain=0.145
cube_data,cube_wavelength,cube=ReadCube(path,cube_name)
flux=FLux(cube_data,cube_wavelength,gain)
continuum=ContinuumEst(flux)
flux_sub=ContinuumSub(flux,continuum)
flux_map=PlotMap(cubedata=flux_sub,cube_wavelength=cube_wavelength,waveprint=False,emissionline='CV')
coordinate_map,boundary=FindSource(flux_map,5.)
# Spectral(flux_sub,cube_wavelength,coordinate_map,image=flux_map,emissionline=False)
# LocateTarget(cube,[220.3520458,40.05268333])
# interpolated_map=MapInterpolation(flux_map)
# kernel=GussianKernel([3.,3.])
# smoothed_map=GussianFilter(interpolated_map,kernel)
# Colormap(smoothed_map,boundary)