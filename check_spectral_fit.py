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
# cube_name='1441+4003_comb_psfs_icubes.fits'

def ReadCube(path,cube_name):
    '''
    load the fits and print the basic information of this cube
    :param path: location of the fits
    :param cube_name: name of the fits
    :return: 3D data and wavelength axis
    '''

    os.chdir(path)
    cube = SpectralCube.read(cube_name)
    wcs=WCS(cube.header,cube.shape)
    print(cube)
    ra1,dec1=(wcs[:,:,0].T)[::-1,:],(wcs[:,:,1].T)[:,::-1]
    # ra=cube.world[0,:,:][2].value
    # dec=cube.world[0,:,:][1].value
    return cube._data,cube.spectral_axis,dec1,ra1,cube


def WCS(header,shape):
    '''
    read the keywords in header which are used to calculate the wcs and do the 3 steps to find the euqarotial coordinate for the pixels in the cube
    :param header: cube's header
    :param shape: cube's shape
    :return: equarotial wcs
    '''

    #read header
    referpoint_physical=np.array([header['CRPIX1'],header['CRPIX2']])
    referpoint_equarotial=np.array([header['CRVAL1'],header['CRVAL2']])
    phi_p=header['LONPOLE']*np.pi/180.
    transfer_matrix=np.array([[header['PC1_1'],header['PC1_2']],[header['PC2_1'],header['PC2_2']]])

    #create a n*m*2 array to hold coordinate for every pixel
    wcs_map=np.zeros((shape[2],shape[1],2))
    # put the physical coordinate of the pixel in corresponding cell which will be used to generate intermediate wcs
    for i in range(shape[2]):
        for j in range(shape[1]):
            wcs_map[i,j]=[i+1,j+1]

    #calculate intermediate wcs
    intermediate_wcs=Physical2Intermediate(wcs_map,transfer_matrix,referpoint_physical)
    #calculate sphere wcs with intermediate wcs
    sphere_wcs=DZenithalProjection(intermediate_wcs)
    #calculate equarotial wcs with sphere wcs
    equarotial_wcs=CoordinateRotation(sphere_wcs,referpoint_equarotial,phi_p)

    return equarotial_wcs

def Physical2Intermediate(wcs_map,transfer_matrix,referpoint_physical):
    '''
    convert physical coordinate to intermediate coordinates by transformation matrix in the header
    :param wcs_map: 3D array contain physical coordinate
    :param transfer_matrix: this matrix used to convert physical coordinate to intermediate coordinate
    :param referpoint_physical: physical coordinate of reference point
    :return: intermediate wcs
    '''
    deta_wcs_map=wcs_map-referpoint_physical
    intermediate_wcs=np.dot(deta_wcs_map,transfer_matrix.T)
    return intermediate_wcs

def DZenithalProjection(intermediate_wcs):
    '''
    convert intermediate coordinate to sphere coordinate
    natue of this step is porjecting 2D plane coordinate system to 2D sphere coordinate system
    :param intermediate_wcs:
    :return: sphere coordinate system
    '''
    intermediate_wcs_shape=np.shape(intermediate_wcs)
    sphere_wcs=np.zeros(intermediate_wcs_shape)

    #do the convertion for every pixels
    #this is zenithal projection which is a very simple projection(rules of projection can be found in wiki)
    for i in range(intermediate_wcs_shape[0]):
        for j in range(intermediate_wcs_shape[1]):
            r_theta=np.sqrt(np.sum(intermediate_wcs[i,j]**2))
            sphere_wcs[i, j] = [np.arctan2(-intermediate_wcs[i, j, 0],intermediate_wcs[i, j, 1]),
                                np.arctan(180. / (r_theta*np.pi))]
            # sphere_wcs[i, j] = [np.arctan2(-intermediate_wcs[i,j,0],intermediate_wcs[i,j,1]), np.arctan(180. / (r_theta * np.pi))]
    # sphere_wcs=sphere_wcs*180./np.pi
    return sphere_wcs

def CoordinateRotation(sphere_wcs,referpoint_equarotial,phi_p):
    '''
    this step rotates the axe of the sphere to the posstion of axes of equarotial coordinate system
    :param sphere_wcs:
    :param referpoint_equarotial: equarotial coordinate of reference point
    :param phi_p:
    :return: equarotial coordinate
    '''

    referpoint_equarotial=referpoint_equarotial*np.pi/180.
    sphere_wcs_shape = np.shape(sphere_wcs)
    equarotial_wcs=np.zeros(sphere_wcs_shape)

    for i in range(sphere_wcs_shape[0]):
        for j in range(sphere_wcs_shape[1]):
            tan=(np.cos(sphere_wcs[i,j,1])*np.sin(sphere_wcs[i,j,0]))/(np.sin(sphere_wcs[i,j,1])*np.cos(referpoint_equarotial[1])+np.cos(sphere_wcs[i,j,1])*np.cos(sphere_wcs[i,j,0])*np.sin(referpoint_equarotial[1]))
            sin=np.sin(sphere_wcs[i,j,1])*np.sin(referpoint_equarotial[1])+np.cos(sphere_wcs[i,j,1])*np.cos(referpoint_equarotial[1])*np.cos(sphere_wcs[i,j,0]-phi_p)
            equarotial_wcs[i,j]=[referpoint_equarotial[0]+np.arctan(tan),np.arcsin(sin)]
    equarotial_wcs[13,47]=referpoint_equarotial

    return equarotial_wcs*180./np.pi


def CubeCut(cube=None,cube_wavelength=None,ra=None,dec=None,emissionline=None,rang=None):
    '''
    cut the cube and return the part we need
    :param cube: cube data
    :param cube_wavelength: wavelength
    :param emissionline: select different part according to this parameter
    :param rang: if set emission line to manual,then this parameter must be given
    :return: part of the initial cube and corresponding wavelength
    '''
    if emissionline=='lyman':
        # cube_cut,wavelength_cut=cube[950:1010,2:67,:],cube_wavelength[950:1010]
        cube_cut, wavelength_cut = cube[940:990, 2:67, :], cube_wavelength[940:990]
    elif emissionline=='CV':
        # print(cube_wavelength[1030:1200])
        # cube_cut,wavelength_cut=cube[1030:1200,2:67,:],cube_wavelength[1030:1200]
        cube_cut, wavelength_cut =cube[1030:1500, 2:67, :], cube_wavelength[1030:1500]
    elif emissionline=='manual':
        cube_cut, wavelength_cut =cube[rang[0,0]:,rang[1,0]:rang[1,1],rang[2,0]:rang[2,1]]
    else:
        cube_cut, wavelength_cut =cube[:,2:67,:],cube_wavelength
    if ra is not None and dec is not None:
        ra_cut = ra[2:67, :]
        dec_cut = dec[2:67, :]

    else:
        ra_cut=None
        dec_cut=None
    return cube_cut,wavelength_cut,ra_cut,dec_cut



def PlotMap(cubedata):
    '''
    plot the 2D image of the cube by summing along the wavelength axis to get the total flux for each pixel
    :param cubedata: 3D cube data
    :return: 2D image array
    '''

    map=np.sum(cubedata,axis=0)

    return (map.T)[::-1,:]

def FLux(cube=None,wavelength_axis=None,ra=None,dec=None,gain=None,emissionline=None):
    """
            convert number of electrons per second of pixel to flux
            :param cube: cube data
            :param emissionline: this parameter determine which part of the given cube we need
            :param wavelength_axis: wavelength used to calculate flux
            :param gain: ccd gain number of electrons / number of photons
            :return: flux and corresponding wavelength
    """
    cube_all=CubeCut(cube=cube,cube_wavelength=wavelength_axis,ra=ra,dec=dec,emissionline=emissionline)
    cubedata,wavelength_cut=cube_all[0],cube_all[1]
    photons=cubedata/gain
    wavelength_axis=np.array((wavelength_cut).to(u.meter).data)
    for i in range(np.shape(wavelength_axis)[0]):
        photons[i]=photons[i]*const.h*(const.c.data/wavelength_axis[i])

    flux=photons

    return flux,wavelength_cut,cube_all[2],cube_all[3]


def ContinuumEst(flux,ran=[1300,1500]):
    '''
    estimate the continuum component of flux
    :param flux: total flux
    :param ran: wavelength range used to calculate continuum
    :return: continuum component
    '''

    continuum=np.mean(flux[ran[0]:ran[1],:,:],axis=0)
    continuum_std=np.std(flux[ran[0]:ran[1],:,:],axis=0)

    return continuum,continuum_std

def ContinuumSub(flux,flux_all):
    '''
    subtract continuum component
    :param flux: total flux
    :param flux_all: used to estimate continuum
    :return: flux and its standard deviation
    '''

    flux_sub=np.zeros(np.shape(flux))
    continuum,continuum_std=ContinuumEst(flux_all)
    for i in range(np.shape(flux)[0]):
        flux_sub[i,:,:]=flux[i,:,:]-continuum
    flux_std=continuum_std
    return flux,flux_std


def BackgroundEstimation(image,region):
    '''
    estimate background value of the given image
    :param image: wait to be analysis
    :param region: which part we use to estimate background
    :return: background and its standard deviation
    '''

    background=np.mean(image[region[0,0]:region[0,1],region[1,0]:region[1,1]])
    std_background=np.std(image[region[0,0]:region[0,1],region[1,0]:region[1,1]])

    return background,std_background

def ImageBackgroundSubtraction(image,background):
    '''
    subtract background from the given image
    :param image: wait to be subtracting
    :param background: background value
    :return: background-subrtacted image
    '''

    image=image-background
    return image


def CubeBackgroundSubtraction(cube):
    '''
    subtract background from the given cube
    :param cube: cube waited to be subtract
    :return: background-subtracted cube
    '''
    cube_bk=np.zeros(np.shape(cube))
    for i in range(len(cube)):
        background,background_std=BackgroundEstimation(cube[i],np.array([[11,20],[5,9]]))
        cube_bk[i]=ImageBackgroundSubtraction(cube[i],background)

    return cube_bk


def FindSource(map,FWHM=3.,sigma=3.):
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
    daofind=sourcefinder(fwhm=FWHM,threshold=sigma*photon_std)
    source=daofind(map-photon_mean)

    y,x=[int(i) for i in np.round(source['xcentroid'].data)],[int(i) for i in np.round(source['ycentroid'].data)]
    coordinate_map=np.vstack((x,y)).T

    return coordinate_map,boundary

def Spectral(cube,source_coordinate,wavelength):
    '''
    :param cube: cube data
    :param source_coordinate: coordinate of the sources
    :param wavelength: wavelength
    :return: None
    '''
    spectral_set=SourceSpectral(cube,source_coordinate)
    SpectralPlot(spectral_set,wavelength,source_coordinate)
    return None

def SourceSpectral(cube,source_coordinate):
    '''
    output the 1D spectral of selected sources
    :param cube: cube data which is a 3D numpy array
    :param source_coordinate: coordinate of selected sources which is a 2D numpy array
    :return: spectrals
    '''
    spectral_set=[]
    for individual in source_coordinate:
        cut=cube[:,individual[0]-1:individual[0]+2,individual[1]-1:individual[1]+2]#add neighboor pixels' value up to imporve SN
        spectral_individual=np.sum(np.sum(cut,axis=1),axis=1)
        spectral_set.append(spectral_individual)
    spectral_set=np.array(spectral_set)

    return spectral_set

def SpectralPlot(spectral_set,wavelength,source_coordinate):
    '''
    accept the spectral from SourceSpectral and plot them
    :param spectral_set: spectral of each sources
    :param wavelength: wavelength used to plot spectral
    :param source_coordinate: mark the the spectral to tell to which source it belongs
    :return: None
    '''
    n_source=len(spectral_set)
    wavelength=wavelength.value
    for i in range(1,n_source+1):
        plt.figure(i)
        plt.plot(wavelength,spectral_set[i-1])
        plt.text(x=0.99*np.median(wavelength),y=np.mean(spectral_set[i-1]),s=str(source_coordinate[i-1,0])+','+str(source_coordinate[i-1,1]))
        plt.xlabel(r'wavelength')
        plt.ylabel(r'flux')
        plt.title(r'spectral')
    plt.show()
    return None

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

def MapInterpolation(map,internum):
    '''
    interpolate value to the map
    :param map: original map
    :return: interpolated map
    '''

    #interploate original map along x axis and y axis
    for i in range(internum):
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
            # i_map=kernel_map[x-2:x+3,y-2:y+3]
            # pixel_value=kernel*kernel_map[x-2:x+3,y-2:y+3]
            # s=np.sum(kernel*kernel_map[x-2:x+3,y-2:y+3])
            new_map[x-2,y-2]=np.sum(kernel*kernel_map[x-2:x+3,y-2:y+3])

    return new_map[:-2,:-2]


def CubeInterpolation(cube,ra,dec,internum):
    '''
    interpolate the image of every wavelength in the cube
    :param cube: cube waiting interpolation
    :param internum: number of interpolation
    :return: interpolated cube
    '''
    shape0 = np.shape(cube)
    shape1 = np.shape(MapInterpolation(cube[0], internum))
    cube_inter = np.zeros((shape0[0], shape1[0], shape1[1]))
    for i in range(len(cube)):
        cube_inter[i]=MapInterpolation(cube[i],internum)
    ra_inter=MapInterpolation(ra,internum)
    dec_inter=MapInterpolation(dec,internum)
    return cube_inter,ra_inter,dec_inter

def CubeSmooth(cube):
    '''
    smooth the image of every wavelength in this cube
    :param cube: cube waiting smoothing
    :return: smoothed cube
    '''
    kernel=GussianKernel([3.,3.])
    shape0=np.shape(GussianFilter(cube[0],kernel))
    shape1=np.shape(cube)
    cube_sm=np.zeros((shape1[0],shape0[0],shape0[1]))
    for i in range(len(cube)):
        cube_sm[i]=GussianFilter(cube[i],kernel)

    return cube


def TwoDSpectral(cube,rang=None):
    '''
    extract the 2D spectra from the cube
    :param cube:
    :return: 2D spectral array
    '''
    cube=CubeCut(cube,emissionline='manual',rang=rang)
    twod_spectral=np.sum(cube,axis=2)
    return twod_spectral


def Colormap(map,**krg):
    '''
    plot the map
    :param map: array of the map
    :return: none
    '''

    #set the axis and color
    ra=np.linspace(krg['ra'][1],krg['ra'][0],np.shape(map)[0])
    dec=np.linspace(krg['dec'][0],krg['dec'][1],np.shape(map)[1])
    X,Y=np.meshgrid(dec,ra)
    Y=Y[::-1,:]
    bounds=np.linspace(np.min(map),np.max(map),25)
    norm = clr.BoundaryNorm(boundaries=bounds, ncolors=310)

    #plot it
    fig=plt.figure(1)
    ax=plt.gca()
    ax.set_aspect(1)
    ax.invert_xaxis()

    pcm=ax.pcolormesh(X,Y,map[::-1,:],norm=norm,cmap='RdBu_r')
    fig.colorbar(pcm, ax=ax, extend='both', orientation='vertical')
    plt.grid(axis='x')
    plt.grid(axis='y')
    # plt.xlabel(r'$wavelength$')
    # plt.ylabel(r'$distance$')
    # plt.axis('off')
    plt.show()
    return None

gain=0.145
cube_data,cube_wavelength,dec,ra,cube=ReadCube(path,cube_name)
flux_all,wavelength_all,ra_cut,dec_Cut=FLux(cube=cube_data,wavelength_axis=cube_wavelength,gain=gain)

# flux_all_sub,continuum_std=ContinuumSub(flux_all,flux_all)
# flux_all_sub_inter=CubeInterpolation(flux_all_sub,2)
# flux_all_sub_inter_sm=CubeSmooth(flux_all_sub_inter)
# flux_all_sub_inter_sm_map=PlotMap(flux_all_sub_inter_sm)
# Colormap(flux_all_sub_inter_sm_map)
flux_CV,wavelength_CV,ra_cut,dec_cut=FLux(cube=cube_data,wavelength_axis=cube_wavelength,ra=ra,dec=dec,gain=gain,emissionline='lyman')
flux_CV_sub,continuum_std=ContinuumSub(flux_CV,flux_all)
map=PlotMap(flux_CV_sub)
plt.imshow(map)
plt.show()
# coordinate_source=np.array([[32,11],[22,15],[39,14]])
# Spectral(flux_CV_sub,coordinate_source,wavelength_CV)
# flux_CV_sub_inter,ra_inter,dec_inter=CubeInterpolation(flux_CV_sub,ra_cut,dec_cut,4)
# flux_CV_sub_inter_sm=CubeSmooth(flux_CV_sub_inter)
# flux_CV_sub_inter_sm_map=PlotMap(flux_CV_sub_inter_sm)
# Colormap((flux_CV_sub_inter_sm_map.T)[::-1,:],dec=[np.mean(ra_inter[-1,:]),np.mean(ra_inter[0,:])],ra=[np.mean(dec_inter[:,0]),np.mean(dec_inter[:,-1])])

# coordinate_map,boundary=FindSource(smoothed_map,FWHM=4.6,sigma=1.2)#4.6 and 1.2 are the best value to select the three source
# coordinate_map,boundary=FindSource(flux_map,FWHM=3.,sigma=5.)
