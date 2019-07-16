from cube import CubeData, Cubegenerator
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from matplotlib import pyplot as plt
from astropy import units as u
from pylab import pcolor
import numpy as np
import Cosmology
import IO


def Cubepre(path,cubename):
    '''
    read the datacube
    :param path: location of the IFU data
    :param cubename: its name
    :return: cube, wavelength axis, wcs
    '''

    #read the cube
    cube_data, cube_wavelength, wcs, cube = CubeData.ReadCube(path, cubename)

    #add unit to this cube
    flux_cube=cube_data*1e-16*u.erg/(u.s*(u.cm**2)*u.AA)
    return flux_cube,cube_wavelength,wcs

def Mapspectral(waverange=[3800.,4200.]):
    '''
    plot the spectra for the whole datacube, choose a proper size and generate the spectra of this position then move to the neighbouring postion,
    step by step then map the whole datacube
    :param waverange: wavelength interval within which plot the spectra
    :return: None
    '''


    #read the datacube
    flux,wavelength,wcs=Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits')

    #generate spectra for each postion, spectra_list is a list contains spectra for each postion
    spectral_list,position_list=SpectralData.Mapspectral(flux,[7,2])

    #plot and save the image of spectra
    for i in range(len(spectral_list)):

        #reduce the noise
        spectral_list[i]=ImgInterSmo.NoiseFilter(spectral_list[i],5,0.2)
        plt.figure(i,figsize=(17,7))
        img=SpectralData.SpectralPlot(onedspectral=spectral_list[i],wavelength=wavelength.value,wavecut=waverange,title=str(position_list[i][1])+','+str(position_list[i][0]))
        lowerline,upperline=SpectralData.WavelengthSelectPlot(spectral_list[i],waverange)
        plt.savefig(str(i)+'.png')
    return None

def Indispectral(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',position=[0,0],size=[5,2],cutinterval=[4020.,4028.],normalize=False):
    '''
    generate the spectra for a certain postion and wavelength interval
    :param position: selected position
    :param size: region size used to generate the spectra, in unit of pixel
    :param cutinterval: wavelength interval
    :return: img of spectra  with two red line parallel to y axis
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre(path, filename)
    ra,dec,=CubeData.WCS(wcs)

    #generate the spectra
    onedspectral=SpectralData.SourceSpectral(flux,size,position,ra.value,dec.value)
    lowlimit,uplimit=SpectralData.Findwavelength(wavelength.value,cutinterval)
    # if normalize is True:
    #     onedspectral=onedspectral/np.max(onedspectral)
    # plt.figure('spectra',figsize=(17, 7))
    # img=SpectralData.SpectralPlot(onedspectral,wavelength.value,cutinterval)
    # # lowerline,upperline=SpectralData.WavelengthSelectPlot(onedspectral.value,cutinterval)
    # return img#lowerline,upperline

    return onedspectral[lowlimit:uplimit]/np.max(onedspectral[lowlimit:uplimit]), wavelength.value[lowlimit:uplimit]

def Indiimg(path=None,filename=None,wavelengthcut=[4020.,4028.],internum=None):
    '''
    plot the 2D image for a selected wavelength interval
    :param wavelengthcut: wavelength interval
    :return: 2D array of the image and its ra and dec
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre(path, filename)
    # flux=CubeData.Cubebadpixelremovor(flux)
    RA,DEC=CubeData.WCS(wcs)
    wavelengthcut=np.array(wavelengthcut)
    wavelengthcut_conti=np.array([wavelengthcut[0]-15.,wavelengthcut[1]-(wavelengthcut[1]-wavelengthcut[0])])

    #generate the 2D image
    twodimg,dec,ra=ImagPlot.Cutmap(flux.value,wavelength.value,wavelengthcut,DEC.value,RA.value)
    twodimg_conti,_,_=ImagPlot.Cutmap(flux.value,wavelength.value,wavelengthcut_conti,DEC.value,RA.value)
    twodimg_conti=twodimg_conti*((wavelengthcut[1]-wavelengthcut[0])/(wavelengthcut_conti[1]-wavelengthcut_conti[0]))

    #subtract the continuum component
    twodimg=twodimg-twodimg_conti

    #interpolate and smooth
    twodimg = ImgInterSmo.ImgSmoothor(twodimg, [1.5, 0.428])
    if internum is not None:
        twodimg=ImgInterSmo.Arrayinterpolation(twodimg,internum)


    return twodimg,ra,dec

def Mapimg(internum,cutrange):
    '''
    plot the images for the whole wavelength, devide the wavelength to n intervals and plot image for each interval
    :param internum: interpolate and smooth the image, this variable control how smooth it is
    :param cutrange: the width of each interval
    :return: None
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits')
    ra,dec=CubeData.WCS(wcs)

    #to match the interpolated image we should also interpolate the ra,dec because these are the x-axis and y-axis
    ra,dec=ImgInterSmo.Arrayinterpolation(ra.value,internum)[1:-1],ImgInterSmo.Arrayinterpolation(dec.value,internum)[1:-1]

    #cut the datacube for each interval
    for i in range(int(1600/cutrange)):
        cutflux, cutwavelength = CubeData.CubeCut(flux, wavelength, 'manual', [cutrange*i,cutrange*(i+1)])

        #generate the 2D image for each interval and smooth it
        twodimg = (ImagPlot.PlotMap(cutflux)).value
        twodimg=ImgInterSmo.InSmimg(twodimg,internum,[3.,3.])
        print(cutwavelength)

        #plot and save this images
        ImagPlot.Twodplotimg(twodimg,dec,ra)
        # plt.show()
        plt.savefig('map'+str(int(np.median(cutwavelength.value)))+'.png')
    return None

def Sourcecheck(position,wavelengthcut,mark):
    '''
    plot the image of selected wavelength interval and plot the spectra of the selected position
    :param position: ra, dec of the selected source
    :param wavelengthcut: wavelength interval
    :param mark: if it's 2 then plot the 2D image, else plot the 3D image
    :return: None
    '''


    #generate the spectral of the position, the variable spectra is an numpy array, each cell stand for a wavelength, the value is the intensity of this wavelength
    spectral,upperline,lowline=Indispectral(position,[5,3],wavelengthcut)

    #generate the image within the wavelength interval
    twodimg,ra,dec=Indiimg(wavelengthcut)

    #mark the selected pixels(check if it's the right position)
    twodimg = CubeData.MarkSource(twodimg, ra, dec, position)
    position = CubeData.Coordinateconvert(position)
    print(position)

    #plot 2D or 3D image
    if mark is '2':
        img=ImagPlot.Twodplotimg(twodimg,dec,ra)
    else:
        img=ImagPlot.Threedplotimg(twodimg,dec,ra)

    plt.show()

def AbsorpFinder(abrange,upperrange,lowerrange):
    flux, wavelength, wcs = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits')
    ra, dec = CubeData.WCS(wcs)
    map,ummap,lmmap,dec,ra=CubeData.Findabsorption(flux.value,wavelength.value,abrange,upperrange,lowerrange,dec.value,ra.value)
    umimg=ImagPlot.Twodplotimg(ummap,dec,ra)
    lmimg=ImagPlot.Twodplotimg(lmmap,dec,ra)
    img=ImagPlot.Twodplotimg(map,dec,ra)
    plt.show()

def Indiimgreadout(path=None,filename=None,waveinterval=None,name=None):
    '''
    generate the fits file
    :param waveinterval: wavelength interval selected from the datacube
    :param name: name of the output fits file
    :return: None
    '''

    #generate the image of the selected wavelength interval
    twodimg,ra,dec=Indiimg(path=path,filename=filename,wavelengthcut=waveinterval)

    #read the header
    header=CubeData.Readheader(path,filename)

    #export the fits file
    CubeData.Imgreadout(twodimg,header,name)
    return None

def SingleContour(baseimgname,contouimgname,levels=None,col='red'):
    '''
    overlay a contour to an image
    :param baseimgname: contour overlay on this image
    :param contouimgname: image used to extract the contour
    :param levels: value of each contour
    :param col: color of the contour
    :return: None
    '''

    #read the baseimage and contouimage
    basimg,basheader,baswcs=IO.Accessfits(baseimgname)
    contouimg,contouheader,contouwcs=IO.Accessfits(contouimgname)


    vscale=600.

    #plot the image
    ax=plt.subplot()
    ax.contour(contouimg,levels=levels,linewidths=.5,color=col)
    plt.imshow(basimg[4:67,:],cmap='jet',vmax=vscale,vmin=-vscale)
    plt.axis('auto')
    plt.axis('off')
    cbar=plt.colorbar()
    cbar.set_label(r'$velocity (km/s) $')

def contouroverlay():
    hstimg,hstheader,hstwcs=IO.Accessfits('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammothicut.fits',0)
    baseimg= np.full(np.shape(hstimg), np.nan)
    lymanimg,lyheader,lywcs=IO.Accessfits('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/lyman.fits')
    lywcs=lywcs.dropaxis(2)
    heiiimg,heheader,hewcs=IO.Accessfits('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/heii.fits')
    heiiwcs=hewcs.dropaxis(2)
    civimg,civheader,civwcs=IO.Accessfits('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/civ.fits')
    civwcs=civwcs.dropaxis(2)
    coaimg,coaheader,coawcs=IO.Accessfits("/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galA.fits")
    coawcs=coawcs.dropaxis(2)
    cobimg, cobheader, cobwcs = IO.Accessfits(
        "/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galB.fits")
    cobwcs = cobwcs.dropaxis(2)
    cocimg, cocheader, cocwcs = IO.Accessfits(
        "/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galC.fits")
    cocwcs = cocwcs.dropaxis(2)
    codimg, codheader, codwcs = IO.Accessfits(
        "/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galD.fits")
    codwcs = codwcs.dropaxis(2)

    norm=ImagPlot.Scaleimgconverter(hstimg)
    ax=plt.subplot(projection=hstwcs)

    # extract the contours and overlay them on the base image
    contourlyman,ax= ImagPlot.Contourgenerator(ax, lymanimg, lywcs,[1.e-17,4.1667e-17,7.1334e-17,9.11001e-17],'red')#
    contourheii,ax = ImagPlot.Contourgenerator(ax, heiiimg, heiiwcs,[4.e-18,9.23278e-18,1.1e-17,2.e-17,3.e-17],'cyan')#
    contourciv,ax = ImagPlot.Contourgenerator(ax, civimg, civwcs,[4.e-18,1.0025e-17,1.305e-17,1.6075e-17,1.91e-17], 'lime')#
    contourco1,ax=ImagPlot.Contourgenerator(ax,coaimg[0],coawcs,[0.0123874,0.0155811,
                                                      0.0187747,0.0219684,0.0251621],'white')#[0.006,0.00919368,0.0123874,0.0155811,0.0187747,0.0219684,0.0251621]
    contourco2,ax = ImagPlot.Contourgenerator(ax, cobimg[0], cobwcs,
                                   [0.00817988,0.0118598,0.0155396,0.0192195], 'white')#[0.0045,0.00817988,0.0118598,0.0155396,0.0192195]
    contourco3,ax = ImagPlot.Contourgenerator(ax, cocimg[0], cocwcs,
                                   [0.00697652,0.00995303,0.0129295,0.0159061], 'white')#[0.004,0.00697652,0.00995303,0.0129295,0.0159061]
    contourc04,ax = ImagPlot.Contourgenerator(ax, codimg[0], codwcs,
                                   [0.007,0.0100019,0.0130039], 'white')

    ax.imshow(hstimg,norm=norm,cmap='gray')
    # plt.tight_layout(rect=[0,0,.92,1])
    plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.1)
    lines=[contourlyman.collections[0],contourheii.collections[0],contourciv.collections[0],contourco1.collections[0]]
    labels=['$Lyman\alpha$','$HeII$','$CIV$','$CO(1-0)$']
    plt.legend(lines,labels)
    plt.xlabel(r'$Right Ascension$')
    plt.ylabel(r'$declination$')
    plt.show()
    return None

def velocitymap(path=None,filename=None,z=2.3093,lamda0=None,waveinterval=[4010.,4050.],internum=None,mask_sig1=0.,mask_sig2=2.):
    '''
    calculate the flux-weighted velocity map and plot it
    :param z: redshift use to calculate the observed wavelength lambda0  (lambda-lambda0)*c/lambda0
    :param waveinterval: wavelength interval used to calculate velocity
    :param veloscale: this is the up limit of the velocity map, velocity higher than this value will be replaced with it, -veloscale is the low limit
    :return: None
    '''
    flux, wavelength, wcs = Cubepre(path, filename)#read the datacube, flux is cube
    header = CubeData.Readheader(path, filename)
    waveinterval=np.array(waveinterval)

    #convert ra,dec to angle distance from the mammoth-1 source
    ra, dec, = CubeData.WCS(wcs)
    ra_dis,dec_dis=CubeData.Angle2distance(ra,dec,[220.3519792,40.052625])
    ra_dis,dec_dis=ra_dis.value,dec_dis.value

    #cut the datacube along the axis of wavelength, only keep the cube within the waveinterval
    rangeflux,rangewavelength=CubeData.CubeCut(flux.value,wavelength.value,'manual',waveinterval)
    badflux,_=CubeData.CubeCut(flux.value,wavelength.value,'manual',waveinterval-100)
    continumflux,_=CubeData.CubeCut(flux.value,wavelength.value,'manual',waveinterval-(waveinterval[1]-waveinterval[0]))


    #convert wavelength to velocity
    velocity=SpectralData.wavelength2velocity(rangewavelength,z,lamda0)

    #subtract continuum component, remove bad pixels
    rangeflux_subtracted=CubeData.Continumsubtractor(continumflux,rangeflux)
    bad_sigma=CubeData.BadvalueEstimator(badflux[:,20:30,4:8])
    rangeflux_badpix = CubeData.Cubebadpixelremovor(rangeflux_subtracted, sigma=bad_sigma)


    # calculate the emission line within seeing, reduce the noise of each slice
    cube_velocity = ImgInterSmo.CubeSmooth(rangeflux_badpix, [3., 0.9])#


    # interpolate and smooth cube
    if internum is not None:
        cube_velocity,ra_dis,dec_dis=ImgInterSmo.CubeInterpolation(cube_velocity,ra_dis,dec_dis,internum)


    # do the optimal extraction, we generate two mask cubes for this data cube
    # the first mask is used to remove the influence of the background, because in the process
    # some pixels who have no emission line will be kept and this will influence the velocity map
    # the second mask is used to select the emission region.
    cube_velocity=CubeData.CubeNoiseFilter(cube_velocity, 3, .2)
    mask_cube1=CubeData.Maskgenerator(cube_velocity,mask_sig1)
    mask_cube2=CubeData.Maskgenerator2(cube_velocity,mask_sig2)
    mask_cube=mask_cube2*mask_cube1

    # convert the flux image to velocity map
    velomap,fluxmap=CubeData.Cubeweightedmean(cube_velocity,velocity,mask_cube)
    velomap=velomap/1e3



    #replace the velocity of pixels which have too large velocity
    # velomap[np.where(abs(velomap)>1500.)]=np.nan
    velomap[np.where(velomap<-2000.)]=np.nan
    velomap[np.where(velomap > 2000.)] =np.nan
    # CubeData.Imgreadout(velomap,header,'velocitymap.fits')

    return velomap,fluxmap,ra_dis,dec_dis

def velocitymap_filter(path=None,filename=None,z=2.3093,lamda0=None,waveinterval=[3990.,4060.],sigma_num=0.25,internum=None,mask_sig1=0.,mask_sig2=2.):
    '''
    show the velocity map which is only within the emission region
    :param path: path of which files used
    :param filename: data cube used to calculate the velocity map, it's transfered to function velocitymap
    :param z: redshift corresponding to the emission region
    :param lamda0: intrinsic wavelength of this emission line
    :param waveinterval: wavelength range within which calculated the velocity map, it's transfered to function velocitymap
    :param sigma_num: used  to keep pxiels whose value is beyond this parameter
    :param internum: used to control the interpolation, it's transfered to function velocitymap
    :return: velocity map only emission region kept
    '''

    # produce the velocity map
    velomap,fluxmap,ra_dis,dec_dis=velocitymap(path,filename,z,lamda0,waveinterval,internum,mask_sig1,mask_sig2)
    SNmap=CubeData.SNmapgenerator(fluxmap,waveinterval,1e3,600)


    ##

    # do the filter, only keep the pixels we need
    velomap_filtered=CubeData.Regionfilter(SNmap,velomap,sigma_num)
    velo_shape=np.shape(velomap_filtered)


    if lamda0==1215.673:
        #for lyman-alpha
        velomap_filtered[:,:int(0.15*velo_shape[0])]=np.nan
        velomap_filtered[:,-1-int(0.07*velo_shape[0]):]=np.nan
        velomap_filtered[:int(0.1*velo_shape[0]),:]=np.nan
        velomap_filtered[-1-int(0.1 * velo_shape[0]):, :] = np.nan
        # fluxmap[:, :int(0.1 * velo_shape[0])] = np.nan
        # fluxmap[:, -1 - int(0.1 * velo_shape[0]):] = np.nan
        # fluxmap[:int(0.01 * velo_shape[0]), :] = np.nan
        # fluxmap[-1 - int(0.1 * velo_shape[0]):, :] = np.nan
    elif lamda0==1640:
        # for heii
        velomap_filtered[:, :int(0.7 * velo_shape[0])] = np.nan
        velomap_filtered[:, -1 - int(0.25 * velo_shape[0]):] = np.nan
        velomap_filtered[:int(0.25 * velo_shape[0]), :] = np.nan
        velomap_filtered[-1 - int(0.15 * velo_shape[0]):, :] = np.nan
        # fluxmap[:, :int(0.7 * velo_shape[0])] = np.nan
        # fluxmap[:, -1 - int(0.25 * velo_shape[0]):] = np.nan
        # fluxmap[:int(0.3 * velo_shape[0]), :] = np.nan
        # fluxmap[-1 - int(0.25 * velo_shape[0]):, :] = np.nan
    else:
        #for civ
        velomap_filtered[:, :int(0.5 * velo_shape[0])] = np.nan
        velomap_filtered[:, -1 - int(0.2 * velo_shape[0]):] = np.nan
        velomap_filtered[:int(0.2 * velo_shape[0]), :] = np.nan
        velomap_filtered[-1 - int(0.2 * velo_shape[0]):, :] = np.nan
        # fluxmap[:, :int(0.75 * velo_shape[0])] = np.nan
        # fluxmap[:, -1 - int(0.27 * velo_shape[0]):] = np.nan
        # fluxmap[:int(0.3 * velo_shape[0]), :] = np.nan
        # fluxmap[-1 - int(0.25 * velo_shape[0]):, :] = np.nan


    return velomap_filtered,fluxmap,SNmap,ra_dis,dec_dis

def slitspectra(position):
    '''
    plot the 2D slit spectra
    :param position: ra,dec of the mammoth-1
    :return: None
    '''

    #read the datacube
    flux, wavelength, wcs = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits')
    ra, dec, = CubeData.WCS(wcs)

    #convert the ra,dec of mammoth-1 to image coordinate
    ra_img, dec_img = CubeData.Findposition(ra.value, dec.value,
                                            position)  # position in the form of [ra, dec] in unit of deg, physical position returned is also in this form

    dec_init=dec_img-9
    # calculate the emission line within seeing
    flux=CubeData.Cubeseeinglimit(flux.value)

    #find the location of up limit and low limit in the wavelength array
    wavelow,waveup=SpectralData.Findwavelength(wavelength.value,[4000.,4050.])

    #convert wavelength to velocity
    velocity = SpectralData.wavelength2velocity(wavelength[wavelow:waveup].value, 2.310, 1215.673)

    #convert ra to angle distance from mammoth-1
    delta_arcsec=(ra.value-position[0])*3600.

    #plot the spectra
    fig=plt.figure()
    for i in range(1,17):
        twodspectra=flux[wavelow:waveup,:,dec_init+i-1]#for each dec plot the 2D spectra
        print(dec_init+i-1)
        ax=fig.add_subplot(4,4,i)
        pcolor(velocity/1e3,(ra.value-position[0])*3600.,twodspectra.T,cmap='jet')

        # ax.imshow(twodspectra.T,cmap=cm.jet)

    fig.subplots_adjust(wspace=0., hspace=0.)
    fig.text(0.5, 0.04, r'$velocity(km/s)$', ha='center', va='center')
    fig.text(0.06, 0.5, r'$\Delta arcsec$', ha='center', va='center', rotation='vertical')

    plt.show(block=True)

    return None


def emissionlinegenerator(wavelength_c,width,peakflux,wavelengthrange):
    '''
    with center wavelength, line width and peak intensity, we use gaussian function to rebuild the
    emission line
    :param wavelength_c: central wavelength of the emission line
    :param width: numpy array, line width for each emission line
    :param peakflux: numpy array, peak intensity for each emission line
    :param wavelengthrange: numpy array, wavelength range
    :return: spectra of a emission line
    '''

    #convert line width to sigma for gaussian function
    sigma=width/(2*np.sqrt(2*np.log(2)))
    emissionline=Cubegenerator.gaussian(wavelengthrange,wavelength_c,sigma,peakflux)
    return emissionline


def emissionspectragenerator(wavelength,width,peakflux):
    '''
    this function generate the sky emission spectra
    :param wavelength: wavelength list for each emission line
    :param width: width list for each emission line
    :param peakflux: peak flux for each emission line
    :return: spectra with all emission line
    '''
    spectra=np.zeros(np.shape(wavelength))
    for i in range(len(wavelength)):
        spectra=emissionlinegenerator(wavelength[i],width[i],peakflux[i],wavelength)+spectra
    return spectra

def Run():
    # Mapspectral(waverange=[4000.,4050.])
    #  Run_indispectral([25,18])
    # Run_indiimg()
    # Indispectral([220.3491,40.0525])
    # Run_img([100,170])
    # Run_img([0,100])
    # Mapimg(internum=3,cutrange=8)
    # Indiimg(mark='2')

    #check the wavelength calibration of this two data cube

    #read the cube data, extract the spectra from it and also smooth the spectra
    spectra1,wavelength1=Indispectral(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/MAMMOHT-1_individual',
                                      filename='1441+4003_00136_icubes_cut.fits',size=[20,10],
                                      position=[220.3499625,40.04884167],cutinterval=[3600.,4300.],normalize=True)
    spectra1=ImgInterSmo.NoiseFilter(spectra1,3,.1)
    spectra2, wavelength2 = Indispectral(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/MAMMOHT-1_individual',
                                         filename='1441+4003_00145_icubes_cut.fits',size=[20,10],
                                         position=[220.3499625, 40.04884167], cutinterval=[5100., 5600.], normalize=True)
    spectra2 = ImgInterSmo.NoiseFilter(spectra2, 3, .1)

    #read the sky emission line template, generate the template spectra and also smooth it.
    skyspec5 = IO.Read_dat(
        '/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/sky emission line/J_A+A_407_1157/table5.dat')
    wavelength5=skyspec5[:-150,1]-4
    peakflux5 = skyspec5[:-150, 3] / np.max(skyspec5[:-150, 3])
    width5=skyspec5[:-150,2]
    skyspectra5=emissionspectragenerator(wavelength5,width5,peakflux5)
    skyspectra5=ImgInterSmo.NoiseFilter(skyspectra5,10,.1)
    skyspec6 = IO.Read_dat(
        '/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/sky emission line/J_A+A_407_1157/table6.dat')
    wavelength6=skyspec6[10:-20,1]
    peakflux6 = skyspec6[10:-20, 3]/np.max(skyspec6[10:-20,3])
    width6=skyspec6[10:-20,2]
    skyspectra6=emissionspectragenerator(wavelength6,width6,peakflux6)

    #plot the results
    fig,AX=plt.subplots(2,1)
    AX=AX.flatten()
    AX[0].plot(wavelength1,spectra1,label='cube spectra')
    AX[0].plot(wavelength5,skyspectra5,c='orange',label='sky emission line')
    AX[1].plot(wavelength2, spectra2,label='cube spectra')
    AX[1].plot(wavelength6, skyspectra6, c='orange',label='sky emission line')
    AX[1].tick_params(labelsize=15.)
    AX[0].tick_params(labelsize=15.)
    AX[0].legend(fontsize=15.)
    AX[1].legend(fontsize=15.)
    fig.text(0.5, 0.03, 'Wavelength($\AA$)', ha='center', fontsize=25.)
    fig.text(0.05, 0.5, 'Normalized Intensity', va='center', rotation='vertical', fontsize=25.)
    plt.show()




    # Sourcecheck([220.3522,40.0529],[4006.,4014.],mark='2')#source 1
    # Sourcecheck([220.3539,40.0536],[4006.,4014.],mark='2')#source 2
    # Sourcecheck([220.3491,40.0522],[4022.,4030.],mark='2')# source 3
    # Sourcecheck([220.34855,40.0525],[4022.,4030.],mark='2')# source 4
    # Sourcecheck([220.3489,40.0522],[4030.,4038.],mark='2')# source 5
    # Sourcecheck([220.3493,40.0525],[4029.,4038.],mark='2')# source 6
    # Sourcecheck([220.3519792,40.052625],[4000.,4030.],mark='2')
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL1',220.351083333+0.0015231999999798518)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL2',40.0515833333+0.0006535899999988715)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL1',220.351083333+0.0015231999999798518)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL2',40.0515833333+0.0006535899999988715)
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_ss_icubes.fits',waveinterval=[4010.,4050.],name='lyman.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5464.,5493.],name='heii.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5160.,5200.],name='civ.fits')
    # contouroverlay()
    # Sourcecheck([220.3505375,40.05144444],[5460.,5500.],mark='2')
    # Sourcecheck([220.3505375, 40.05144444], [5163., 5203.], mark='2')
    # Sourcecheck([220.3528625, 40.05370556], [5100., 5550.], mark='2')
    # lymanvelomap, lymanimg,lymanSN,ra_dis, dec_dis=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits',
    #                                                  z=2.310,lamda0=1215.673,waveinterval=[3990.,4060],sigma_num=6,internum=[2,8],mask_sig1=1.2e-19,mask_sig2=.35)#1.6e-19
    #
    # lymanimg=lymanimg*1e19
    # ImagPlot.Twodplotimg(map=[lymanvelomap,lymanimg],x=dec_dis,y=ra_dis,subclo=2,subrow=1,xlabel=r'arcsec',
    #                      ylabel=r'arcsec',cbarlabel=[r'$velocity(km/s)$','$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$'],
    #                      subtitle=['velocity map','flux map'],title='z=2.311',contourmap=lymanimg,contourlevel=[2,4.5,7.7,11.8,17])#[2,4.5,7.7,11.8,17]
    #
    # heiivelomap,heiiimg,heiiSN ,ra_dis, dec_dis=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',
    #                                      z=2.340,lamda0=1640,waveinterval=[5459.,5498.],sigma_num=3.8,internum=[2,8],mask_sig1=4.5e-20,mask_sig2=0.2)#.6e-19
    #
    # heiiimg=heiiimg*1e19
    # ImagPlot.Twodplotimg([heiivelomap, heiiimg], dec_dis, ra_dis, subclo=2, subrow=1, xlabel=r'arcsec',
    #                      ylabel=r'arcsec', cbarlabel=[r'$velocity(km/s)$','$flux(10^{-19} \ erg/s/cm^{2}/\AA)$'],
    #                      subtitle=['velocity map', 'flux map'],title='z=2.340',contourmap=heiiimg,contourlevel=[0.58,.9,1.3,1.8,2])
    # #
    # #
    # civvelomap, civimg,civSN,ra_dis, dec_dis=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',
    #                                     z=2.34385,lamda0=1549,waveinterval=[5160.,5200.],sigma_num=4.,internum=[2,8], mask_sig1=5.3e-20,mask_sig2=.2)
    #
    # civimg=civimg*1e19
    # ImagPlot.Twodplotimg([civvelomap, civimg], dec_dis, ra_dis, subclo=2, subrow=1, xlabel=r'arcsec',
    #                      ylabel=r'arcsec', cbarlabel=[r'$velocity(km/s)$','$flux(10^{-19} \ erg/s/cm^{2}/\AA)$'],
    #                      subtitle=['velocity map', 'flux map'],title='z=2.344',contourmap=civimg,contourlevel=[.653,.823,1.05,1.24,1.46,1.74])


    # datacube=Cubegenerator.Cubegenerator((100,130,140),[4000.,4060.])
    # mask=CubeData.Maskgenerator(datacube,2.)
    # velocity = SpectralData.wavelength2velocity(np.linspace(4000.,4060.,100), 2.311, 1215.673)
    # velomap,fluxmap=CubeData.Cubeweightedmean(datacube,velocity,mask)
    # velomap_filtered = CubeData.Regionfilter(fluxmap, velomap, 1000)
    # plt.imshow(velomap,cmap='jet')
    # # plt.imshow(fluxmap,cmap='jet')
    # plt.show()
    # scale=Cosmology.Angle_distance(2.307)*Cosmology.Arcsec2rad(20.)
    # print(scale)
    # a=1/(1+2.30)
    # scale=scale/a
    # print(scale)




    return None

Run()