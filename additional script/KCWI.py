from cube import CubeData, Cubegenerator
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from matplotlib import pyplot as plt
from astropy import units as u
from pylab import pcolor
import numpy as np
import IO
import matplotlib.gridspec as grs
import pyspeckit


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
    ra,dec=CubeData.WCS(wcs)
    return flux_cube,cube_wavelength,ra,dec


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

def Indiimg(path=None,filename=None,wavelengthcut=[4020.,4028.],wavelengthcut_conti=[4123,4269],internum=None,smoothmark=None):
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

    #generate the 2D image
    twodimg,dec,ra=ImagPlot.Cutmap(flux.value,wavelength.value,wavelengthcut,DEC.value,RA.value)

    #subtract the continuum component
    if wavelengthcut_conti is not None:
        twodimg_conti, _, _ = ImagPlot.Cutmap(flux.value, wavelength.value, wavelengthcut_conti, DEC.value, RA.value)
        twodimg_conti = twodimg_conti * (
                    (wavelengthcut[1] - wavelengthcut[0]) / (wavelengthcut_conti[1] - wavelengthcut_conti[0]))
        twodimg=twodimg-twodimg_conti

    #interpolate and smooth
    if smoothmark is not None:
        twodimg = ImgInterSmo.ImgSmoothor(twodimg, [1.5, 0.428])
    if internum is not None:
        twodimg=ImgInterSmo.Arrayinterpolation(twodimg,internum)


    return twodimg,ra,dec

def Indiimgreadout(path=None,filename=None,waveinterval=None,wavelengthcut_conti=None,name=None,interpara=None,smoothpara=None):
    '''
    generate the fits file
    :param waveinterval: wavelength interval selected from the datacube
    :param name: name of the output fits file
    :return: None
    '''

    #generate the image of the selected wavelength interval
    twodimg,ra,dec=Indiimg(path=path,filename=filename,wavelengthcut=waveinterval,wavelengthcut_conti=wavelengthcut_conti,internum=interpara,smoothmark=smoothpara)

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
    lysourceimg,_,_=IO.Accessfits('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/lyman_source2.fits')
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
    contoursource, ax = ImagPlot.Contourgenerator(ax, lysourceimg, lywcs, [2.78194e-19,4.29945e-19,5.11695e-19,5.93445e-19 ],
                                                 'darkorange')  #
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
    lines=[contourlyman.collections[0],contoursource.collections[0],contourheii.collections[0],contourciv.collections[0],contourco1.collections[0]]
    labels=['$Lyman\alpha$','Source','$HeII$','$CIV$','$CO(1-0)$']
    plt.legend(lines,labels)
    plt.xlabel(r'$Right Ascension$')
    plt.ylabel(r'$declination$')
    plt.show()
    return None


def velopre(flux,wavelength,ra,dec,waveinterval=[4010.,4050.],internum=None,mask_sig2=2.):

    waveinterval = np.array(waveinterval)

    # convert ra,dec to angle distance from the mammoth-1 source
    ra_dis, dec_dis = CubeData.Angle2distance(ra, dec, [220.3519792, 40.052625])
    ra_dis, dec_dis = ra_dis.value, dec_dis.value

    # cut the datacube along the axis of wavelength, only keep the cube within the waveinterval
    rangeflux, rangewavelength = CubeData.CubeCut(flux.value, wavelength.value, 'manual', waveinterval)
    badflux, _ = CubeData.CubeCut(flux.value, wavelength.value, 'manual', waveinterval - 100)
    continumflux, _ = CubeData.CubeCut(flux.value, wavelength.value, 'manual',
                                       waveinterval - (waveinterval[1] - waveinterval[0]))

    # subtract continuum component, remove bad pixels
    rangeflux_subtracted = CubeData.Continumsubtractor(continumflux, rangeflux)
    # bad_sigma = CubeData.BadvalueEstimator(badflux[:, 20:30, 4:8])
    # rangeflux_badpix = CubeData.Cubebadpixelremovor(rangeflux_subtracted, sigma=bad_sigma)
    # rangeflux_badpix=rangeflux_subtracted


    if internum is not None:
        cube_velocity,ra_dis,dec_dis=ImgInterSmo.CubeInterpolation(rangeflux_subtracted,ra_dis,dec_dis,internum)
    cube_velocity = ImgInterSmo.CubeSmooth(cube_velocity, [3., 3.])  # [3., 0.9]
    map = np.mean(cube_velocity, axis=0)


    # do the optimal extraction, we generate two mask cubes for this data cube
    # the first mask is used to remove the influence of the background, because in the process
    # some pixels who have no emission line will be kept and this will influence the velocity map
    # the second mask is used to select the emission region.
    cube_velocity=CubeData.CubeNoiseFilter(cube_velocity, 3, .2)
    mask_cube=CubeData.Maskgenerator2(cube_velocity,mask_sig2)
    cube=cube_velocity*mask_cube
    return cube_velocity,mask_cube,rangewavelength,ra_dis,dec_dis


def velocitymap(flux,wavelength,ra,dec,z=2.3093,lamda0=None,waveinterval=[4010.,4050.],internum=None,mask_sig2=2.):
    '''
    calculate the flux-weighted velocity map and plot it
    :param z: redshift use to calculate the observed wavelength lambda0  (lambda-lambda0)*c/lambda0
    :param waveinterval: wavelength interval used to calculate velocity
    :param veloscale: this is the up limit of the velocity map, velocity higher than this value will be replaced with it, -veloscale is the low limit
    :return: None
    '''
    cube_velocity,mask_cube,rangewavelength,ra_dis,dec_dis=velopre(flux,wavelength,ra,dec,
                                                                   waveinterval=waveinterval,
                                                                   internum=internum,mask_sig2=mask_sig2)


    # convert wavelength to velocity
    # mask_cube=np.ones(np.shape(mask_cube))
    velocity = SpectralData.wavelength2velocity(rangewavelength, z, lamda0)
    # convert the flux image to velocity map and velocity dispersion map
    velomap=CubeData.Cubeweightedmean(cube_velocity,velocity,mask_cube)
    velodispmap=CubeData.Cubeweightstd(velocity,cube_velocity,mask_cube,velomap)/1e3
    velomap = velomap / 1e3

    #do the optimal extraction for data cube
    fluxmap=CubeData.OptimalextractImg(cube_velocity,mask_cube)



    #replace the velocity of pixels which have too large velocity
    # velomap[np.where(abs(velomap)>1500.)]=np.nan
    velomap[np.where(velomap<-2000.)]=np.nan
    velomap[np.where(velomap > 2000.)] =np.nan
    velodispmap[np.where(velodispmap > 1000.)] = np.nan
    # CubeData.Imgreadout(velomap,header,'velocitymap.fits')

    return velomap,velodispmap,fluxmap,ra_dis,dec_dis


def velocitymap_filter(flux,wavelength,ra,dec,z=2.3093,lamda0=None,
                       waveinterval=[3990.,4060.],
                       sigma_num=0.25,internum=None,mask_sig2=2.):
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
    velomap,velodisp,fluxmap,ra_dis,dec_dis=velocitymap(flux,wavelength,ra,dec,z,lamda0,waveinterval,internum,mask_sig2)
    SNmap=CubeData.SNmapgenerator(fluxmap,waveinterval,1e3,600)

    # do the filter, only keep the pixels we need
    velomap_filtered=CubeData.Regionfilter(SNmap,velomap,sigma_num)
    velo_shape=np.shape(velomap_filtered)


    if lamda0==1215.673:
        #for lyman-alpha
        velomap_filtered[:,:int(0.15*velo_shape[0])]=np.nan
        velomap_filtered[:,-1-int(0.07*velo_shape[0]):]=np.nan
        velomap_filtered[:int(0.1*velo_shape[0]),:]=np.nan
        velomap_filtered[-1-int(0.1 * velo_shape[0]):, :] = np.nan
    elif lamda0==1640:
        # for heii
        velomap_filtered[:, :int(0.7 * velo_shape[0])] = np.nan
        velomap_filtered[:, -1 - int(0.25 * velo_shape[0]):] = np.nan
        velomap_filtered[:int(0.25 * velo_shape[0]), :] = np.nan
        velomap_filtered[-1 - int(0.15 * velo_shape[0]):, :] = np.nan
    else:
        #for civ
        velomap_filtered[:, :int(0.5 * velo_shape[0])] = np.nan
        velomap_filtered[:, -1 - int(0.2 * velo_shape[0]):] = np.nan
        velomap_filtered[:int(0.2 * velo_shape[0]), :] = np.nan
        velomap_filtered[-1 - int(0.2 * velo_shape[0]):, :] = np.nan
    velodisp[np.where(np.isnan(velomap_filtered))]=np.nan



    return velomap_filtered,velodisp,fluxmap,SNmap,ra_dis,dec_dis

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

def OptimalExtractimg(path=None,filename=None,waveinterval=[4010.,4050.],continuuminterval=[4123.,4269.],internum=None,mask_sig1=0.):

    flux, wavelength, wcs = Cubepre(path, filename)  # read the datacube, flux is cube
    header = CubeData.Readheader(path, filename)
    waveinterval = np.array(waveinterval)

    # convert ra,dec to angle distance from the mammoth-1 source
    ra, dec, = CubeData.WCS(wcs)
    ra_dis, dec_dis = CubeData.Angle2distance(ra, dec, [220.3522, 40.0522])
    ra_dis, dec_dis = ra_dis.value, dec_dis.value

    # cut the datacube along the axis of wavelength, only keep the cube within the waveinterval


    rangeflux, rangewavelength = CubeData.CubeCut(flux.value, wavelength.value, 'manual', waveinterval)
    continumflux, _ = CubeData.CubeCut(flux.value, wavelength.value, 'manual',
                                       continuuminterval)

    # rangeflux_subtracted = CubeData.Continumsubtractor(continumflux, rangeflux)
    # rangeflux_badpix = CubeData.Cubebadpixelremovor(rangeflux_subtracted, sigma=bad_sigma)

    # calculate the emission line within seeing, reduce the noise of each slice
    # cube_velocity = ImgInterSmo.CubeSmooth(rangeflux_subtracted, [3., 0.9])  #
    cube_velocity, ra_dis, dec_dis = ImgInterSmo.CubeInterpolation(rangeflux, ra_dis, dec_dis, internum)
    cube_velocity = CubeData.CubeNoiseFilter(cube_velocity, 3, .2)
    # mask_cube = CubeData.Maskgenerator(cube_velocity, mask_sig1)
    mask_cube=np.ones(np.shape(cube_velocity))
    rangevelocity=SpectralData.wavelength2velocity(rangewavelength,2.31,1215.673)/1e3


    # img=np.mean(cube_velocity,axis=0)
    # img=ImagPlot.ImgCut(img, [0.15, 0.07, 0.1, 0.1])
    # SNmap=CubeData.SNmapgenerator(img,[waveinterval[0],waveinterval[1]],1e3,600)
    # img=CubeData.Regionfilter(SNmap,img,4.)
    # ImagPlot.Twodplotimg(img*1e19,dec_dis,ra_dis,subclo=1,subrow=1,xlabel='arcsec',ylabel='arcsec',
    #                      cbarlabel=['$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$'],subtitle=['$3700 \AA - 4300 \AA$'],markpoint=None)
    # scatter(1.8, -1.5, marker='*', color='midnightblue', s=100. ,label='mammoth-1')
    # scatter(3.3, -2.9, marker='*', color='gray', s=100., label='source-0')
    # scatter(2.1, -9, marker='*', color='magenta', s=100. ,label='source-1')
    # scatter(-2.0, -4.3, marker='*', color='darkolivegreen', s=100., label='source-2')
    # scatter(-7.5, -9.9, marker='*', color='ghostwhite', s=100., label='source-3')
    # plt.legend(fontsize=10.,loc='lower left')
    # plt.show()




    # velocity_list = []
    # map_list = []
    # for i in range(0,np.shape(cube_velocity)[0]-20,20):
    #
    #     map=CubeData.OptimalextractImg(cube_velocity[i:i+20,:,:],mask_cube[i:i+20,:,:])
    #     velocity_list.append(rangevelocity[i])
    #     SNmap=CubeData.SNmapgenerator(map,[rangewavelength[i],rangewavelength[i+1]],1e3,600)
    #     map=ImagPlot.ImgCut(map,[0.15,0.07,0.1,0.1])
    #     map=map*1e19
    #     map=CubeData.Regionfilter(SNmap,map,5.)
    #     map_list.append(map)
    # fig,AX=plt.subplots(2,round(len(map_list)/2),sharex=True, sharey=True)
    # fig.text(0.5, 0.05, 'arcsec', ha='center', fontsize=25.)
    # fig.text(0.08, 0.5, 'arcsec', va='center', rotation='vertical', fontsize=25.)
    # AX=AX.flatten()
    # for i in range(len(map_list)):
    #     AX[i].text(-15,7,'%d km/s'%velocity_list[i])
    #     img=AX[i].pcolor(dec_dis,ra_dis,map_list[i],cmap='gist_ncar')
    #     AX[i].scatter(2.1, -2, marker='*', color='midnightblue', s=100. ,label='mammoth-1')
    #     AX[i].scatter(2.1, -9, marker='*', color='magenta', s=100. ,label='source 1')
    #     AX[i].scatter(-2.0, -4.3, marker='*', color='darkolivegreen',s=100. , label='source 2')
    #     AX[i].legend(fontsize=10.,loc='lower left')
    #     AX[i].set_aspect('equal')
    #     # ImagPlot.Twodplotimg(map,dec_dis,ra_dis,subrow=1,subclo=1,xlabel='arcsec',ylabel='arcsec',cbarlabel=['$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$'],subtitle=[str(int(rangewavelength[i]))],contourmap=map,contourlevel=3)
    # cbar_ax = fig.add_axes([0.93, 0.15, 0.01, 0.7])
    # cbar = fig.colorbar(img, cax=cbar_ax,aspect=10)
    # cbar.set_label('$flux(10^{-19} \ erg/s/cm^{2}/\AA)$', fontsize=25.)
    # cbar.ax.tick_params(labelsize=15.)
    # fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()

    img=CubeData.OptimalextractImg(cube_velocity,mask_cube)
    SNmap=CubeData.SNmapgenerator(img,waveinterval,1e3,600)
    img_filtered=CubeData.Regionfilter(SNmap,img,5.)

    return img_filtered,ra_dis,dec_dis


def Run():

    # Indispectral([220.3491,40.0525])
    # Indiimg(mark='2')
    #
    # #check the wavelength calibration of this two data cube
    #
    #read the cube data, extract the spectra from it and also smooth the spectra
    # spectra1,wavelength1=Indispectral(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/MAMMOHT-1_individual',
    #                                   filename='1441+4003_00136_icubes_cut.fits',size=[20,10],
    #                                   position=[220.3499625,40.04884167],cutinterval=[3600.,4300.],normalize=True)
    # spectra1=ImgInterSmo.NoiseFilter(spectra1,3,.1)
    # spectra2, wavelength2 = Indispectral(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/MAMMOHT-1_individual',
    #                                      filename='1441+4003_00145_icubes_cut.fits',size=[20,10],
    #                                      position=[220.3499625, 40.04884167], cutinterval=[5100., 5600.], normalize=True)
    # spectra2 = ImgInterSmo.NoiseFilter(spectra2, 3, .1)
    #
    # #read the sky emission line template, generate the template spectra and also smooth it.
    # skyspec5 = IO.Read_dat(
    #     '/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/sky emission line/J_A+A_407_1157/table5.dat')
    # wavelength5=skyspec5[:-150,1]-4
    # peakflux5 = skyspec5[:-150, 3] / np.max(skyspec5[:-150, 3])
    # width5=skyspec5[:-150,2]
    # skyspectra5=emissionspectragenerator(wavelength5,width5,peakflux5)
    # skyspectra5=ImgInterSmo.NoiseFilter(skyspectra5,10,.1)
    # skyspec6 = IO.Read_dat(
    #     '/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/MMAMOTH1/sky emission line/J_A+A_407_1157/table6.dat')
    # wavelength6=skyspec6[10:-20,1]
    # peakflux6 = skyspec6[10:-20, 3]/np.max(skyspec6[10:-20,3])
    # width6=skyspec6[10:-20,2]
    # skyspectra6=emissionspectragenerator(wavelength6,width6,peakflux6)
    #
    # #plot the results
    # fig,AX=plt.subplots(2,1)
    # AX=AX.flatten()
    # AX[0].plot(wavelength1,spectra1,label='cube spectra')
    # AX[0].plot(wavelength5,skyspectra5,c='orange',label='sky emission line')
    # AX[1].plot(wavelength2, spectra2,label='cube spectra')
    # AX[1].plot(wavelength6, skyspectra6, c='orange',label='sky emission line')
    # AX[1].tick_params(labelsize=15.)
    # AX[0].tick_params(labelsize=15.)
    # AX[0].legend(fontsize=15.)
    # AX[1].legend(fontsize=15.)
    # fig.text(0.5, 0.03, 'Wavelength($\AA$)', ha='center', fontsize=25.)
    # fig.text(0.05, 0.5, 'Normalized Intensity', va='center', rotation='vertical', fontsize=25.)
    # plt.show()


    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL1',220.351083333+0.0015231999999798518-0.0002083999999911157)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL2',40.0515833333+0.0006535899999988715+0.000358470000001887-0.00024711999999738055)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL1',220.351083333+0.0015231999999798518-0.0002083999999911157)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL2',40.0515833333+0.0006535899999988715+0.000358470000001887-0.00024711999999738055)
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_ss_icubes.fits',
    #                waveinterval=[3700.,4300.],wavelengthcut_conti=None,name='coadd.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', filename='1441+4003_comb_ss_icubes.fits',
    #                waveinterval=[4037., 4047.], wavelengthcut_conti=None, name='ly3747.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5424,5434],name='heii424434.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5130.,5250.],name='civ130250.fits')
    # contouroverlay()

    # flux,wavelength,ra,dec=Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    # lymanvelomap, lymandisp,lymanimg,lymanSN,ra_dis, dec_dis=velocitymap_filter(flux,wavelength,ra,dec,z=2.308,
    #                                                                   lamda0=1215.673,waveinterval=[3990.,4060],
    #                                                                   sigma_num=6,internum=[2,8],mask_sig2=.35)#.18
    #
    # lymanimg=lymanimg*1e19
    # fig=plt.figure(constrained_layout=True)
    # gs=grs.GridSpec(2,2,fig,wspace=0,hspace=0)
    # ax0=fig.add_subplot(gs[:,0])
    # AX=[fig.add_subplot(gs[i,1])  for i in range(2)]
    # AX[-1].yaxis.set_ticks_position('right')
    # imglist=[lymanvelomap,lymandisp,lymanimg]
    # cbarlabel=[r'$velocity(km/s)$',r'$velocity(km/s)$','$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$']
    # for i in range(len(AX)-1):
    #     AX[i].get_shared_x_axes().join(AX[i], AX[-1])
    #     AX[i].get_shared_y_axes().join(AX[i], AX[-1])
    #     AX[i].set_xticklabels([])
    #     AX[i].set_yticklabels([])
    # AX.append(ax0)
    # ImagPlot.Gimgplot(fig,AX,imglist,dec_dis,ra_dis,lymanimg,[2,4.5,7.7,11.8,17],'arcsec','arcsec',
    #                   cbarlabel,['velocity','velocity dispersion','NB'])
    #
    # plt.show()


    flux, wavelength, ra, dec = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits')
    # heiivelomap,heiidisp,heiiimg,heiiSN ,ra_dis, dec_dis=velocitymap_filter(flux,wavelength,ra,dec,z=2.340,
    #                                                                         lamda0=1640,waveinterval=[5459.,5498.],
    #                                                                         sigma_num=3.8,internum=[2,8],mask_sig2=0.2)#.6e-19
    #
    # heiiimg=heiiimg*1e19
    # fig = plt.figure(constrained_layout=True)
    # gs = grs.GridSpec(2, 2, fig, wspace=0, hspace=0)
    # ax0 = fig.add_subplot(gs[:, 0])
    # AX = [fig.add_subplot(gs[i, 1]) for i in range(2)]
    # AX[-1].yaxis.set_ticks_position('right')
    # imglist = [heiivelomap, heiidisp, heiiimg]
    # cbarlabel = [r'$velocity(km/s)$', r'$velocity(km/s)$', '$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$']
    # for i in range(len(AX) - 1):
    #     AX[i].get_shared_x_axes().join(AX[i], AX[-1])
    #     AX[i].get_shared_y_axes().join(AX[i], AX[-1])
    #     AX[i].set_xticklabels([])
    #     AX[i].set_yticklabels([])
    # AX.append(ax0)
    # ImagPlot.Gimgplot(fig, AX, imglist, dec_dis, ra_dis, heiiimg, [0.58,.9,1.3,1.8,2], 'arcsec', 'arcsec',
    #                   cbarlabel, ['velocity', 'velocity dispersion', 'NB'])
    # plt.show()


    civvelomap, civdisp,civimg,civSN,ra_dis, dec_dis=velocitymap_filter(flux,wavelength,ra,dec,z=2.34385,lamda0=1549,
                                                                waveinterval=[5160.,5200.],sigma_num=3.8,internum=[2,8],mask_sig2=.2)

    civimg=civimg*1e19
    fig = plt.figure(constrained_layout=True)
    gs = grs.GridSpec(2, 2, fig, wspace=0, hspace=0)
    ax0 = fig.add_subplot(gs[:, 0])
    AX = [fig.add_subplot(gs[i, 1]) for i in range(2)]
    AX[-1].yaxis.set_ticks_position('right')
    imglist = [civvelomap, civdisp, civimg]
    cbarlabel = [r'$velocity(km/s)$', r'$velocity(km/s)$', '$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$']
    for i in range(len(AX) - 1):
        AX[i].get_shared_x_axes().join(AX[i], AX[-1])
        AX[i].get_shared_y_axes().join(AX[i], AX[-1])
        AX[i].set_xticklabels([])
        AX[i].set_yticklabels([])
    AX.append(ax0)
    ImagPlot.Gimgplot(fig, AX, imglist, dec_dis, ra_dis, civimg, [.603,.823,1.05,1.46,1.74], 'arcsec', 'arcsec',
                      cbarlabel, ['velocity', 'velocity dispersion', 'NB'])
    plt.show()



    # fig=plt.figure(constrained_layout=True)
    # gs=grs.GridSpec(2,4,fig,wspace=0,hspace=0)
    # ax0=fig.add_subplot(gs[:,:2])
    # AX=[fig.add_subplot(gs[i,j]) for j in range(2,4,1) for i in range(2)]
    # AX[-1].yaxis.set_ticks_position('right')
    # for i in range(len(AX)-1):
    #     AX[i].get_shared_x_axes().join(AX[i], AX[-1])
    #     AX[i].get_shared_y_axes().join(AX[i], AX[-1])
    #     AX[i].set_xticklabels([])
    #     AX[i].set_yticklabels([])
    #
    # AX.append(ax0)
    # fig.text(0.5, 0.06, 'arcsec', ha='center', fontsize=25.)
    # fig.text(0.08, 0.5, 'arcsec', va='center', rotation='vertical', fontsize=25.)
    #
    #
    # wavelist=[(4007+i*10,4007+(i+1)*10) for i in range(4)]
    # wavelist.append((3700,4300))
    # for i in range(len(wavelist)):
    #     img,ra,dec=OptimalExtractimg('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits',
    #                       waveinterval=[wavelist[i][0],wavelist[i][1]],internum=[2,8],mask_sig1=0.)
    #     img = ImagPlot.ImgCut(img, [0.1, 0.1, 0.1, 0.1])
    #     img=AX[i].pcolor(dec,ra,img,cmap='gist_ncar')
    #     AX[i].text(-15, 7, str(wavelist[i][0])+'$\AA$-'+str(wavelist[i][1])+'$\AA$',fontsize=17.)
    #
    #
    #     if i==1:
    #         AX[i].scatter(1.8, -1.5, marker='*', color='midnightblue', s=100., label='source-B')
    #         AX[i].scatter(3.3, -2.9, marker='*', color='gray', s=100., label='source-0')
    #         AX[i].scatter(2.1, -9, marker='*', color='magenta', s=100., label='source-1')
    #         AX[i].scatter(-2.0, -4.3, marker='*', color='darkolivegreen', s=100., label='source-2')
    #         AX[i].scatter(-7.5, -9.9, marker='*', color='ghostwhite', s=100., label='source-3')
    #     else:
    #         AX[i].scatter(1.8, -1.5, marker='*', color='midnightblue', s=100.)
    #         AX[i].scatter(3.3, -2.9, marker='*', color='gray', s=100.)
    #         AX[i].scatter(2.1, -9, marker='*', color='magenta', s=100.)
    #         AX[i].scatter(-2.0, -4.3, marker='*', color='darkolivegreen', s=100.)
    #         AX[i].scatter(-7.5, -9.9, marker='*', color='ghostwhite', s=100.)
    #     AX[i].tick_params(labelsize=15.)
    #
    # AX[-1].text(3.1, -1.8, 'source-B', fontsize=17.)
    # AX[-1].text(4.9, -3.2, 'source-0', fontsize=17.)
    # AX[-1].text(3.7, -9.3, 'source-1', fontsize=17.)
    # AX[-1].text(-.4, -5., 'source-2', fontsize=17.)
    # AX[-1].text(-5.9, -10.9, 'source-3', fontsize=17.)
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7],aspect=20)
    # cbar=fig.colorbar(img, cax=cbar_ax)
    # cbar.set_label('flux $(10^{-19} \ erg/s/cm^{2}/\AA)$', fontsize=25.)
    # cbar.ax.tick_params(labelsize=15.)
    #
    # plt.show()


    # fig,AX=plt.subplots(2,3,sharey=True)
    # AX=AX.flatten()
    # fig.text(0.08, 0.5, 'Normalized', va='center', rotation='vertical', fontsize=20.)
    # fig.text(0.5, 0.05, 'velocity km/s', ha='center', fontsize=20.)
    # fig.suptitle('Lyman', fontsize=20.)
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
    #                          wavecut=[3990,4060],fitmodel='gaussian',guesses=[0.8,-120,400,0.8,0,400],rewave=1215.673*3.3125,position=[14, 42, 3],
    #                          axis=AX[0],ylabel='Normalized',annotate='mammoth-1')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
    #                          wavecut=[3990, 4060], fitmodel='gaussian', guesses=[0.64, -510, 190, 0.47, 300, 150],
    #                          rewave=1215.673 * 3.3118, position=[15, 37, 3.5],
    #                          axis=AX[1], ylabel='Normalized', annotate='source-0')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
    #                          wavecut=[3990, 4060], fitmodel='gaussian', guesses=[1,-530,270],
    #                          rewave=1215.673 * 3.3118, position=[11, 34, 4],
    #                          axis=AX[2], ylabel='Normalized', annotate='source-1')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
    #                          wavecut=[3990, 4060], fitmodel='gaussian', guesses=[0.81, -800,240,0.55,-41,100],
    #                          rewave=1215.673 * 3.3118, position=[14, 22, 4],
    #                          axis=AX[3], ylabel='Normalized', annotate='source-2')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_ss_icubes.fits',
    #                          wavecut=[3990, 4060], fitmodel='gaussian', guesses=[1,-800,270],
    #                          rewave=1215.673 * 3.3118, position=[7, 19, 6],
    #                          axis=AX[4], ylabel='Normalized', annotate='source-3')
    # fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()

    # fig, AX = plt.subplots(1, 2, sharex=True, sharey=True)
    # AX = AX.flatten()
    # fig.text(0.08, 0.5, 'Normalized', va='center', rotation='vertical', fontsize=20.)
    # fig.text(0.5, 0.05, 'velocity km/s', ha='center', fontsize=20.)
    # fig.suptitle('HeII', fontsize=20.)
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
    #                          wavecut=[5420,5550], fitmodel='gaussian', guesses=[0.8, -120, 400, 0.8, 0, 400],
    #                          rewave=1640 * 3.3405, position=[14, 42, 3.5],
    #                          axis=AX[0], ylabel='Normalized', annotate='mammoth-1')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
    #                          wavecut=[5420,5550], fitmodel='gaussian', guesses=[0.64, -510, 190, 0.47, 300, 150],
    #                          rewave=1640 * 3.339, position=[15, 37, 4.],
    #                          axis=AX[1], ylabel='Normalized', annotate='source-0')
    # fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()

    # fig, AX = plt.subplots(1, 2, sharex=True, sharey=True)
    # AX = AX.flatten()
    # fig.text(0.08, 0.5, 'Normalized', va='center', rotation='vertical', fontsize=20.)
    # fig.text(0.5, 0.05, 'velocity km/s', ha='center', fontsize=20.)
    # fig.suptitle('CIV', fontsize=20.)
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
    #                          wavecut=[5130.,5250.], fitmodel='gaussian', guesses=[0.8, -120, 400, 0.8, 0, 400],
    #                          rewave=1549 * 3.3452, position=[14, 42, 3.5],
    #                          axis=AX[0], ylabel='Normalized', annotate='mammoth-1')
    # SpectralData.Specplotfit(filename='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/1441+4003_comb_psfs_icubes.fits',
    #                          wavecut=[5130.,5250.], fitmodel='gaussian', guesses=[0.64, -510, 190, 0.47, 300, 150],
    #                          rewave=1549 * 3.34301, position=[15, 37, 4.],
    #                          axis=AX[1], ylabel='Normalized', annotate='source-0')
    # fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()

Run()
