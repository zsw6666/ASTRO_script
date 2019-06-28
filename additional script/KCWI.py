from cube import CubeData
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from datareduction import IO
from matplotlib import pyplot as plt
from astropy import units as u
from pylab import pcolor
import numpy as np


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

def Indispectral(position,size=[5,2],cutinterval=[4020.,4028.]):
    '''
    generate the spectra for a certain postion and wavelength interval
    :param position: selected position
    :param size: region size used to generate the spectra, in unit of pixel
    :param cutinterval: wavelength interval
    :return: img of spectra  with two red line parallel to y axis
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits')
    ra,dec,=CubeData.WCS(wcs)

    #generate the spectra
    onedspectral=SpectralData.SourceSpectral(flux,size,position,ra.value,dec.value)
    plt.figure('spectra',figsize=(17, 7))
    img=SpectralData.SpectralPlot(onedspectral,wavelength)
    lowerline,upperline=SpectralData.WavelengthSelectPlot(onedspectral.value,cutinterval)
    return img,lowerline,upperline

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
    if internum is not None:
        twodimg=ImgInterSmo.MapInterpolation(twodimg,internum)

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
    ra,dec=ImgInterSmo.Onedinterpolation(ra.value,internum)[1:-1],ImgInterSmo.Onedinterpolation(dec.value,internum)[1:-1]

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
    contourlyman,ax= ImagPlot.Contourgenerator(ax, lymanimg, lywcs,[2.2e-17,5.1667e-17,8.1334e-17,1.11001e-16],'red')#[1.4e-19, 1.8271e-19, 2.25421e-19, 2.68131e-19, 3.10841e-19, 3.53552e-19,3.96262e-19]
    contourheii,ax = ImagPlot.Contourgenerator(ax, heiiimg, heiiwcs,[5e-18,1.23278e-17,1.96556e-17,2.69833e-17,3.43111e-17],'cyan')
    contourciv,ax = ImagPlot.Contourgenerator(ax, civimg, civwcs, [7e-18,1.0025e-17,1.305e-17,1.6075e-17,1.91e-17], 'lime')
    contourco1,ax=ImagPlot.Contourgenerator(ax,coaimg[0],coawcs,[0.0123874,0.0155811,
                                                      0.0187747,0.0219684,0.0251621],'white')#[0.006,0.00919368,0.0123874,0.0155811,0.0187747,0.0219684,0.0251621]
    contourco2,ax = ImagPlot.Contourgenerator(ax, cobimg[0], cobwcs,
                                   [0.00817988,0.0118598,0.0155396,0.0192195], 'white')#[0.0045,0.00817988,0.0118598,0.0155396,0.0192195]
    contourco3,ax = ImagPlot.Contourgenerator(ax, cocimg[0], cocwcs,
                                   [0.00697652,0.00995303,0.0129295,0.0159061], 'white')#[0.004,0.00697652,0.00995303,0.0129295,0.0159061]
    contourc04,ax = ImagPlot.Contourgenerator(ax, codimg[0], codwcs,
                                   [0.007,0.0100019,0.0130039], 'white')

    ax.imshow(hstimg,norm=norm,cmap='gray')
    lines=[contourlyman.collections[0],contourheii.collections[0],contourciv.collections[0],contourco1.collections[0]]
    labels=['$Lyman\alpha$','$HeII$','$CIV$','CO']
    plt.legend(lines,labels)
    plt.xlabel(r'$Right Ascension$')
    plt.ylabel(r'$declination$')
    plt.show()
    return None

def velocitymap(path=None,filename=None,z=2.3093,lamda0=None,waveinterval=[4010.,4050.],internum=None):
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

    #calculate the emission line within seeing
    # rangeflux=CubeData.Cubeseeinglimit(rangeflux,[3.,1.])

    #convert wavelength to velocity
    velocity=SpectralData.wavelength2velocity(rangewavelength,z,lamda0)

    #subtract continuum component, remove bad pixels
    rangeflux_subtracted=CubeData.Continumsubtractor(continumflux,rangeflux)
    bad_sigma=CubeData.BadvalueEstimator(badflux[:,20:30,4:8])
    rangeflux_badpix = CubeData.Cubebadpixelremovor(rangeflux_subtracted, sigma=bad_sigma)

    #interpolate and smooth cube
    if internum is not None:
        cube_velocity,ra_dis,dec_dis=ImgInterSmo.CubeInterpolation(rangeflux_badpix,ra_dis,dec_dis,internum)
        cube_velocity=ImgInterSmo.CubeSmooth(cube_velocity,[1.5,0.428])


    # convert the flux image to velocity map
    velomap=CubeData.Cubeweightedmean(cube_velocity,velocity)/1e3



    #replace the velocity of pixels which have too large velocity
    # velomap[np.where(velomap>veloscale)]=veloscale
    # velomap[np.where(velomap<-veloscale)]=-veloscale
    # velomap[np.where(abs(velomap) > veloscale)]=0.
    # CubeData.Imgreadout(velomap,header,'velocitymap.fits')

    return velomap,ra_dis,dec_dis



def velocitymap_filter(path=None,filename=None,z=2.3093,lamda0=None,waveinterval=[4000.,4050.],sigma_num=0.25,internum=None):

    velomap,ra_dis,dec_dis=velocitymap(path,filename,z,lamda0,waveinterval,internum)
    fluxmap,_,_=Indiimg(path,filename,waveinterval,internum)
    velomap_filtered=CubeData.Regionfilter(fluxmap,velomap,sigma_num)
    return velomap_filtered,ra_dis,dec_dis


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


def Run():
    # Mapspectral(waverange=[4000.,4050.])
    #  Run_indispectral([25,18])
    # Run_indiimg()
    # Indispectral([220.3491,40.0525])
    # Run_img([100,170])
    # Run_img([0,100])
    # Mapimg(internum=3,cutrange=8)
    # Indiimg(mark='2')
    # Sourcecheck([220.3522,40.0529],[4006.,4014.],mark='2')#source 1
    # Sourcecheck([220.3539,40.0536],[4006.,4014.],mark='2')#source 2
    # Sourcecheck([220.3491,40.0522],[4022.,4030.],mark='2')# source 3
    # Sourcecheck([220.34855,40.0525],[4022.,4030.],mark='2')# source 4
    # Sourcecheck([220.3489,40.0522],[4030.,4038.],mark='2')# source 5
    # Sourcecheck([220.3493,40.0525],[4029.,4038.],mark='2')# source 6
    # Sourcecheck([220.3519792,40.052625],[4000.,4030.],mark='2')
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL1',220.351083333+0.0025202666666643836)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL2',40.0515833333+0.0017461500000024444)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRPIX1',14.0)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRPIX2',48.0)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL1',220.351083333+0.0025202666666643836)
    # CubeData.Editheader('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI','1441+4003_comb_psfs_icubes.fits','CRVAL2',40.0515833333+0.0017461500000024444)
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_ss_icubes.fits',waveinterval=[4010.,4050.],name='lyman.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5464.,5493.],name='heii.fits')
    # Indiimgreadout(path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',filename='1441+4003_comb_psfs_icubes.fits',waveinterval=[5160.,5200.],name='civ.fits')
    # contouroverlay()
    # Sourcecheck([220.3505375,40.05144444],[5460.,5500.],mark='2')
    # Sourcecheck([220.3505375, 40.05144444], [5163., 5203.], mark='2')
    # Sourcecheck([220.3528625, 40.05370556], [5100., 5550.], mark='2')
    lymanvelomap, ra_dis, dec_dis=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits',
                                                     z=2.310,lamda0=1215.673,sigma_num=7.e-18,internum=[1,4])
    # plt.imshow(lymanvelomap,cmap='jet')
    # plt.show()
    lymanimg,_,_=Indiimg('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits',[4000.,4050.],[1,4])
    lymanvelomap[0:10, :], lymanvelomap[60:, :] = np.nan, np.nan
    lymanvelomap[:, 0:8], lymanvelomap[:, 20:] = np.nan, np.nan
    lymanimg[0:2, :], lymanimg[67:, :] = np.nan, np.nan

    heiivelomap, _, _=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',
                                         z=2.340,lamda0=1640,waveinterval=[5464.,5493.],sigma_num=2e-18)
    heiiimg,_,_=Indiimg('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',[5464.,5493.],[1,4])
    heiivelomap[0:5, :], heiivelomap[60:, :] = np.nan, np.nan
    heiivelomap[:, 0:8], heiivelomap[:, 21:] = np.nan, np.nan
    heiiimg[0:2, :], heiiimg[67:, :] = np.nan, np.nan
    #
    #
    civvelomap, _, _=velocitymap_filter('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',
                                        z=2.34385,lamda0=1549,waveinterval=[5160.,5200.],sigma_num=3.2e-18)
    civimg,_,_=Indiimg('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_psfs_icubes.fits',[5160.,5200.],[1,4])
    civvelomap[0:5, :], civvelomap[60:, :] = np.nan, np.nan
    civvelomap[:, 0:8], civvelomap[:, 21:] = np.nan, np.nan
    civimg[0:2, :], civimg[67:, :] = np.nan, np.nan

    ImagPlot.Twodplotimg([lymanvelomap,heiivelomap,civvelomap],dec_dis,ra_dis,subclo=2,subrow=2,xlabel=r'arcsec',
                         ylabel=r'arcsec',cbarlabel='$velocity(km/s))$',
                         subtitle=['$Lyman \alpha$','$HeII$','$CIV$'])
    ImagPlot.Twodplotimg([lymanimg, heiiimg, civimg], dec_dis, ra_dis, subclo=2, subrow=2, xlabel=r'arcsec',
                         ylabel=r'arcsec', cbarlabel='$flux(erg/s/cm^{2})$',
                         subtitle=['$Lyman \alpha$', '$HeII$', '$CIV$'])
    return None

Run()