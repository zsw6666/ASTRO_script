from cube import CubeData
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from datareduction import IO
from matplotlib import pyplot as plt
from matplotlib import cm
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
    flux,wavelength,wcs=Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits')

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
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    ra,dec,=CubeData.WCS(wcs)

    #generate the spectra
    onedspectral=SpectralData.SourceSpectral(flux,size,position,ra.value,dec.value)
    plt.figure('spectra',figsize=(17, 7))
    img=SpectralData.SpectralPlot(onedspectral,wavelength)
    lowerline,upperline=SpectralData.WavelengthSelectPlot(onedspectral.value,cutinterval)
    return img,lowerline,upperline

def Indiimg(wavelengthcut=[4020.,4028.]):
    '''
    plot the 2D image for a selected wavelength interval
    :param wavelengthcut: wavelength interval
    :return: 2D array of the image and its ra and dec
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    ra,dec=CubeData.WCS(wcs)

    #generate the 2D image
    twodimg,dec,ra=ImagPlot.Cutmap(flux.value,wavelength.value,wavelengthcut,dec.value,ra.value)
    return twodimg,ra,dec

def Mapimg(internum,cutrange):
    '''
    plot the images for the whole wavelength, devide the wavelength to n intervals and plot image for each interval
    :param internum: interpolate and smooth the image, this variable control how smooth it is
    :param cutrange: the width of each interval
    :return: None
    '''

    #read the datacube
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
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
    flux, wavelength, wcs = Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    ra, dec = CubeData.WCS(wcs)
    map,ummap,lmmap,dec,ra=CubeData.Findabsorption(flux.value,wavelength.value,abrange,upperrange,lowerrange,dec.value,ra.value)
    umimg=ImagPlot.Twodplotimg(ummap,dec,ra)
    lmimg=ImagPlot.Twodplotimg(lmmap,dec,ra)
    img=ImagPlot.Twodplotimg(map,dec,ra)
    plt.show()

def Indiimgreadout(waveinterval,name):
    '''
    generate the fits file
    :param waveinterval: wavelength interval selected from the datacube
    :param name: name of the output fits file
    :return: None
    '''

    #generate the image of the selected wavelength interval
    twodimg,ra,dec=Indiimg(waveinterval)

    #read the header
    header=CubeData.Readheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits')

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


    vscale=np.max(basimg)

    #plot the image
    ax=plt.subplot()
    ax.contour(contouimg,levels=levels,linewidths=.5,color=col)
    plt.imshow(basimg,cmap=cm.jet,vmax=vscale,vmin=-vscale)
    plt.axis('auto')
    plt.axis('off')
    cbar=plt.colorbar()
    cbar.set_label(r'$velocity (km/s) $')

def contouroverlay():
    hstimg,hstheader,hstwcs=IO.Accessfits('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammothicut.fits',0)
    lymanimg,lyheader,lywcs=IO.Accessfits('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/lyman.fits')
    lywcs=lywcs.dropaxis(2)
    coaimg,coaheader,coawcs=IO.Accessfits("/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galA.fits")
    coawcs=coawcs.dropaxis(2)
    cobimg, cobheader, cobwcs = IO.Accessfits(
        "/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galB.fits")
    cobwcs = cobwcs.dropaxis(2)
    cocimg, cocheader, cocwcs = IO.Accessfits(
        "/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galC.fits")
    cocwcs = cocwcs.dropaxis(2)
    codimg, codheader, codwcs = IO.Accessfits(
        "/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/MMAMOTH1/mammoth-1_HST_NB_CO/mammoth_CO_1-0/mammoth.galD.fits")
    codwcs = codwcs.dropaxis(2)

    norm=ImagPlot.Scaleimgconverter(hstimg)
    ax=plt.subplot(projection=hstwcs)
    ax,lycontour= ImagPlot.Contourgenerator(ax, lymanimg, lywcs,[2.25421e-19, 2.68131e-19,
                                                        3.10841e-19, 3.53552e-19,3.96262e-19],'red')#[1.4e-19, 1.8271e-19, 2.25421e-19, 2.68131e-19, 3.10841e-19, 3.53552e-19,3.96262e-19]
    ax,coacontour=ImagPlot.Contourgenerator(ax,coaimg[0],coawcs,[0.0123874,0.0155811,
                                                      0.0187747,0.0219684,0.0251621],'white')#[0.006,0.00919368,0.0123874,0.0155811,0.0187747,0.0219684,0.0251621]
    ax,cobcontour = ImagPlot.Contourgenerator(ax, cobimg[0], cobwcs,
                                   [0.00817988,0.0118598,0.0155396,0.0192195], 'white')#[0.0045,0.00817988,0.0118598,0.0155396,0.0192195]
    ax,coccontour = ImagPlot.Contourgenerator(ax, cocimg[0], cocwcs,
                                   [0.00697652,0.00995303,0.0129295,0.0159061], 'white')#[0.004,0.00697652,0.00995303,0.0129295,0.0159061]
    ax,codcontour = ImagPlot.Contourgenerator(ax, codimg[0], codwcs,
                                   [0.007,0.0100019,0.0130039], 'white')

    ax.imshow(hstimg,norm=norm)
    plt.xlabel(r'$Right Ascension$')
    plt.ylabel(r'$declination$')
    plt.show()
    return None

def velocitymap(z=2.308,waveinterval=[4000.,4050.],veloscale=500.):
    '''
    calculate the flux-weighted velocity map and plot it
    :param z: redshift use to calculate the observed wavelength lambda0  (lambda-lambda0)*c/lambda0
    :param waveinterval: wavelength interval used to calculate velocity
    :param veloscale: this is the up limit of the velocity map, velocity higher than this value will be replaced with it, -veloscale is the low limit
    :return: None
    '''
    flux, wavelength, wcs = Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')#read the datacube, flux is cube
    # header = CubeData.Readheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')

    #convert ra,dec to angle distance from the mammoth-1 source
    ra, dec, = CubeData.WCS(wcs)
    ra_dis,dec_dis=CubeData.Angle2distance(ra,dec)
    ra_dis,dec_dis=ra_dis.value,dec_dis.value

    #cut the datacube along the axis of wavelength, only keep the cube within the waveinterval
    rangeflux,rangewavelength=CubeData.CubeCut(flux.value,wavelength.value,'manual',waveinterval)

    #calculate the emission line within seeing
    rangeflux=CubeData.Cubeseeinglimit(rangeflux)

    #convert wavelength to velocity
    velocity=SpectralData.wavelength2velocity(rangewavelength,z,1215.673)

    #convert the flux image to velocity map
    velomap=CubeData.Cubeweightedmean(rangeflux,velocity)/1e3

    #replace the velocity of pixels which have too large velocity
    # velomap[np.where(velomap>veloscale)]=veloscale
    # velomap[np.where(velomap<-veloscale)]=-veloscale
    # CubeData.Imgreadout(velomap,header,'velocitymap.fits')

    #plot the velocity map
    ax,cbar=ImagPlot.Twodplotimg(velomap,dec_dis,ra_dis,vmax=veloscale,vmin=-veloscale)
    cbar.set_label('$velocity \ km/s)$')
    plt.xlabel(r'arcsec')
    plt.ylabel(r'arcsec')
    # img1=ImagPlot.Twodplotimg(subvelomap,subdec.value,subra.value,vmax=800.,vmin=-800.)
    return None

def slitspectra(position):
    '''
    plot the 2D slit spectra
    :param position: ra,dec of the mammoth-1
    :return: None
    '''

    #read the datacube
    flux, wavelength, wcs = Cubepre('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
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
        pcolor(velocity/1e3,(ra.value-position[0])*3600.,twodspectra.T,cmap=cm.jet)

        # ax.imshow(twodspectra.T,cmap=cm.jet)

    fig.subplots_adjust(wspace=0., hspace=0.)
    fig.text(0.5, 0.04, r'$velocity(km/s)$', ha='center', va='center')
    fig.text(0.06, 0.5, r'$\Delta arcsec$', ha='center', va='center', rotation='vertical')

    plt.show(block=True)

    return None






# Mapspectral(waverange=[4000.,4050.])
# Run_indispectral([25,18])
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
# CubeData.Editheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL1',220.351083333+0.0025202666666643836)
# CubeData.Editheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRVAL2',40.0515833333+0.0017461500000024444)
# CubeData.Editheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRPIX1',14.0)
# CubeData.Editheader('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits','CRPIX2',48.0)
# Indiimgreadout([4000.,4040.],'lyman.fits')
# Indiimgreadout([3800.,4300.],'bw.fits')
# contouroverlay()
velocitymap()
plt.show()
# SingleContour('/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/velocitymap.fits','/Users/shiwuzhang/work&study/ASTRO/MAMMOTH_KCWI/lyman.fits',levels=[1.48181e-19,2.59317e-19,3.70453e-19,4.81588e-19,5.92724e-19],col='red')
# plt.show()
# slitspectra([220.35197916666667,40.052625])