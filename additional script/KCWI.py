from cube import CubeData, Cubegenerator
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from mathfunction import statistic
from matplotlib import pyplot as plt
import matplotlib.patches as pth
from astropy import units as u
import numpy as np
import IO
import matplotlib.gridspec as grs
import Cosmology


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

def Specmap(velocube,maskcube,velomap,dec_dis,ra_dis,velocityrange
            ,patchsize=[15,15]):
    '''
    map the extended emission region of IFU data
    and plot the spectra of some region with large
    SNR and significant emission line
    :param velocube: datacube from which we extract spectra
    :param maskcube: maskcube
    :param velomap: flux-weighted velocity map, we use this map to
    determine which region we choose to extract the spectra
    :param dec_dis: dec distance from source-B
    :param ra_dis: ra distance from source-B
    :param velocityrange: velocity range converted from wavelength range
    :param patchsize: size of patch from which we extract the spectra
    :return: None, plot the extracted spectra
    '''

    #create some list to store the information
    #of extracted spectra
    shape=np.shape(velocube)
    declist,ralist,speclist,\
    specflist,noiselist,\
    modellist,paralist=[],[],[],[],[],[],[]

    #divide image into some small patch
    #and traverse all of the patchs, extract
    #spectra from the patch satisified some conditions
    for i in range(0,shape[1],patchsize[0]):
        for j in range(0,shape[2],patchsize[1]):

            #extract the patch from velocity map,
            #find how many pixels in the patch with np.nan
            #if the ratio of the np.nan pixel is more than half of
            #the total number of pixels in the patch, then neglect
            #this patch, else extract the spectra from it and
            #calculate it noise level
            patch=velomap[i:i+patchsize[0],j:j+patchsize[1]]
            if ImgInterSmo.Imgnan(patch):
                spec=SpectralData.RegionSpectrum(velocube[:, i:i + patchsize[0],
                                                 j:j + patchsize[1]])*1e19
                # calculate the noise level, before calculating
                # filter spectra with high pass filter
                noise,specn=SpectralData.Specnoiseestor(spec)

                #if there's more than 30% of pixels with value
                #large than 3 times of the noise level, we think
                #there's significant emission line and extract the
                #spectra to plot,before fit it smooth the spectra
                # to remove the noise with lowpass filter and
                # return the fitted spectra
                spec_filterd = ImgInterSmo.NoiseFilter(spec, 4, 0.2, 'lowpass')
                if ImgInterSmo.Isarraysignificance(spec_filterd,specn):
                    #fit the extracted spectra with gaussian function
                    spec_model,fitpara=statistic.Gaussianfit(velocityrange,spec_filterd,4.5*noise,prominence=.45)
                    paraset=np.array([[fitpara[i],fitpara[i+1],fitpara[i+2]] for i in range(len(fitpara))
                                   if i%3==0 and fitpara[i]<30.])
                    para=paraset[np.where(paraset[:,0]==np.max(paraset[:,0]))][0]
                    print(para)

                    #store the extracted spectra and
                    #fitted model spectra and model parameter
                    specflist.append(spec_filterd)
                    modellist.append(spec_model)
                    speclist.append(spec)
                    declist.append(dec_dis[j])
                    ralist.append(ra_dis[i])
                    noiselist.append(noise)
                    paralist.append(para)
            else:
                continue
    #plot the spectra
    fig1, AX1 = plt.subplots(1)
    img = AX1.pcolor(dec_dis, ra_dis, velomap, cmap='jet')
    cbar = fig1.colorbar(img, ax=AX1, orientation="vertical", aspect=20)
    AX1.set_xlabel('$\Delta$y(kpc)', fontsize=15)
    AX1.set_ylabel('$\Delta$x(kpc)', fontsize=15)
    cbar.set_label('velocity(km/s)', fontsize=25.)
    cbar.ax.tick_params(labelsize=15.)

    fig2,AX=plt.subplots(3,3,sharex=True)
    AX=AX.flatten()
    for i in range(len(AX)):
        ax = AX[i]
        if i<len(speclist):
            v=np.sum(modellist[i]*velocityrange)/np.sum(modellist[i])
            rect = pth.Rectangle((declist[i], ralist[i]), 18, 18,
                                 linewidth=1, edgecolor='black', facecolor='none')
            AX1.add_patch(rect)
            AX1.text(declist[i], ralist[i],str(i))
            AX1.text(declist[i]+7, ralist[i]+7, str(int(paralist[i][2])))
            ax.step(velocityrange,speclist[i])
            ax.plot(velocityrange, modellist[i])
            ax.vlines(0., ymax=np.max(speclist[i]),
                      ymin=np.min(speclist[i]),
                      linestyles='dashed')
            ax.vlines(v,ymax=np.max(speclist[i]),
                      ymin=np.min(speclist[i]),
                      linestyles='dashed',colors='gray')
            ax.hlines(1.*noiselist[i],xmax=np.max(velocityrange),
                      xmin=np.min(velocityrange),
                      linestyles='dashed',colors='red')
            ax.title.set_text('v:'+str(int(v))+'km/s  '+
                              '$\sigma$:'+str(int(paralist[i][2])) + 'km/s')
            ax.text(-1000,0.8*np.max(speclist[i]),str(i))
        else:
            ax.axis('off')
    fig2.text(0.5, 0.05, 'velocity(km/s)', ha='center', fontsize=15.)
    fig2.text(0.05, 0.5, 'Intensity($10^{-19} erg/s/cm^{2}/\AA$)',
              va='center', rotation='vertical', fontsize=15.)
    plt.subplots_adjust(hspace=.35)
    return AX,AX1

def Imgmap(datacube,maskcube,interval_size,wavelengthrange,
           z=2.308,intrinsicwavelength=1215.673):
    '''
    do optimal extraction for the input datacube and generate a series of
    psudo nb image within the velocity range and with fix interval size
    :param datacube: input data cube
    :param maskcube: mask cube
    :param interval_size: velocity size for each psudo nb image, in wavelngth unit
    :param velocityrange: total wavelength range
    :param dec: dec coordiante
    :param ra: ra coordinate
    :return:  a series of psudo nb image array
    '''
    velocityrange = SpectralData.wavelength2velocity(wavelengthrange, z, intrinsicwavelength)/1e3
    cube_shape=np.shape(datacube)
    img_list=[]
    slices_velo_list=[]
    for i in range(0,len(datacube),int(2*interval_size)):
        img_array=CubeData.OptimalextractImg(datacube[i:i+int(2*interval_size),:,:],
                                   maskcube[i:i + int(2 * interval_size), :, :])
        img_array[:, :int(0.15 * cube_shape[0])] = np.nan
        img_array[:, -1 - int(0.07 * cube_shape[0]):] = np.nan
        img_array[:int(0.1 * cube_shape[0]), :] = np.nan
        img_array[-1 - int(0.1 * cube_shape[0]):, :] = np.nan
        slices_velo=np.mean(velocityrange[i:i+int(2*interval_size)])

        slices_velo_list.append(slices_velo)
        img_list.append(img_array)
    return img_list,slices_velo_list

def Indiimg(path=None,filename=None,wavelengthcut=[4020.,4028.],
            wavelengthcut_conti=[4123,4269],internum=None,smoothmark=None):
    '''
    plot the 2D image for a selected wavelength interval
    :param wavelengthcut: wavelength interval
    :return: 2D array of the image and its ra and dec
    '''

    #read the datacube
    flux, wavelength ,RA,DEC= Cubepre(path, filename)
    # flux=CubeData.Cubebadpixelremovor(flux)
    wavelengthcut=np.array(wavelengthcut)

    #generate the 2D image
    twodimg,DEC,RA=ImagPlot.Cutmap(flux.value,wavelength.value,wavelengthcut,DEC.value,RA.value)

    #subtract the continuum component
    if wavelengthcut_conti is not None:
        twodimg_conti, _, _ = ImagPlot.Cutmap(flux.value, wavelength.value, wavelengthcut_conti, DEC, RA)
        twodimg=twodimg-twodimg_conti

    #interpolate and smooth
    if smoothmark is not None:
        twodimg = ImgInterSmo.ImgSmoothor(twodimg, [1.5, 0.428])
    if internum is not None:
        twodimg=ImgInterSmo.Arrayinterpolation(twodimg,internum)


    return twodimg*1e19

def Indiimgreadout(path=None,filename=None,waveinterval=None,wavelengthcut_conti=None,
                   name=None,interpara=None,smoothpara=None):
    '''
    generate the fits file
    :param waveinterval: wavelength interval selected from the datacube
    :param name: name of the output fits file
    :return: None
    '''

    #generate the image of the selected wavelength interval
    twodimg,ra,dec=Indiimg(path=path,filename=filename,wavelengthcut=waveinterval,
                           wavelengthcut_conti=wavelengthcut_conti,internum=interpara,
                           smoothmark=smoothpara)

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

def contouroverlay(baseimgname,contourimgname,xlabel,ylabel,conlabel,baseext=0,contourext=0,
                   contourstd=2e-18,contourlevel=None):
    '''
    overlay contour on the base image
    :param baseimgname: name of base image
    :param baseext: extension number of base image
    :param contourimgname: name of contour image
    :param contourext: extension number of contour image
    :return: img
    '''

    baseimg,baseheader,basewcs=IO.Accessfits(baseimgname,baseext)
    contourimg,contourheader,contourwcs=IO.Accessfits(contourimgname,contourext)
    contourwcs=contourwcs.dropaxis(2)
    norm=ImagPlot.Scaleimgconverter(baseimg,interval='manual')
    ax=plt.subplot(projection=basewcs)
    ax.imshow(baseimg, norm=norm)
    ax.tick_params(labelsize=15.)
    contour = plt.contour(contourimg, transform=ax.get_transform(contourwcs),
                         levels=np.array(contourlevel)*contourstd, linewidths=1.,colors='white')

    lines = [contour.collections[0]];labels=[conlabel]
    plt.legend(lines,labels,fontsize=25)
    plt.xlabel(xlabel,fontsize=20.)
    plt.ylabel(ylabel,fontsize=20.)

    return ax,contour
def velopre(flux,wavelength,ra,dec,waveinterval=[4010.,4050.],z=2.308,
            conti_interval=None,internum=None,mask_sig1=1.,mask_sig2=2.,consub=True,
            coordinate_wcs=True):

    waveinterval = np.array(waveinterval)
    if conti_interval is None:
        conti_interval=waveinterval-(waveinterval[1] - waveinterval[0])
    else:
        conti_interval=np.array(conti_interval)
    # convert ra,dec to angle distance from the mammoth-1 source
    if coordinate_wcs:
        ra_dis, dec_dis = CubeData.Coord2angledis(ra, dec, [0,0],unit=u.deg)
    else:
        ra_dis, dec_dis = CubeData.Coord2angledis(ra, dec, [220.3519792, 40.052625])
        ra_dis,dec_dis=CubeData.Angledis2comovindis(ra_dis,dec_dis,z)
    ra_dis, dec_dis = ra_dis.value, dec_dis.value

    # cut the datacube along the axis of wavelength, only keep the cube within the waveinterval
    rangeflux, rangewavelength = CubeData.CubeCut(flux.value, wavelength.value, 'manual', waveinterval)
    badflux, _ = CubeData.CubeCut(flux.value, wavelength.value, 'manual', waveinterval - 100)
    continumflux, _ = CubeData.CubeCut(flux.value, wavelength.value, 'manual',conti_interval)

    # subtract continuum component, remove bad pixels
    if consub:
        rangeflux_subtracted = CubeData.Continumsubtractor(continumflux, rangeflux)
    else:
        rangeflux_subtracted=rangeflux


    if internum is not None:
        cube_velocity,ra_dis,dec_dis=ImgInterSmo.CubeInterpolation(rangeflux_subtracted,ra_dis,dec_dis,internum)
    cube_velocity = ImgInterSmo.CubeSmooth(cube_velocity, [3., 3.])  # [3., 0.9]
    map = np.mean(cube_velocity, axis=0)


    # do the optimal extraction, we generate two mask cubes for this data cube
    # the first mask is used to remove the influence of the background, because in the process
    # some pixels who have no emission line will be kept and this will influence the velocity map
    # the second mask is used to select the emission region.
    cube1=CubeData.CubeNoiseFilter(cube_velocity,5,.1,'highpass')
    cube_velocity=CubeData.CubeNoiseFilter(cube_velocity, 5, .3)

    mask_cube=CubeData.Maskgenerator2(cube_velocity,cube1,mask_sig2)
    mask_cube1=CubeData.Maskgenerator(cube_velocity,cube1,mask_sig2)
    return cube_velocity,mask_cube,mask_cube1,rangewavelength,ra_dis,dec_dis,cube1

def velocitymap(cube_velocity,mask_cube,rangewavelength,ra_dis,dec_dis,z=2.3093,lamda0=None):
    '''
    calculate the flux-weighted velocity map and plot it
    :param z: redshift use to calculate the observed wavelength lambda0  (lambda-lambda0)*c/lambda0
    :param waveinterval: wavelength interval used to calculate velocity
    :param veloscale: this is the up limit of the velocity map, velocity higher than this value will be replaced with it, -veloscale is the low limit
    :return: None
    '''


    # convert wavelength to velocity
    # mask_cube=np.ones(np.shape(mask_cube))
    velocity = SpectralData.wavelength2velocity(rangewavelength, z, lamda0)
    # convert the flux image to velocity map and velocity dispersion map
    velomap=CubeData.Cubeweightedmean(cube_velocity,velocity,mask_cube[0])
    velodispmap=CubeData.Cubeweightstd(velocity,cube_velocity,mask_cube[0],velomap)/1e3
    velomap = velomap / 1e3

    #do the optimal extraction for data cube
    fluxmap=CubeData.OptimalextractImg(cube_velocity,mask_cube[1])



    #replace the velocity of pixels which have too large velocity
    # velomap[np.where(abs(velomap)>1500.)]=np.nan
    velomap[np.where(velomap<-1600.)]=np.nan
    velomap[np.where(velomap > 1500.)] =np.nan
    velodispmap[np.where(velodispmap > 1000.)] = np.nan
    # CubeData.Imgreadout(velomap,header,'velocitymap.fits')

    return velomap,velodispmap,fluxmap,ra_dis,dec_dis

def velocitymap_filter(cube_velocity,mask_cube,rangewavelength,ra_dis,dec_dis,
                       z=2.3093,lamda0=None,waveinterval=[3990.,4060.],
                       sigma_num=0.25):
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
    velomap,velodisp,fluxmap,ra_dis,dec_dis=velocitymap(cube_velocity,mask_cube,rangewavelength,
                                                        ra_dis,dec_dis,z,lamda0)
    SNmap=CubeData.SNmapgenerator(fluxmap,waveinterval,1e3,600)

    # do the filter, only keep the pixels we need
    velomap_filtered=CubeData.Regionfilter(SNmap,velomap,sigma_num)
    velomap_filtered=velomap
    velo_shape=np.shape(velomap_filtered)


    # if lamda0==1215.673:
    #     #for lyman-alpha
    #     velomap_filtered[:,:int(0.15*velo_shape[0])]=np.nan
    #     velomap_filtered[:,-1-int(0.07*velo_shape[0]):]=np.nan
    #     velomap_filtered[:int(0.1*velo_shape[0]),:]=np.nan
    #     velomap_filtered[-1-int(0.1 * velo_shape[0]):, :] = np.nan
    # elif lamda0==1640:
    #     # for heii
    #     velomap_filtered[:, :int(0.65 * velo_shape[0])] = np.nan
    #     velomap_filtered[:, -1 - int(0.2 * velo_shape[0]):] = np.nan
    #     velomap_filtered[:int(0.3 * velo_shape[0]), :] = np.nan
    #     velomap_filtered[-1 - int(0.15 * velo_shape[0]):, :] = np.nan
    # else:
    #     #for civ
    #     velomap_filtered[:, :int(0.5 * velo_shape[0])] = np.nan
    #     velomap_filtered[:, -1 - int(0.2 * velo_shape[0]):] = np.nan
    #     velomap_filtered[:int(0.33 * velo_shape[0]), :] = np.nan
    #     velomap_filtered[-1 - int(0.2 * velo_shape[0]):, :] = np.nan
    velodisp[np.where(np.isnan(velomap_filtered))]=np.nan
    fluxmap[np.where(np.isnan(velomap_filtered))]=0.





    return velomap_filtered,velodisp,fluxmap,SNmap,ra_dis,dec_dis

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
    ra_dis, dec_dis = CubeData.Coord2angledis(ra, dec, [0,0],unit=u.deg)#220.3522, 40.0522
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




    img=CubeData.OptimalextractImg(cube_velocity,mask_cube)
    SNmap=CubeData.SNmapgenerator(img,waveinterval,1e3,600)
    img_filtered=CubeData.Regionfilter(SNmap,img,5.)

    return img_filtered,ra_dis,dec_dis

# def RegionSpec()
def Sourcemark(AX,sourcecoordinatelist,sourcenamelist,referecepoint=None,z=2.310):

    sourcecoordinatelist=np.array(sourcecoordinatelist)
    if referecepoint is not None:
        sourcecoordinatelist=((((sourcecoordinatelist-referecepoint)*u.deg).to(u.rad)).value)
        sourcecoordinatelist=(sourcecoordinatelist*Cosmology.Angle_distance(z)).to(u.kpc).value
    for i in range(len(AX)):
        for j in range(len(sourcecoordinatelist)):
            AX[i].scatter(sourcecoordinatelist[j][0], sourcecoordinatelist[j][1],
                          marker='+', s=180.,linewidth=2.3, color='black')
            if i == 0:
                AX[i].text(sourcecoordinatelist[j][0] + 2.,
                           sourcecoordinatelist[j][1] + 2.,
                           sourcenamelist[j], color='red',
                           fontsize=15.)
    return None

def Run_sky():
    '''
    check the wavelength calibration of this two data cube
    :return: None
    '''
    # read the cube data, extract the spectra from it and also smooth the spectra
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
    return None

def Run_rawimg():
    '''
    plot raw image of the three emission line
    :return: None
    '''

    fluxss, wavelengthss, rass, decss = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
                                                '1441+4003_comb_ss_icubes.fits')

    fluxpfs, wavelengthpfs, rapfs, decpfs = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
                                                    '1441+4003_comb_psfs_icubes.fits')
    ra,dec=CubeData.Coord2angledis(rass,decss,[220.3519792, 40.052625],u.arcsec)
    ra_dis, dec_dis = CubeData.Angledis2comovindis(ra, dec, 2.310)
    ra_dis, dec_dis=ra_dis.value,dec_dis.value

    interval_list=[[3990,4060],[5459.,5498.],[5160.,5200.]]
    contiinterval_list=[[3800,3850],[5409.,5427.],[5141,5159]]
    img_list,wavecut_list=[],[]
    img_name=['Lyman','HeII','CIV']
    path='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI'
    filename_list=['1441+4003_comb_ss_icubes.fits',
                   '1441+4003_comb_psfs_icubes.fits',
                   '1441+4003_comb_psfs_icubes.fits']
    for i in range(len(filename_list)):
        img=Indiimg(path,filename_list[i],wavelengthcut=interval_list[i],
                    wavelengthcut_conti=contiinterval_list[i],
                    internum=None,smoothmark='smooth')
        img_list.append(img)

    for i in range(len(img_list)):
        fig,ax=plt.subplots(1,1)
        img=ax.pcolor(dec_dis,ra_dis,img_list[i])
        circle = pth.Rectangle((3.5, -10), 10,10,linewidth=1,
                            edgecolor='red', facecolor='none')
        ax.add_patch(circle)
        plt.vlines(2,-130,75,linestyles='dashed',colors='red')
        plt.vlines(15, -130, 75, linestyles='dashed',colors='red')
        ax.text(-150, 60, img_name[i], fontsize=25., color='red')
        cbar = fig.colorbar(img, ax=ax, orientation="vertical", aspect=20)  #
        cbar.set_label('$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$', fontsize=25.)
        cbar.ax.tick_params(labelsize=15.)
        ax.tick_params(labelsize=15.)
        plt.xlabel('$\Delta$x(kpc)',fontsize=20)
        plt.ylabel('$\Delta$y(kpc)',fontsize=20)
        plt.show()



    return None

def Run_spectrum():
    '''
    plot 1d and 2d spectrum
    :return:
    '''

    fluxss, wavelengthss, rass, decss = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
                                                '1441+4003_comb_ss_icubes.fits')

    fluxpfs, wavelengthpfs, rapfs, decpfs = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
                                                    '1441+4003_comb_psfs_icubes.fits')
    ra, dec = CubeData.Coord2angledis(rass, decss, [220.3519792, 40.052625], u.arcsec)
    ra_dis, dec_dis = CubeData.Angledis2comovindis(ra, dec, 2.310)
    ra_dis, dec_dis = ra_dis.value, dec_dis.value
    ra_c,dec_c=CubeData.Findposition(ra_dis,dec_dis,[-7, 0],ra_min=5,dec_min=5)
    interval_list=[[3980,4070],[5439,5518],[5100,5260]]
    flux_list=[fluxss.value,fluxpfs.value,fluxpfs.value]
    wave_list=[wavelengthss.value,wavelengthpfs.value,wavelengthpfs.value]
    name_list=['Lyman','HeII','CIV']
    for i in range(3):
        cube,wave=CubeData.CubeCut(flux_list[i],wave_list[i],'manual',interval_list[i])
        cube_patch=cube[:,ra_c-1:ra_c+2,dec_c:dec_c+1]
        spec=SpectralData.RegionSpectrum(cube_patch)
        noise, specn = SpectralData.Specnoiseestor(spec)
        spec_filterd = ImgInterSmo.NoiseFilter(spec, 2, 0.2, 'lowpass')
        spec_fit,para=statistic.Gaussianfit(wave,spec_filterd,2.5*noise,None)


        cube_slit=cube[:,:,dec_c:dec_c+1]
        tspec,_,spatial_axis=SpectralData.Slitspectrum(cube_slit,None,dec_dis,ra_dis,wave,2.3,4000,False)
        tspec=np.transpose(tspec)[::-1,:]

        fig,AX=plt.subplots(2,1)
        AX=AX.flatten()
        AX[0].step(wave,spec*1e19)
        AX[0].plot(wave, spec_fit * 1e19)
        AX[0].step(wave,specn*1e19,'gray')
        AX[0].set_ylabel('$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$',fontsize=20.)
        img=AX[1].pcolor(wave,spatial_axis[:-1],tspec[:-1,:]*1e19,cmap='jet')
        AX[1].set_ylabel('$\Delta$y(kpc)',fontsize=20.)
        AX[1].set_xlabel('wavelength($\AA$)',fontsize=20.)
        AX[0].tick_params(labelsize=15.)
        AX[1].tick_params(labelsize=15.)
        AX[0].text(interval_list[i][0]+5,np.max(spec*1e19)*.9,name_list[i],color='red',fontsize=25.)
        AX[0].text(interval_list[i][1] - 16, np.max(spec * 1e19) * .9,'$\lambda$:'+str(int(para[1]))+'$\AA$' ,
                   color='black', fontsize=20.)
        cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.35])
        cbar=plt.colorbar(img,ax=AX[1],cax=cbaxes)
        cbar.set_label('$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$', fontsize=20.)
        plt.show()




    return None



def Run_velo():
    SN_sig = 1.
    fluxss, wavelengthss, rass, decss = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
                                                '1441+4003_comb_ss_icubes.fits')
    cube1, maskcube2, maskcube1, wavelength1, ra_dis1, dec_dis1, cuben = velopre(fluxss, wavelengthss, rass, decss,
                                                                                 [3990, 4060],
                                                                                 internum=[2, 8],
                                                                                 mask_sig2=3., consub=True,
                                                                                 coordinate_wcs=False)  # [3990,4060],[4010,4035]

    lymanvelomap, lymandisp, lymanimg, lymanSN, ra_dis, dec_dis = velocitymap_filter(cube1, [maskcube2, maskcube1],
                                                                                     wavelength1, ra_dis1, dec_dis1,
                                                                                     z=2.308,
                                                                                     lamda0=1215.673,
                                                                                     waveinterval=[3990., 4060],
                                                                                     sigma_num=SN_sig)  # .18
    lymanimg = lymanimg * 1e19
    # maskcube2 = np.ones(np.shape(cube1))
    # velocityrange=SpectralData.wavelength2velocity(wavelength1,2.308,1215.673)/1e3
    # cube1,maskcube2=cube1[:,10:-10,10:-10],maskcube2[:,10:-10,10:-10]
    # ra_dis,dec_dis=ra_dis[10:-10],dec_dis[10:-10]
    # lymanvelomap = lymanvelomap[10:-10, 10:-10]
    # lymandisp = lymandisp[10:-10, 10:-10] #* 2.36
    # AX,AX1=Specmap(cube1,maskcube2,lymanvelomap,dec_dis,ra_dis,velocityrange)
    fig = plt.figure(constrained_layout=False)
    gs = grs.GridSpec(1, 3, fig)
    AX = [fig.add_subplot(gs[i]) for i in range(3)]
    imglist = [lymanimg, lymanvelomap, lymandisp]
    cbarlabel = ['$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$', r'$velocity(km/s)$', r'$velocity(km/s)$']
    ImagPlot.Gimgplot(fig, AX, imglist, dec_dis, ra_dis, lymanimg, [2, 4.5, 7.7, 11.8, 17],  '$\Delta$x(kpc)','$\Delta$y(kpc)',
                      cbarlabel, ['Lyman', 'velocity', 'dispersion'], contourmark=[True, None, None])

    # mask_datacube=cube1*maskcube2
    # maskcube3=CubeData.Maskgenerator3(mask_datacube,cuben,.7)
    # img_list,velo_list=Imgmap(cube1*1e19,maskcube3,3,wavelength1)
    #
    # fig,AX=ImagPlot.Mapimgplot(img_list,velo_list,dec_dis1,ra_dis1,'DEC','RA',
    #                         '$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$')
    Sourcemark(AX, sourcecoordinatelist=[[40.05268589,220.3520946],[40.05221366,220.3519858],
                                         [40.05297755,220.3520401],[40.05222735,220.3493004],
                                         [40.0534915,220.3529473],[40.05224151,220.3531833],
                                         [40.04747144,220.3524812],[40.0535024,220.349928],
                                         [40.05106072,220.3497904]],
                   sourcenamelist=['', '', '', '', '',
                                   '','','',''], referecepoint=[40.052625,220.3519792])
    # plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()

    # fluxpfs, wavelengthpfs, rapfs, decpfs = Cubepre('/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI',
    #                                                 '1441+4003_comb_psfs_icubes.fits')
    # SN_sig = 3.8
    # mask_sig1 = SN_sig / 25.
    # cubeheii, maskcubeheii, maskcubeheii2, wavelengthheii, ra_dis, dec_dis ,cuben= velopre(fluxpfs, wavelengthpfs, rapfs, decpfs,
    #                                                                       [5459.,5498.],conti_interval=[5397.,5460.],
    #                                                                       internum=[2, 8], mask_sig1=mask_sig1,
    #                                                                       mask_sig2=.2,coordinate_wcs=False)
    # heiivelomap,heiidisp,heiiimg,heiiSN ,ra_dis, dec_dis=velocitymap_filter(cubeheii,[maskcubeheii,maskcubeheii2],
    #                                                                         wavelengthheii,ra_dis,dec_dis,z=2.340,
    #                                                                         lamda0=1640,waveinterval=[5459.,5498.],
    #                                                                         sigma_num=SN_sig)#.6e-19
    #
    # # maskcubeheii = np.ones(np.shape(cubeheii))
    # # velocityrange = SpectralData.wavelength2velocity(wavelengthheii, 2.340, 1640) / 1e3
    # # cubeheii_cut, maskcubeheii_cut = cubeheii[:, 10:-10, 10:-10], maskcubeheii[:, 10:-10, 10:-10]
    # # ra_dis, dec_dis = ra_dis[10:-10], dec_dis[10:-10]
    # # heiivelomap_cut = heiivelomap[10:-10, 10:-10] * 2.36
    # # Specmap(cubeheii_cut, maskcubeheii_cut, heiivelomap_cut, dec_dis, ra_dis, velocityrange)
    # heiiimg = heiiimg * 1e19
    # fig = plt.figure(constrained_layout=False)
    # gs = grs.GridSpec(1, 3, fig)
    # AX = [fig.add_subplot(gs[i]) for i in range(3)]
    # AX[-1].yaxis.set_ticks_position('right')
    # imglist = [heiiimg,heiivelomap, heiidisp]
    # cbarlabel=['$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$',r'$velocity(km/s)$',r'$velocity(km/s)$']
    # ImagPlot.Gimgplot(fig, AX, imglist, dec_dis, ra_dis, heiiimg, [0.58,.9,1.3,1.8,2],'$\Delta$x(kpc)','$\Delta$y(kpc)',
    #                   cbarlabel, ['HeII','velocity', 'dispersion'],contourmark=[True,None,None])
    # Sourcemark(AX, sourcecoordinatelist=[[40.05268589, 220.3520946], [40.05221366, 220.3519858],
    #                                      [40.05297755, 220.3520401], [40.05222735, 220.3493004],
    #                                      [40.0534915, 220.3529473], [40.05224151, 220.3531833],
    #                                      [40.04747144, 220.3524812], [40.0535024, 220.349928],
    #                                      [40.05106072, 220.3497904]],
    #            sourcenamelist=['', '', '', '', '',
    #                            '', '', '', ''], referecepoint=[40.052625, 220.3519792])
    # plt.show()

    # SN_sig = 3.8
    # mask_sig1 = SN_sig / 8.#8.
    # cubeciv, maskcubeciv, maskcubeciv2, wavelengthciv, ra_dis, dec_dis ,cuben= velopre(fluxpfs, wavelengthpfs, rapfs,
    #                                                                                  decpfs,
    #                                                                                  [5160.,5200.],conti_interval=[5137,5160],
    #                                                                                  internum=[2, 8],
    #                                                                                  mask_sig1=mask_sig1,
    #                                                                                  mask_sig2=.2,coordinate_wcs=False)
    # civvelomap, civdisp,civimg,civSN,ra_dis, dec_dis=velocitymap_filter(cubeciv,[maskcubeciv,maskcubeciv2],
    #                                                                     wavelengthciv,ra_dis,dec_dis,z=2.34385,
    #                                                                     lamda0=1549,waveinterval=[5160.,5200.],
    #                                                                     sigma_num=SN_sig)
    # # maskcubeciv = np.ones(np.shape(cubeciv))
    # # velocityrange = SpectralData.wavelength2velocity(wavelengthciv, 2.34385, 1549) / 1e3
    # # cubeciv_cut, maskcubeciv_cut = cubeciv[:, 10:-10, 10:-10], maskcubeciv[:, 10:-10, 10:-10]
    # # ra_dis, dec_dis = ra_dis[10:-10], dec_dis[10:-10]
    # # civvelomap_cut = civvelomap[10:-10, 10:-10] * 2.36
    # # Specmap(cubeciv_cut, maskcubeciv_cut, civvelomap_cut, dec_dis, ra_dis, velocityrange)
    # civimg=civimg*1e19
    # fig = plt.figure(constrained_layout=False)
    # gs = grs.GridSpec(1, 3, fig, )
    # AX = [fig.add_subplot(gs[i]) for i in range(3)]
    # AX[-1].yaxis.set_ticks_position('right')
    # imglist = [civimg,civvelomap, civdisp]
    # cbarlabel=['$intensity(10^{-19} \ erg/s/cm^{2}/\AA)$',r'$velocity(km/s)$',r'$velocity(km/s)$']
    # ImagPlot.Gimgplot(fig, AX, imglist, dec_dis, ra_dis, civimg, [.603,.823,1.05,1.46,1.74], '$\Delta$x(kpc)', '$\Delta$y(kpc)',
    #                   cbarlabel, ['CIV','velocity', 'dispersion'],contourmark=[True,None,None])
    # Sourcemark(AX, sourcecoordinatelist=[[40.05268589, 220.3520946], [40.05221366, 220.3519858],
    #                                      [40.05297755, 220.3520401], [40.05222735, 220.3493004],
    #                                      [40.0534915, 220.3529473], [40.05224151, 220.3531833],
    #                                      [40.04747144, 220.3524812], [40.0535024, 220.349928],
    #                                      [40.05106072, 220.3497904]],
    #            sourcenamelist=['', '', '', '', '',
    #                            '', '', '', ''], referecepoint=[40.052625, 220.3519792])
    # plt.show()
    return None

def Run_contour():


    # ax, contour = contouroverlay(baseimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/HST_WFC3_MAMMOTH.fits',
    #                              baseext=1,
    #                              contourimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/lyman.fits', contourext=0,
    #                              conlabel='Lyman', ylabel='DEC(J2000)', xlabel='RA(J2000)',
    #                              contourstd=2e-18, contourlevel=[4, 8, 12, 16, 24, 32, 40, 48, 56])
    # ax, contour = contouroverlay(baseimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/HST_WFC3_MAMMOTH.fits',baseext=1,
    #                              contourimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/civ.fits', contourext=0,
    #                              conlabel='CIV', ylabel='DEC(J2000)', xlabel='RA(J2000)',contourstd=1e-18,
    #                              contourlevel=[2,4,6,8,10,12])
    # ax, contour = contouroverlay(baseimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/HST_WFC3_MAMMOTH.fits',baseext=1,
    #                              contourimgname='/Users/shiwuzhang/W&S/ASTRO/MAMMOTH_KCWI/heii.fits', contourext=0,
    #                              conlabel='HeII', ylabel='DEC(J2000)', xlabel='RA(J2000)',contourstd=8e-19,
    #                              contourlevel=[2,4,6,8,10,12])
    # plt.show()
    return None

def Run_edit():

    # Indispectral([220.3491,40.0525])
    # Indiimg(mark='2')


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
    return None

# Run_rawimg()
# Run_spectrum()
Run_velo()