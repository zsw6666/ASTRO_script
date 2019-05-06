from cube import CubeData
from spectral import SpectralData
from imag import ImagPlot,ImgInterSmo
from matplotlib import pyplot as plt
from astropy import units as u
import numpy as np
from astropy import constants as const


def Cubepre(path,cubename):
    cube_data, cube_wavelength, wcs, cube = CubeData.ReadCube(path, cubename)
    flux_cube=cube_data*1e-16*u.erg/(u.s*(u.cm**2)*u.AA)
    wavelength=cube_wavelength
    return flux_cube,wavelength,wcs


def Run_Mapspectral():
    flux,wavelength,wcs=Cubepre('/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI','1441+4003_comb_ss_icubes.fits')
    spectral_list,position_list=SpectralData.Mapspectral(flux,[5,2])
    for i in range(len(spectral_list)):
        plt.figure(i,figsize=(17,7))
        img=SpectralData.SpectralPlot(spectral_list[i],wavelength,str(position_list[i][1])+','+str(position_list[i][0]))
        plt.savefig(str(i)+'.png')
        # plt.show()
        # if input()=='a':
        #     continue
    return None

def Run_indispectral():
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    onedspectral=SpectralData.SourceSpectral(flux,[18,3],[5,2])

def Run_indiimg(cut_range):
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    cutflux,cutwavelength=CubeData.CubeCut(flux,wavelength,'manual',cut_range)
    twodimg=(ImagPlot.PlotMap(cutflux)).value
    print(cutwavelength)
    plt.imshow(twodimg)
    plt.show()
    return None
def Run_mapimg():
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    ra,dec=CubeData.WCS(wcs)
    ra,dec=ImgInterSmo.Onedinterpolation(ra.value,4)[1:-1],ImgInterSmo.Onedinterpolation(dec.value,4)[1:-1]
    for i in range(100):
        cutflux, cutwavelength = CubeData.CubeCut(flux, wavelength, 'manual', [16*i,16*(i+1)])
        twodimg = (ImagPlot.PlotMap(cutflux)).value
        twodimg=ImgInterSmo.InSmimg(twodimg,4,[3.,3.])
        print(cutwavelength)
        plt.figure()
        ImagPlot.Colormap(twodimg,dec,ra)
        # plt.show()
        plt.savefig('map'+str(int(np.median(cutwavelength.value)))+'.png')
    return None

def Run_3Dmapimg(cut_range):
    flux, wavelength ,wcs= Cubepre('/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI', '1441+4003_comb_ss_icubes.fits')
    cutflux, cutwavelength = CubeData.CubeCut(flux, wavelength, 'manual', cut_range)
    twodimg = (ImagPlot.PlotMap(cutflux)).value
    ra,dec=CubeData.WCS(wcs)
    ImagPlot.ThreeDmap(twodimg,dec,ra)

# Run_Mapspectral()
# Run_indispectral()
# Run_img([100,170])
# Run_img([0,100])
# Run_indiimg([650,690])
Run_mapimg()
# Run_3Dmapimg([940,990])