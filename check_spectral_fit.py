import numpy as np
import matplotlib.pyplot as plt
from cube import CubeData
from imag import ImagPlot, ImgInterSmo
from photutils import DAOStarFinder as sourcefinder


path = '/Users/shiwuzhang/ASTRO/MAMMOTH_KCWI'
cube_name = '1441+4003_comb_ss_icubes.fits'
# cube_name='1441+4003_comb_psfs_icubes.fits'


def FindSource(map, FWHM=3., sigma=3.):
    '''
    this function is used to find the point-like source
    :param map: image's 2-dimensional array
    :param sigma: limit beyond which we consider as the source
    :return: source's coordinate in image
    '''
    # convert 2D array to 1D for convenience
    photon = map.flatten()

    # calculate the standard deviation and mean value of pixels
    photon_std = float(np.std(photon))
    photon_mean = np.mean(photon)
    boundary = photon_mean+sigma*photon_std

    # use the module already exist
    daofind = sourcefinder(fwhm=FWHM, threshold=sigma*photon_std)
    source = daofind(map-photon_mean)

    y, x = [int(i) for i in np.round(source['xcentroid'].data)], [
        int(i) for i in np.round(source['ycentroid'].data)]
    coordinate_map = np.vstack((x, y)).T

    return coordinate_map, boundary


def LocateTarget(cube, target_coordinate):

    map = ImagPlot.PlotMap(cube._data)
    wavelength, dec, ra = cube.world[:]
    deta_dec, deta_ra = np.abs(
        dec.value[0, :, :]-target_coordinate[1]), np.abs(ra.value[0, :, :]-target_coordinate[0])

    x, y = np.where(deta_ra == np.min(deta_ra))[
        0], np.where(deta_dec == np.min(deta_dec))[1]

    map[x, y] = 10

    plt.imshow(map)
    plt.show()



gain = 0.145
cube_data, cube_wavelength, dec, ra, cube = CubeData.ReadCube(path, cube_name)
flux_all, wavelength_all, ra_cut, dec_Cut = CubeData.FLux(
    cube=cube_data, wavelength_axis=cube_wavelength, gain=gain)

# flux_all_sub,continuum_std=CubeData.ContinuumSub(flux_all,flux_all)
# flux_all_sub_inter=ImgInterSmo.CubeInterpolation(flux_all_sub,2)
# flux_all_sub_inter_sm=ImgInterSmo.CubeSmooth(flux_all_sub_inter)
# flux_all_sub_inter_sm_map=ImagPlot.PlotMap(flux_all_sub_inter_sm)
# ImagPlot.Colormap(flux_all_sub_inter_sm_map)
flux_CV, wavelength_CV, ra_cut, dec_cut = CubeData.FLux(
    cube=cube_data, wavelength_axis=cube_wavelength, ra=ra, dec=dec, gain=gain, emissionline='lyman')
flux_CV_sub, continuum_std = CubeData.ContinuumSub(flux_CV, flux_all)
# map=PlotMap(flux_CV_sub)
# plt.imshow(map)
# plt.show()
# coordinate_source=np.array([[32,11],[22,15],[39,14]])
# Spectral(flux_CV_sub,coordinate_source,wavelength_CV)
flux_CV_sub_inter, ra_inter, dec_inter = ImgInterSmo.CubeInterpolation(
    flux_CV_sub, ra_cut, dec_cut, 4)
flux_CV_sub_inter_sm = ImgInterSmo.CubeSmooth(flux_CV_sub_inter)
flux_CV_sub_inter_sm_map = ImagPlot.PlotMap(flux_CV_sub_inter_sm)
coordinate_source, boundary = FindSource(
    flux_CV_sub_inter_sm_map, FWHM=150., sigma=7.)
for coordinate in coordinate_source:
    flux_CV_sub_inter_sm_map[coordinate[0], coordinate[1]] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0]-1, coordinate[1]-1] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0]+1, coordinate[1]+1] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0], coordinate[1]-1] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0]-1, coordinate[1]] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0], coordinate[1]+1] = 1e-17
    flux_CV_sub_inter_sm_map[coordinate[0]+1, coordinate[1]] = 1e-17
ImagPlot.Colormap((flux_CV_sub_inter_sm_map.T)[::-1, :], dec=[np.mean(ra_inter[-1, :]), np.mean(
    ra_inter[0, :])], ra=[np.mean(dec_inter[:, 0]), np.mean(dec_inter[:, -1])])

# coordinate_map,boundary=FindSource(smoothed_map,FWHM=4.6,sigma=1.2)#4.6 and 1.2 are the best value to select the three source
# coordinate_map,boundary=FindSource(flux_map,FWHM=3.,sigma=5.)
