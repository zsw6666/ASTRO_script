import numpy as np

__all__=['AxiesInterpolation','MapInterpolation','GussianKernel','GussianFilter','CubeInterpolation','CubeSmooth']
def AxiesInterpolation(map, axies=0):
    '''
    interpolate value to the map
    :param map: original image
    :param axies: axies along which interpolated
    :return: image interpolated along chosen axis
    '''

    # choose which axis to interpolate along
    if axies == 1:
        map = map.T

    # create new array waited values
    interpolated_map = np.zeros((2*np.shape(map)[0]-1, np.shape(map)[1]))
    # apply the mean value of the neighbor pixel
    inter_axies = (map[:-1, :]+map[1:, :])/2.

    # put the original image and interpolated image to the new image
    for axiesi in range(np.shape(interpolated_map)[0]):
        if axiesi % 2 == 0:
            interpolated_map[axiesi, :] = map[int(axiesi/2), :]
        else:
            interpolated_map[axiesi, :] = inter_axies[int((axiesi-1)/2), :]

    # transverse the image
    if axies == 1:
        final_map = interpolated_map.T
    else:
        final_map = interpolated_map

    return final_map


def MapInterpolation(map, internum):
    '''
    interpolate value to the map
    :param map: original map
    :return: interpolated map
    '''

    # interploate original map along x axis and y axis
    for i in range(internum):
        map = AxiesInterpolation(map)
        map = AxiesInterpolation(map, axies=1)
    map = map[::-1]

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

    mux, muy = 0, 0
    scalex, scaley = sigma
    kernel = np.zeros((5, 5))  # initilize the image
    # N photons obey standard distribution
    Px = np.random.normal(loc=mux, scale=scalex, size=10000)
    # statistic number of photons in each X bins
    Xcounts, Xbins = np.histogram(Px, np.linspace(-5, 5, 6))

    # for n in each x bins,generate n points obey normal distribution and save the value in image varaible
    Ybins = np.linspace(-5, 5, 6)
    for i in range(len(Xcounts)):
        Py = np.random.normal(loc=muy, scale=scaley, size=Xcounts[i])
        Ycounts, Ybins = np.histogram(Py, Ybins)
        kernel[:, i] = Ycounts  # rewrite the image
    kernel = kernel/np.sum(kernel)

    return kernel


def GussianFilter(map, kernel):
    '''
    This function use gussian kernel to smooth the image

    :param map: original image
    :param kernel: gussian kernel
    :return: smoothed image
    '''

    # to calculate some pixels on the boundary of image,it first expands the original image
    kernel_map = np.zeros(np.array(np.shape(map))+np.array((4, 4)))
    kernel_map[2:np.shape(kernel_map)[0]-2, 2:np.shape(kernel_map)[1]-2] = map

    # create new map
    map_shape = np.shape(map)
    new_map = np.zeros(map_shape)

    # calculate the value for each pixel
    for x in range(2, map_shape[0]):
        for y in range(2, map_shape[1]):
            # i_map=kernel_map[x-2:x+3,y-2:y+3]
            # pixel_value=kernel*kernel_map[x-2:x+3,y-2:y+3]
            # s=np.sum(kernel*kernel_map[x-2:x+3,y-2:y+3])
            new_map[x-2, y-2] = np.sum(kernel*kernel_map[x-2:x+3, y-2:y+3])

    return new_map[:-2, :-2]


def CubeInterpolation(cube, ra, dec, internum):
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
        cube_inter[i] = MapInterpolation(cube[i], internum)
    ra_inter = MapInterpolation(ra, internum)
    dec_inter = MapInterpolation(dec, internum)
    return cube_inter, ra_inter, dec_inter


def CubeSmooth(cube):
    '''
    smooth the image of every wavelength in this cube
    :param cube: cube waiting smoothing
    :return: smoothed cube
    '''
    kernel = GussianKernel([3., 3.])
    shape0 = np.shape(GussianFilter(cube[0], kernel))
    shape1 = np.shape(cube)
    cube_sm = np.zeros((shape1[0], shape0[0], shape0[1]))
    for i in range(len(cube)):
        cube_sm[i] = GussianFilter(cube[i], kernel)

    return cube