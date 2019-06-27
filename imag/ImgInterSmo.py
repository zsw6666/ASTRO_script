import numpy as np
from scipy import signal, ndimage, interpolate


def cellInterpolation(onedarray):
    former=onedarray[:-1]
    latter=onedarray[1:]
    intercell=(former+latter)/2.
    newarray=np.zeros((1,np.size(intercell)+np.size(onedarray)))[0]
    for i in range(len(newarray)):
        if i%2==0:
            newarray[i]=onedarray[int(i/2)]
        else:
            newarray[i]=intercell[int((i-1)/2)]
    return newarray

def Onedinterpolation(onedarray,internum):
    for i in range(internum):
        onedarray=cellInterpolation(onedarray)
    return onedarray

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

def CubeSmooth(cube,sigma=(3.,)):
    '''
    smooth the image of every wavelength in this cube
    :param cube: cube waiting smoothing
    :return: smoothed cube
    '''
    for i in range(len(cube)):
        cube[i,2:67,:]=ndimage.gaussian_filter(cube[i,2:67,:],sigma)
    return cube

def InSmimg(map,internum,sigma):
    '''
    interpolate and smooth the image
    :param map: map waited to be interpolated and smoothed
    :param internum: paramter used to control the interpolation
    :param sigma: parameter used to control the smooth
    :return: interpolated and smoothed image
    '''

    map=MapInterpolation(map,internum)
    map=ndimage.gaussian_filter(map,sigma)
    return map

def NoiseFilter(data,N,wn):
    '''
    reduce the noise
    :param data: data waited to be rduced
    :param N: filter-control parameter
    :param wn: filter-control parameter
    :return: filtered data
    '''

    b, a = signal.butter(N, wn,'lowpass')
    filtereddata = signal.filtfilt(b, a, data)
    return filtereddata

def Imgseeinglimit(img,size=[3,1]):
    '''
    consider the seeing when access the spectra
    :param img: image used to do the seeing limit(我不知道这里应该叫啥......就是在用到光谱的时候应该考虑一个seeing之内的)
    :param size: size in unit of pixel within a seeing
    :return:
    '''


    img_shape=np.shape(img)

    #for each cell(within a seeing), use the mean intensity as the intensity of each pixel within the cell
    for i in range(0,int(img_shape[0]/size[0])):
        for j in range(0,int(img_shape[1]/size[1])):
            img[int((i-1)*size[0]):int(i*size[0]),int((j-1)*size[1]):int(j*size[1])]=np.mean(img[int((i-1)*size[0]):int(i*size[0]),int((j-1)*size[1]):int(j*size[1])])

    return img
