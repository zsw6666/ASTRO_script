import numpy as np
from scipy import signal, ndimage


def Arrayinterpolation(array, internum):
    '''
    interpolation the nD numpy array
    :param array: numpy array waited for interpolation
    :param internum: control the interpolation
    :return: interpolated array
    '''
    array_intered=ndimage.zoom(array,internum,order=1)
    return array_intered

def ImgSmoothor(img,sigma):
    '''
    smooth the 2D image
    :param img: 2D numpy array
    :param sigma: parameter controling the smooth
    :return: smoothed image
    '''
    img_smooth=ndimage.gaussian_filter(img,sigma)
    return img_smooth

def CubeInterpolation(cube, ra=None, dec=None, internum=None):
    '''
    interpolate the image of every wavelength in the cube
    :param cube: cube waiting interpolation
    :param internum: number of interpolation
    :return: interpolated cube
    '''
    # create the new array used to carry the new interpolated array
    cube_shape=np.shape(cube)
    inter_shape=(cube_shape[0],cube_shape[1]*internum[0],cube_shape[2]*internum[1])
    cube_inter=np.zeros(inter_shape)

    if (ra is not None) and (dec is not None):
        # interpolate the coordinate
        ra=Arrayinterpolation(ra,internum[0])
        dec=Arrayinterpolation(dec,internum[1])


    # for each slice in cube, interpolate the image
    for i in range(np.shape(cube)[0]):
        cube_inter[i,:,:]=Arrayinterpolation(cube[i,:,:],internum)

    return cube_inter,ra,dec

def CubeSmooth(cube,sigma=(3.,)):
    '''
    smooth the image of every wavelength in this cube
    :param cube: cube waiting smoothing
    :return: smoothed cube
    '''

    # for each slice in cube, smooth the image
    for i in range(len(cube)):
        cube[i,:,:]=ImgSmoothor(cube[i,:,:],sigma)
    return cube

def InSmimg(map,internum,sigma):
    '''
    interpolate and smooth the image
    :param map: map waited to be interpolated and smoothed
    :param internum: paramter used to control the interpolation
    :param sigma: parameter used to control the smooth
    :return: interpolated and smoothed image
    '''

    map=Arrayinterpolation(map,internum)
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

# if __name__=='__main__':
#     import matplotlib.pyplot as plt
#     x = np.linspace(0, 10*np.pi, 20)
#     y = np.cos(x)
#     yinter=Arrayinterpolation(y,1000)
#     xinter=np.linspace(0,10*np.pi,1000)
#     plt.plot(xinter,yinter)
#     plt.plot(x,y)
#     plt.show()

