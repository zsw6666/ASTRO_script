import numpy as np
from matplotlib import pyplot as plt
import os


def ReadImage(path,imagename):
    '''
    This function is used to read image
    :param path: path of the image file
    :param imgname: name of image
    :return: image (n*n array)
    '''

    os.chdir(path)
    image=plt.imread(imagename)

    return image


def Vsibility(image):
    '''
    This function used to do fourier transformation
    :param image: original image
    :return: image after fourier transfprmation
    '''
    fft2=np.fft.fft2(image)
    return fft2

def SampleFunction(image_F):
    '''
    This function is used to resample the image(2D guassian function)
    :param image_F: visibility image
    :param u: mean value of guassian function
    :param scale: scale of guassian function
    :return: scatter resample function
    '''
    Xsize,Ysize=np.shape(image_F)
    ux,uy=[256,256]
    SF=np.zeros((Xsize,Ysize))
    scalex,scaley=[1000.,1000.]
    for i in range(Xsize):
        for j in range(Ysize):
            deta=(((i-ux)**2)/scalex)+(((j-uy)**2)/scaley)
            SF[i,j]=(1/(2*scalex*scaley*np.pi))*np.exp((-0.5)*deta)
    SF=SF/np.sum(SF)

    SF=np.fft.ifft2(SF)

    return SF


def MV_image(V_image,S_image):
    '''
    This function used to multiply visibility image with sample image
    :param V_image: visibility image
    :param S_image: sample image
    :return: measurement image
    '''

    image_MV=np.fft.fft2(V_image*S_image)
    return abs(image_MV)


def Dirty(image):
    '''
    This function produce dirty image from visibility image by inverse fourier transform
    :param image: visibility image
    :return: dirty image
    '''

    ifft=np.fft.ifft2(image)
    ifft=np.abs(ifft)
    return ifft

if __name__=='__main__':
    image=ReadImage('/Users/shiwuzhang/ASTRO','lena.jpg')
    image_F=Vsibility(image)
    image_V=abs(image_F)
    image_iSF=SampleFunction(image_F)
    image_SF=abs(image_iSF)
    image_DB=abs(np.fft.fft2(image_iSF))
    image_MV=MV_image(image_F,image_iSF)
    image_D=Dirty(image_MV)
    title=['origin','visibility','sample function','sampled visibility','dirty','dirty beam']
    imageset=[image,image_V,image_SF,image_MV,image_D,image_DB]
    for i in range(6):
        plt.figure(i)
        plt.title(title[i])
        plt.imshow(imageset[i],'gray')
    plt.show()




