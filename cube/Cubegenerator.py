import numpy as np
import  matplotlib.pyplot as plt


'''
This script is used to generate the test data cube
author: shiwu
date: 2019.07.12 
'''


def Cubegenerator(shape,wavelengthrange):

    rawcube=np.zeros(shape)
    spectral_length=shape[0]
    row_length=shape[1]
    for i in range(row_length):
        central_wavelength=(i*(wavelengthrange[1]-wavelengthrange[0])/row_length)+wavelengthrange[0]
        width=0.2*(wavelengthrange[1]-wavelengthrange[0])+np.random.normal(0,0.02*(wavelengthrange[1]-wavelengthrange[0]),1)
        rawcube[:,i,:]=rawcube[:,i,:]+(Spectragenerator(spectral_length,wavelengthrange,central_wavelength,width)).reshape(-1,1)

    return rawcube

def Spectragenerator(spectra_length,wavelengthrange,central_wavelength,width):

    wavelength=np.linspace(wavelengthrange[0],wavelengthrange[1],spectra_length)
    spectra_pure=gaussian(wavelength,central_wavelength,width)
    noise=np.random.normal(0.,0.1*np.max(spectra_pure),spectra_length)
    spectra=spectra_pure+noise
    return spectra


def gaussian(x, mu, sig,A=1):
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


if __name__ =='__main__':
    datacube=Cubegenerator((100,300,100),[4000,4060])
    map=np.sum(datacube,axis=0)
    plt.imshow(map,cmap='jet')
    plt.show()