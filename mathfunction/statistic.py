import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from matplotlib import pyplot as plt

def Nd2one1(ndarray):

    onedarray=ndarray.flatten()
    return onedarray

def Distributiongenerator(sample,binnum='auto'):

    counts,bins=np.histogram(sample,bins=binnum)
    return counts,bins[1:]

def Gaussian(x,A,mu,sig):
    '''
    Gaussian function
    :param x: independent variable
    :param A: amplitude
    :param mu: center of gaussian function
    :param sig: sigma
    :return: array of dependent variable of x
    '''
    return A*np.exp(-.5*((x-mu)**2)/sig**2)

def Multigaussian(n=1):
    '''
    dynamicall multi-gaussian function
    :param n: number of gaussian function convovled
    :return: function improved gaussian function
    '''
    def wrapper(x,*para):
        mgaussian=0.
        for i in range(0,3*n,3):
            mgaussian=mgaussian+Gaussian(x,para[i],para[i+1],para[i+2])
        return mgaussian
    return wrapper

def Curvefit(fitfunc,x,y,initpara=None):
    '''
    fit curve with scipy.optimize.curve_fit
    :param fitfunc: function uesd to fit the data
    :param x: independent variable array
    :param y: dependent variable array(fitting data)
    :param initpara: initial parameter input to fit function
    :return: fit parameter
    '''
    popt, pcov = curve_fit(fitfunc, x, y, p0=initpara,method='lm',maxfev=60000)
    return x,popt

def Peakfinder(x,height=None):
    '''
    find the peak of the data above a certain noise level
    :param x: data
    :param height: noise level beyond which we find peak
    :return: indecis of the peak in the data
    '''
    return find_peaks(x,height)[0]

def GaussianParaestimator(x,y,noise):
    '''
    estimate the initial parameter of the gaussian function
    :param x: independent variable
    :param y: dependent variable
    :param noise: noise
    :return: parameter array
    '''
    #using Peakfinder to find the peak and
    #using the value of the peak as the
    #amplitude of gaussian function
    #using its x position as the center
    #the default sigma use the value of 100 km/s
    peak_indecis = Peakfinder(y, height=4 * noise)
    para=[[y[i],x[i],100] for i in peak_indecis]
    return np.array(para).flatten()

def Gaussianfit(x,y_init,noise):
    '''
    fit data with multi-gaussian function
    :param x: independent variable
    :param y_init: dependent variable
    :param noise: noise
    :return: model data and model paramter
    '''
    y=y_init.copy()

    #estimate the initial parameter
    initpara = GaussianParaestimator(x, y, noise)

    #remove the noise
    y[np.where(y<3*noise)]=.2*noise
    y[-5:]=.2*noise
    #fit the data and apply the fit paramter to function
    #to calculate the model data
    _, popt = Curvefit(Multigaussian(int(len(initpara) / 3)), x,
                                 y, initpara)
    y_model = Multigaussian(int(len(initpara) / 3))(x, *popt)
    return y_model,popt

def test():
    x=np.linspace(-10,10)
    y=Multigaussian(2)(x,10,0,2,10,5,3)
    y_noise=y+np.random.normal(0,3,len(y))
    x,popt=Curvefit(Multigaussian(2),x,y_noise)
    plt.plot(x,y_noise)
    plt.plot(x,Multigaussian(2)(x,*popt))
    plt.show()


# test()
