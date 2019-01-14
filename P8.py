import numpy as np
import matplotlib.pyplot as plt


def Signal(origin=100,m=1):
    '''
    This function produce the real signal
    :param origin: original value of this signal
    :param m: index
    :return: value of signal with index of m
    '''

    return origin*0.6**m


def Noise(mean,sigma):
    '''
    This function produce noise
    :param mean: mean value of the noise
    :param sigma: scale of the noise
    :return: value of noise
    '''

    return np.random.normal(mean,sigma,1)

def Filter(w):
    '''
    This is the filter which remove noise from the real signal
    :param w: real signal(with noise)
    :return: signal after filtered
    '''

    x=np.array([])
    #value of parameters of filter
    a=np.array([0.451,0.165])

    for n in range(len(w)-1):
        x=np.append(x,np.sum(a*w[n:n+2]))

    return x


def Experiment():
    '''
    firstly input the ideal signal
    then add it with noise to produce the real signal
    next filter the real signal
    finally plot them
    '''

    sig,noi=np.array([]),np.array([])

    #generate the sequence of noise and ideal signal
    for i in range(200):
        sig=np.append(sig,Signal(m=i))
        noi=np.append(noi,Noise(mean=0,sigma=5.))

    w=sig+noi

    sig_f=Filter(w)

    plt.figure(1)
    plt.plot(range(200),sig,color='b',label='ideal signal')
    plt.plot(range(200),w,color='g',label='real signal')
    plt.plot(range(199),sig_f,color='r',label='signal after filter')
    plt.legend(loc=0,)
    plt.show()


if __name__=='__main__':

    Experiment()



