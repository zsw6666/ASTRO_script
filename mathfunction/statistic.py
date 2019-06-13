import numpy as np

def Nd2one1(ndarray):

    onedarray=ndarray.flatten()
    return onedarray

def Distributiongenerator(sample,binnum='auto'):

    counts,bins=np.histogram(sample,bins=binnum)
    return counts,bins[1:]

def Sigmaclip(sample,n=4.):
    for i in range(30):
        sample_sigma=np.std(sample)
        sample_expectation=np.median(sample)
        sample[np.where(abs(sample-sample_expectation)>n*sample_sigma)]=np.median(sample)
    return sample