import numpy as np
import matplotlib.pyplot as plt
from datareduction import IO

def io(path,name,mark=None):
    '''
    read data
    :param path: where
    :param name: name of data file
    :return: data
    '''
    data=IO.txt_io(path,name,mark)
    return data
def run():
    '''
    for each channel, we should use the flux at this channel to generate the psf and superpose the phtons at each channel then generate the observed spectrum
    '''
    data=io('/Users/shiwuzhang/ASTRO/HW/高能天体物理/极光计划作业/02','crab.txt')
    data2=io('/Users/shiwuzhang/ASTRO/HW/高能天体物理/极光计划作业/01','effective erea.txt')
    energylist=[]
    effecarealist=[]
    for i in data2:
        energy=i[0]*1e-3
        effecarea=i[1]
        energylist.append(energy)
        effecarealist.append(effecarea)
    energylist=np.array(energylist)
    effecarealist=np.array(effecarealist)
    energylist,effecarealist=Interpolation(energylist,effecarealist)
    energy,fwhm,intensity=data[:,0],data[:,1],data[:,2]
    effectiveerea = Fitfunc(energylist, effecarealist, energy)
    effectiveerea[np.where(effectiveerea<0)]=effectiveerea[np.where(effectiveerea<0)[0][-1]+1]
    source_spectrum=Spectrumgenerator(energy,fwhm,intensity)
    plt.plot(energy, source_spectrum, c='b', label='source')
    source_spectrum=source_spectrum*effectiveerea
    bkg_spectrum=Bkggenerator(energy)
    # bkg_spectrum=Spectrumgenerator(energy,fwhm,bkg_intensity)
    spectrum=bkg_spectrum+source_spectrum
    energy_distribution=normalization(spectrum,fwhm)
    photonstime_list=Photonsgenerator(np.sum(spectrum*2*fwhm))
    photonnum_total=np.sum(photonstime_list)
    photondis=energy_distribution*photonnum_total
    new_spectrum=photondis/1800.
    plt.figure(1)
    # plt.plot(energy,intensity,c='y',label='raw')
    # plt.plot(energy,spectrum,c='r',label='source+bkg')
    # plt.plot(energy,source_spectrum,c='b',label='source')
    # plt.plot(energy,bkg_spectrum,c='g',label='bkg')
    plt.title(r'$spectrum$')
    plt.xlabel(r'$energy \ (KeV)$')
    plt.ylabel(r'$Intensity \ (photons \cdot s^{-1} \cdot KeV^{-1})$')
    plt.legend()
    plt.figure(2)
    plt.title(r'$spectrum$')
    plt.xlabel(r'$energy \ (KeV)$')
    plt.ylabel(r'$Intensity \ (photons \cdot s^{-1} \cdot KeV^{-1})$')
    plt.plot(energy,new_spectrum)
    # plt.plot(energy,intensity)
    plt.show()
    a=1
    return None
def Guassianfunc(u,sigma,x):
    '''
    define gaussian function
    :param u: mean value
    :param sigma: standard deviation
    :param x: input value
    :return: output value
    '''
    return np.exp(-0.5*((x-u)/sigma)**2)
def Intergrator(gfunc,u,sigma,interval):
    '''
    define intergrator to calculate the parameter A of gaussian function
    :param gfunc: gaussian function
    :param u: mean value of gaussian function
    :param sigma: standard deviation of gaussian function
    :param interval: domain of defination
    :return: result of intergration in this interval
    '''
    cellarray=np.linspace(interval[0],interval[1],5000)
    deta=cellarray[1]-cellarray[0]
    for i in range(len(cellarray)):
        cellarray[i]=gfunc(u,sigma,cellarray[i])*deta
    intergration=np.sum(cellarray)
    return intergration,deta
def Agenerator(func,u,sigma,interval,area):
    '''
    calculate parameter A of gaussian distribution
    :param func: gaussian function
    :param u: mean value of gaussian function
    :param sigma:  standard deviation of gaussian function
    :param interval: domain of defination
    :param area: area under the gaussian curve, we should keep the number of photons conserved
    :return: parameter A
    '''
    intergration, deta = Intergrator(func, u, sigma, interval)
    A=area/intergration
    return A,deta
def Gaussexpendor(u,sigma,A,disarray):
    '''
    convert deta function to gaussian function
    :param u: mean value of gaussian function
    :param sigma: standard deviation of gaussian function
    :param A: normalization parameter of gaussian function
    :param disarray: bins
    :return: result
    '''
    gaussdis=A*Guassianfunc(u,sigma,disarray)
    return gaussdis
def Sigmagenerator(energy):
    '''
    calculate the standard deviation of gaussian funcrtion
    :param energy: the certain energy
    :return: standard deviation
    '''
    fwhm=0.44*np.sqrt(energy)
    sigma=fwhm/2.355
    return sigma
def Bkggenerator(energy):
    '''
    background
    :param energy: the certain energy
    :return: background counts at this energy
    '''
    bkg=0.15*(energy**-2.6)/(1.4**2.)
    return bkg
def Spectrumgenerator(energy,fwhm,intensity):
    '''
    generate the gaussian extended spectrum
    :param energy: energy channels
    :param fwhm: half width of each energy channel
    :param intensity: intensity for each energy channel
    :return: extended spectrum
    '''
    photondis_list=[]
    flux=2.*fwhm*intensity
    for i in range(len(energy)):
        sigma=Sigmagenerator(energy[i])#calculate sigma for this energy channel
        A,deta=Agenerator(Guassianfunc,energy[i],sigma,[energy[0],energy[-1]],flux[i])#calculate the normalization parameter
        photondis=Gaussexpendor(energy[i],sigma,A,energy)
        photondis_list.append(photondis)#generate the extended spectrum for this energy channel
    photondis_list=np.array(photondis_list)
    spectrum=np.sum(photondis_list,axis=0)
    return spectrum
def Photonsgenerator(lamda):
    '''
    simulate the poisson process of photons arrivation
    :param lamda: mean rate of photons arrivation
    :return: array whose cell represent the photon number in this second
    '''
    photontime_list=np.random.poisson(lamda,1800)
    return photontime_list
def normalization(spectrum,fwhm):
    '''
    normalize the energy spectrum
    :param spectrum: input spectrum
    :param fwhm: half width of each channel
    :return: normalized spectrum
    '''
    flux=spectrum*2*fwhm
    Luminosity=np.sum(flux)
    distribution=spectrum/Luminosity
    return distribution
def Interpolation(x,y):
    k1=(y[1]-y[0])/(x[1]-x[0])
    k2=(y[-1]-y[-2])/(x[-1]-x[-2])
    deta_x=0.08
    for i in range(25):
        x=np.append(x[0]-deta_x,x)
        x=np.append(x,x[-1]+deta_x)
        y=np.append(y[0]-k1*deta_x,y)
        y=np.append(y,y[-1]+k2*deta_x)
    return x,y
def Fitfunc(x,y,z,n=14):
    fitpara=np.polyfit(x,y,n)
    s=np.zeros(np.shape(z))
    for i in range(len(fitpara)):
        s+=fitpara[i]*(z**(n-i))
    return s

run()



# Photonsgenerator(1000)