import numpy as np
import matplotlib.pyplot as plt

def Gussian(u,x,apha=np.ones(4),scale=1,mean=0):
    '''
    This is the normalized guassian function
    :param u: image's position(numpy array)
    :param x: object's position(numpy array)
    :param apha: parameter list(numpy array)
    :return: value of gussian function
    :param scale: standard deviation
    :param mean: mean value
    :return: value of gussian function
    '''

    U=np.mat(np.array([1,u,1,u])*apha)
    X=np.mat(np.vstack([x,x,x**2,x**2]))
    variable=np.array(u-x+U*X)#calculate variable in PSF
    gussian=(1/(scale*np.sqrt(2*np.pi)))*np.exp(-((variable-mean)**2)/(2*(scale**2)))

    return gussian


def Object(x):
    '''
    This is the object function
    :param x: position(numpy array)
    :return: value at each position
    '''
    #return x**2
    return np.sin(x)
    #return np.ones(len(x))


def plot(x,y,c='r',title=None,xlabel=None,ylabel=None):
    '''
    This function used to plot k(u,x),o(x),i(u)
    :param x: independent variable
    :param y: dependent variable
    :param c: color
    :param title: title of the image
    :param xlabel: xlabel
    :param ylabel: ylabel
    :return:image
    '''
    plt.plot(x,y,c)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def Image(kernel,object,x,u,apha,scale=1):
    '''
    This function generate the image
    :param kernel: kernel function
    :param object: object function
    :param x: object position
    :param u: image position
    :param apha: parameter list
    :param scale: scale of the kernel function
    :return: image value
    '''

    objvalue=object(x)
    ivalue=np.ones(len(u))#init i(u) the default value of i(u) is 1
    for i in range(len(u)):
        ivalue[i]=np.sum(objvalue*kernel(u[i],x,apha,scale=scale))#calculate image value at u pisition

    return ivalue

if __name__=='__main__':
    u,x=np.linspace(-20,20,200),np.linspace(-10,10,200)
    apha=np.array([1e-2,2e-2,3e-2,4e-2])#define the parameter list
    ivalue=Image(Gussian,Object,x,u,apha,scale=2)#ca
    plot(u,ivalue,title='image',xlabel='u',ylabel='image value')
    plot(x,Object(x),title='object',xlabel='x',ylabel='object value')
    imagez=np.ones(len(u))
    for i in range(len(u)):
        imagez[i]=Gussian(u[i],np.array([0]),apha,scale=10)
    plot(u,imagez,title='kernel',xlabel='u',ylabel='kernel value')



