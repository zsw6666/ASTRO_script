import numpy as np
import matplotlib.pyplot as plt



def source(mus,g):
    '''
    This function define the point source
    :param mus: mean value of possion distribution
    :param g: correct factor
    :return: corrected N
    '''


    N=np.random.poisson(lam=mus)
    N=mus+g*(N-mus)

    return N



def PSFC(mu,sigma,N,narray):
    '''
    This is point spread function (standard 2D gussian distribution)
    first generate N random numbers which obey gussian distribution
    and counts the number of points(Xcounts) in each x bins
    then define counts in each x bins as Ny and counts the number of points(Ycounts) in each y bins
    finally save this number in image

    :param coordinate: point's coordinate in image
    :param mu: mean value
    :param sigma: standard deviation
    :param N:total photon number
    :param narray: number of cloumn and row
    :return: 2D standard guassian distribution function
    '''

    mux,muy=mu
    scalex,scaley=sigma
    image=np.zeros((narray,narray)) # initilize the image
    Px=np.random.normal(loc=mux,scale=scalex,size=N)  # N photons obey standard distribution
    Xcounts,Xbins=np.histogram(Px,np.linspace(-5,5,narray+1))  # statistic number of photons in each X bins

    #for n in each x bins,generate n points obey normal distribution and save the value in image varaible
    Ybins=np.linspace(-5,5,narray+1)
    for i in range(len(Xcounts)):
        Py=np.random.normal(loc=muy,scale=scaley,size=Xcounts[i])
        Ycounts,Ybins=np.histogram(Py,Ybins)
        image[:,i]=Ycounts # rewrite the image


    return image



def IMSHOW(sources,loc,scale,g,n):
    '''
    this function generate mean,variance and covariance image
    :param n : number of images
    :param g : correct factor
    :param source : mean value of the possion distribution
    :param loc : mean of the two-dimensional normal distribution
    :param scale standard variation of the two-dimensional distribution
    :return: mean image,variancd image,co-variance image
    '''

    #save the image in imageset and then calculate mean and variance
    imageset=[]
    for i in range(n):
        Nsource=source(sources,g)
        try:
            imagei=PSFC(mu=loc,sigma=scale,N=Nsource,narray=30)
        except ValueError: # there are unexpected error
            print i
            continue
        imageset.append(imagei)

    #calculate mean and variance
    imageset=np.array(imageset)
    mean_image=np.mean(imageset,axis=0)
    var_image=np.var(imageset,axis=0)/10.


    #plot the mean and variance image
    # plt.figure(1)
    # im_mean=plt.imshow(mean_image)
    # plt.colorbar(im_mean)
    # plt.figure(2)
    # im_var=plt.imshow(var_image)
    # plt.colorbar(im_var)



    return mean_image,var_image



if __name__=='__main__':
    mean_mean=[]
    mean_max=[]
    mean_min=[]
    var_mean=[]

    # save the mean max and min value of mean and variance image for different g
    for g in range(0,150,10):
        mean_image,var_image,=IMSHOW(sources=50000, loc=[0, 0], scale=[2, 2], g=g, n=3000)
        mean_mean.append(np.mean(mean_image))
        mean_max.append(np.max(mean_image))
        mean_min.append(np.min(mean_image))
        var_mean.append(np.mean(var_image))


    #plot these images
    plt.figure(1)
    plt.title('g-mean')
    plt.xlabel('g')
    plt.ylabel('mean value')
    plt.scatter(range(1,151,10),mean_mean,color='b',label='mean')
    plt.scatter(range(1, 151, 10), mean_max, color='r',label='max')
    plt.scatter(range(1, 151, 10), mean_min, color='y',label='min')
    plt.legend()
    plt.figure(2)
    plt.title('g-variance')
    plt.xlabel('g')
    plt.ylabel('variance value')
    plt.scatter(range(1, 151, 10), var_mean, color='b')
    plt.show()








