'''
This script is used to plot the point spread function of collimator
this is a simple simulation to plot the point spread function of scan imaging
put the source at original point and move the collimator(just change its x,y but keep its z)
then calculate photons arriving at detector
'''

import numpy as np
import matplotlib.pyplot as plt


class PointSource:
    '''
    The point source with a certain power
    '''


    def __init__(self,power,Ephoto,coordinate):
        '''
        initialize the point source with power coordinate and energy of photons
        :param power: power of the source(luminosity)
        :param Ephoto: energy for each photons(assume the photons is monochromatic)
        :param coordinate: coordinate of the source
        '''

        self.p=power
        self.x,self.y,self.z=coordinate
        self.coor=coordinate
        self.Ephoto=Ephoto


    def flux(self,r):
        '''
        This method calculate the number of photons at a certaiu distance per area
        :param r: distance from the source
        :return: number of photons
        '''

        flux=self.p/((4*np.pi)*(r**2))
        nphoto=flux/self.Ephoto

        return nphoto



class Collimator:
    '''
    This is the collimator used for scan imaging with certain property(length,width and hight)
    '''


    def __init__(self,length,width,hight,coordinate):
        '''
        initialize the collimator with length,width hight and coordiante
        :param length: x
        :param width: y
        :param hight: z
        :param coordinate: coordinate of its acme
        '''

        self.deta_x=length
        self.deta_y=width
        self.deta_z=hight
        self.x,self.y,self.z=coordinate


    def collimate(self,source,theta,phi):
        '''
        this mrthod limit photons with a certain angular of incidence can pass it
        :param source: point source
        :param theta: angular of incidence on xz plane
        :param phi: angular of incidence on yz plane
        :return: set of angular which can pass the collimator
        '''

        #the larger angular is in behind of the samll one
        theta.sort()
        phi.sort()

        #calculate the distance and the number of photons collimators can get
        r=np.sqrt(((self.x-source.x)**2)+((self.y-source.y)**2)+((self.z-source.z)**2))
        Nphoto=int(source.flux(r)*self.deta_y*self.deta_x)

        #we assume that if the photons encounter the wall of the collimator,it will be absorbed
        #only photons don't contact with the wall will be detected by detector
        AnguSet=[]
        for photoi in range(Nphoto):
            thetai,phii=np.random.uniform(theta[0],theta[1]),np.random.uniform(phi[0],phi[1])
            if ((self.z-source.z+self.deta_z)*np.tan(thetai)<=self.x+0.5*self.deta_x and (self.z-source.z+self.deta_z)*\
                np.tan(thetai)>=self.x-0.5*self.deta_x) and ((self.z-source.z+self.deta_z)*np.tan(phii)<=self.y+0.5*self.deta_y and (self.z-source.z+self.deta_z)*\
                                             np.tan(phii)>=self.y-0.5*self.deta_y):
                AnguSet.append([thetai,phii])

        #AnguSet include the theta and phi of photons who don't contact with the wall
        #these two angular will be used to calculate the position photons hit the detector
        return AnguSet


    def move(self,x,y):
        '''
        This method is used to move the collimator
        :param x: new position
        :param y: new posioton
        :return: no return
        '''
        self.x=x
        self.y=y



class Detector:
    '''
    This class is detector which detect the photons arrive at the plate
    '''


    def __init__(self,N,z):
        '''
        This method initialize the detector with its size and position
        :param N: size
        :param z: position
        '''
        self.N=N
        self.detectplate=np.zeros((N,N))
        self.z=z


    def detect(self,AnguSet,SourceCoor,ColliCoor):
        '''
        This method calculate the position coming photons hit on detector
        :param AnguSet: coming photons' angular position
        :param SourceCoor: point source's position
        :param ColliCoor:  collimator's position
        :return: no return
        '''

        #tell if there are photons
        if len(AnguSet)!=0:
            #calculate the position photons hit(assume the center of the array is (0,0,z))
            for angu_coor in AnguSet:
                deta_x1,deta_y1=(self.z-SourceCoor[2])*np.tan(angu_coor[0]),(self.z-SourceCoor[2])*np.tan(angu_coor[1])
                deta_x2,deta_y2=(ColliCoor[2]-SourceCoor[2])*np.tan(angu_coor[0])-ColliCoor[0], \
                                (ColliCoor[2] - SourceCoor[2]) * np.tan(angu_coor[1]) - ColliCoor[1]
                pixelx,pixely=int(round(self.N/2+SourceCoor[0]+deta_x2+deta_x1)),int(round(self.N/2+SourceCoor[1]+deta_y2+deta_y1))

                #if this pixel get one photons,the value increase by 1
                self.detectplate[pixelx,pixely]+=1



def ScanImaging(sourcep,ephoton,sizeC,sizeD):
    '''
    this function used to do the scan imaging simulation
    move the collimator and calculate the number of photons hit on detector
    :param sourcep: power of the source
    :param ephoton: energy of each photon
    :param sizeC: size of collimator
    :param sizeD: size of detector
    :return: the image
    '''

    #define source,collimator,detector
    source=PointSource(sourcep,ephoton,np.array([0.,0.,0.]))
    collimator=Collimator(sizeC[0],sizeC[1],sizeC[2],np.array([0.,0.,2.]))
    detector=Detector(sizeD,30)

    #define how the collimator moves
    X=np.linspace(-50,50,300)
    Y=np.linspace(-50,50,300)
    Xv,Yv=np.meshgrid(X,Y)

    #for each position,calculate the photons arrive at the detector
    for coor in zip(Xv.flat,Yv.flat):
        collimator.move(coor[0],coor[1])
        theta=[np.arctan((collimator.x-0.5*collimator.deta_x-source.x)/(collimator.z-source.z)),
                np.arctan((collimator.x+0.5*collimator.deta_x-source.x)/(collimator.z-source.z))]
        phi=[np.arctan((collimator.y-0.5*collimator.deta_y-source.y)/(collimator.z-source.z)),
                np.arctan((collimator.y+0.5*collimator.deta_y-source.x)/(collimator.z-source.z))]
        angset=collimator.collimate(source,theta,phi)
        detector.detect(angset,[source.x,source.y,source.z],[collimator.x,collimator.y,collimator.z])

    #return the image which is the PSF
    return detector.detectplate


if __name__=='__main__':
    PSF=ScanImaging(sourcep=5000,ephoton=0.01,sizeC=[6,6,5],sizeD=220)
    plt.figure(1)
    a=plt.imshow(PSF)
    colbar=plt.colorbar(a)
    colbar.set_ticks(np.linspace(np.min(PSF),np.max(PSF),6))
    plt.show()












