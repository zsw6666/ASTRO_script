import numpy as np
import matplotlib.pyplot as plt


def PartialDerivative(f,x0,y0,step=1e-5):
    '''
    this function is used to derivative
    :param f: function detivatived
    :param x0: varaible x
    :param y0: varaibel y
    :param step: small step used to calculate the derivative
    :return: value of derivative at this location
    '''

    derivative=np.array([(f(x0+step,y0)-f(x0,y0))/step,(f(x0,y0+step)-f(x0,y0))/step])

    return derivative

def Potential(R,z):
    '''
    potential of this grivitational field
    :param R: coordiante r
    :param z: coordinate z
    :return: potential at a certain location
    '''
    v0=1.
    q=.7
    Lz=0.1
    PHI=0.5*(v0**2)*np.log((R**2)+((z/q)**2))+0.5*(Lz/R**2)

    return PHI

def maxVR(E,R):
    P=Potential(R,0.)
    return np.sqrt(2.*(E-P))

def RungeKutta(w,step):
    '''
    Runge-Kutta method used to calculate the new location of point in phase-space
    :param w:
    :param step:
    :return: new location in phase-space
    '''
    k1=step*np.array([w[1],-PartialDerivative(Potential,w[0][0],w[0][1])])
    k2=step*np.array([w[1]+0.5*k1[1],-PartialDerivative(Potential,w[0][0]+0.5*k1[0][0],w[0][1]+0.5*k1[0][1])])
    k3=step*np.array([w[1]+0.5*k2[1],-PartialDerivative(Potential,w[0][0]+0.5*k2[0][0],w[0][1]+0.5*k2[0][1])])
    k4 = step * np.array([w[1] + k3[1], -PartialDerivative(Potential,w[0][0] + k3[0][0], w[0][1] + k3[0][1])])
    w=w+(1/6.)*(k1+2.*k2+2.*k3+k4)

    return w

def Iteration(R,vR,E,time):
    '''
    generate the orbit and surface of section by Runge-kutta method
    :param q: spatial coordiante
    :param p:momentum coordinate
    :param time: number of steps
    :return: orbit in phase space
    '''
    h=0.05
    z=0.
    vz=np.sqrt(2*(E-Potential(R,z)-(vR**2)/2.))
    q=np.array([R,z])
    p=np.array([vR,vz])
    Q, P = q, p
    for i in range(time):
        w=np.array([q,p])
        w=RungeKutta(w,h)
        q,p=w

        Q = np.vstack((Q, q))
        P = np.vstack((P, p))

    return Q,P

def SurfaceoSection(Q,P):
    index=np.array([])
    for i in range(len(Q[:,1])-1):
        if Q[i,1]<0 and Q[i+1,1]>0:
            index=np.append(index,i+1)

    indexR=[]
    for i in index:
        indexR.append(int(i))
    R=Q[indexR,0]
    vR=P[indexR,0]
    Z=Q[indexR,0]
    return R,vR,Z

def plot():
    E=1.25
    R=0.8
    z=0.
    VR=np.linspace(0.,0.9,3)
    i=0

    plt.figure(1)
    for vR in VR:
        Q, P = Iteration(R, vR, E, 15000)
        r, vr, z = SurfaceoSection(Q, P)
        # plt.figure(i)
        plt.scatter(r,vr,s=1.,c='b')
        i+=1
    plt.title(r'surface of section')
    plt.xlabel(r'$R\,(\mathrm{kpc})$')
    plt.ylabel(r'$v_R\,(\mathrm{km\,s}^{-1})$')
    plt.show()
    return None

# plot()
a=Iteration(0.8,0.9,1.25,15000)