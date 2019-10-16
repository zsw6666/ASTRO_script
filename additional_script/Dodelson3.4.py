import numpy as np
import astropy.constants as const
import astropy.units as u
from matplotlib import  pyplot as plt
def lamda(x):
    s=(255/(886.7*(x**5)))*(12.+6.*x+x**2)
    return s
def H(x):
    Q=((1.293*u.MeV).to(u.J)).value
    g_star=10.75
    s=(2*np.pi/3)*np.sqrt((g_star*np.pi*const.G.value)/(5.))*((Q**2)/(x**2))
    return s
def Mu(x):
    tao=886.7
    Q=((1.293*u.MeV).to(u.J)).value
    s1=(-255./(tao*Q))*(4*(np.pi**3)*const.G.value*(Q**2)*10.75/45.)**(-0.5)
    s2=(4/x**3)+(3/x**2)+(1/x)+np.exp(-x)*((4/x**3)+(3/x**2))
    return s1*s2
def Intfunc(x,x0):
    s1=(lamda(x)*np.exp(-x))/(x*H(x))
    s2=np.exp(Mu(x)-Mu(x0))
    return s1*s2
def Integrator(func,interval,x0):
    x=np.linspace(interval[0],interval[1],1000)
    deta=x[1]-x[0]
    y=Intfunc(x,x0)
    s=np.sum(y*deta)
    return s
def Kai_n(x):
    s=Integrator(Intfunc,[0.1,x],x)
    return s

def Run():
    xlist=np.linspace(1,1000,1000)
    ylist=[]
    for x in xlist:
        ylist.append(Kai_n(x))
    ylist=np.array(ylist)
    plt.plot(xlist,ylist)
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$X_{n}$')
    plt.show()
    return None

Run()
