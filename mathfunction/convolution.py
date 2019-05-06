import numpy as np
import matplotlib.pyplot as plt
def Convolution(func1,func2,high,low):

    x=np.linspace(low,high,500)
    solution=0
    deta=x[1]-x[0]
    for xi in x:
        solution+=func1(xi)*func2(xi)*deta
    return solution

def FuncMulti(func1,func2,high,low):
    x=np.linspace(low,high,500)
    yval=[]
    I=[]
    T=[]
    for xi in x:
        I.append(func1(xi))
        T.append(func2(xi))
        yval.append(func1(xi)*func2(xi))
    T=np.array(T)
    plt.plot(x,T)
    I = np.array(I)
    yval=np.array(yval,dtype='float32')
    # plt.plot(x,yval)
    plt.show()
    return np.array(yval)