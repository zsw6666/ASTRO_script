import numpy as np
from scipy.linalg import hilbert

def Croutdecomp(arr):
    '''
    use crout fraction to decompose arr to
    upper triangular matrix and lower triangular matrix
    :param arr: inout matrix
    :return: upper triangular matrix and lower triangular matrix
    '''

    #initial upper and lower matrix
    shape=np.shape(arr)
    L,U=np.identity(shape[0]),np.zeros(shape)
    U[0,:]=arr[0,:]
    L[:,0]=arr[:,0]/U[0,0]

    #compute each cell in upper and lower matrix
    for i in range(shape[0]):
        for j in range(shape[0]):
            if i<=j and i>=1:
                frac=np.sum(U[:i,j]*L[i,:i])
                U[i,j]=arr[i,j]-frac
            elif i>j and j>=1:
                frac=np.sum(U[:j,j]*L[i,:j])
                L[i,j]=(arr[i,j]-frac)/U[j,j]

    return L,U

def QRHouseholder(arr):
    '''
    QR decomposition with householder transformation
    :param arr: input array
    :return: Q and R
    '''

    #initial R,v and H, v is used to calculate H
    arr_c=arr.copy()
    R=arr.copy()
    v=arr_c[:,0]
    N=np.shape(v)[0]
    I=np.eye(N)
    H=I.copy()
    i=1


    while np.min(np.shape(arr_c))>=2:

        #calculate H
        v=np.reshape(v,(np.shape(v)[0],1))
        e=np.zeros(np.shape(v))
        e[0]=1.
        v=v-np.linalg.norm(v)*e
        vt=np.reshape(v,(1,np.shape(v)[0]))
        Hi=I-2*np.dot(v,vt)/np.dot(vt,v)
        Hi_shape=np.shape(Hi)
        HI=np.eye(N)
        HI[-Hi_shape[0]:,-Hi_shape[1]:]=Hi

        #multiply H for calculating Q
        H=np.matmul(HI,H)
        #calculate R
        R=np.matmul(HI,R)
        arr_c=R[i:,i:]
        i+=1
        v=arr_c[:,0]
        I=np.eye(np.shape(v)[0])
    Q=H.T
    # R[np.where(np.abs(R)<1e-10)]=0.
    return Q,R

def Msolve(matrix,b,QR):

    if QR:
        m1,m2=QRHouseholder(matrix)
    else:
        m1,m2=Croutdecomp(matrix)

    m1_in,m2_in=np.linalg.inv(m1),np.linalg.inv(m2)
    x=np.dot(m2_in,np.dot(m1_in,b))

    return x

def HMequation(n):
    '''
    intial the function
    :param n: dimension of the matrix
    :return: hm matrix and b
    '''

    hm=hilbert(n)
    shape = np.shape(hm)
    b = np.zeros((1, shape[0]))
    for i in range(shape[0]):
        for j in range(1, n+1):
            b[0, i] = b[0, i] + 1 / (i + j)
    b = np.reshape(b, (n, 1))
    return hm,b

def Run():

    hm, b = HMequation(5)
    x_QR = Msolve(hm, b, True)
    x_LU = Msolve(hm, b, False)
    x_initial = np.repeat(1, np.shape(x_QR)[0])
    print(np.mean(np.abs(x_QR-x_initial)))
    print(np.mean(np.abs(x_LU - x_initial)))
    for i in range(5,20):
        hm, b = HMequation(i)
        x_QR=Msolve(hm,b,True)
        x_LU=Msolve(hm,b,False)
        x_initial=np.repeat(1,np.shape(x_QR)[0])

        #calculate the relative error
        rerr_QR=np.abs(x_QR-x_initial)/np.abs(x_QR)
        rerr_LU = np.abs(x_LU - x_initial) / np.abs(x_LU)

        #tell if there's relative error more than 0.5
        if np.max(rerr_QR)>=0.5:
            print('QR,%d'%i)
            break
        if np.max(rerr_LU) >= 0.5:
            print('LU,%d' % i)
            break


    return None

Run()