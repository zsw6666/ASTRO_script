import numpy as np
import matplotlib.pyplot as plt
qr,qz=[0.1,0.1]
pr,pz=[0.1,0.1]
q=np.array([1,(np.exp(1.56)-1)*0.81])
p=np.array([0,0])
q=np.array([0.3,2.448550987])
p=np.array([0.1,0.1])
h=0.05
Q,P=q,p
for i in range(400):
    deta_phi=np.array([(q[0]/((q[0]**2)+(q[1]/0.9)**2))-(0.04/q[0]**3),(q[1]/0.81)/((q[0]**2)+(q[1]/0.9)**2)])
    p=p-deta_phi*h
    q=q+p*h
    Q=np.vstack((Q,q))
    P=np.vstack((P,p))

plt.figure(1)
plt.plot(Q[:,0],Q[:,1])
plt.xlabel('R')
plt.ylabel('z')
plt.title('orbit')
plt.figure(2)

plt.plot(Q[:,0],P[:,0])
plt.xlabel('R')
plt.ylabel('PR')
plt.title('surface of section')
plt.show()

print(Q)
print(P)