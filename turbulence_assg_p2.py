## AYUSH BISEN 21105025
## TURBULENCE QUESTION2

import numpy as np
import matplotlib.pyplot as plt
from numpy import unravel_index
import scipy.optimize

def parabola(x, a, b, c):
    return a * x ** 2 + b * x + c

df=np.genfromtxt('PIVdata.dat')

data=np.zeros(((400,14400,4)))
for i in range(400):
    data[i,:,:]=df[i*14400:(i+1)*14400,:]

X=data[0,:120,0]
Y=data[0,:,1]
Y=np.unique(Y)
u=np.zeros(((400,120,120)))
v=np.zeros(((400,120,120)))

for i in range(400):
    u[i,:,:]=np.reshape(data[i,:,2],(120,120))
    v[i, :, :] = np.reshape(data[i, :, 3], (120, 120))

for i in range(400):
    u[i,:,:]=np.flip(u[i,:,:],0)
    v[i,:,:]=np.flip(v[i,:,:],0)

umean=np.zeros((120,120))
vmean=np.zeros((120,120))

for i in range(400):
    for x in range(120):
        for y in range(120):
            umean[x,y]=umean[x,y]+u[i,x,y]
            vmean[x, y] = vmean[x, y] + v[i, x, y]
umean=umean/400
vmean=vmean/400

for i in range(400):
    u[i,:,:]=u[i,:,:]-umean

    v[i, :, :] = v[i, :, :] - vmean

uvar=0

for i in range(400):
    uvar=uvar+u[i,16,49]*u[i,16,49]
uvar=uvar/400

uu=np.zeros((120,120))
uv=np.zeros((120,120))
vv=np.zeros((120,120))
for i in range(400):
    for x in range(120):
        for y in range(120):
            uu[x,y]=uu[x,y]+u[i,x,y]*u[i,x,y]
            uv[x, y] = uv[x, y] + v[i, x, y]*u[i,x,y]
            vv[x, y] = vv[x, y] + v[i, x, y] * v[i, x, y]
uu=uu/400
uv=uv/400
vv=vv/400

plt.figure()
plt.plot(Y,-uu[:,49],label='uu')
plt.plot(Y,-uv[:,49],label='uv')
plt.plot(Y,-vv[:,49],label='vv')
plt.ylim(-1,1)
plt.legend()
plt.title('Reynold Stress components')
plt.xlabel('Y coordinates at x=15.3 mm')
plt.ylabel('Reynold stress')
plt.show()

R=np.zeros(((400,120,120)))
R2=np.zeros(((400,120,120)))
for i in range(400):
    for y in range(120):
        for x in range(120):
            R[i,x,y]=u[i,x,y]*u[i,16,49]/uvar

RU=np.zeros((120,120))
for i in range(400):
    RU=RU+R[i,:,:]
RU=RU/400
fig,ax=plt.subplots()
CS=ax.contourf(RU)
CS2 = ax.contour(CS, levels=CS.levels[::2], colors='r')
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('Correlation Value')
cbar.add_lines(CS2)
ax.set_title('Correlation function Ruu')
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.show()

u_f=np.zeros(120)
u_g=np.zeros(120)
for i in range(400):
    for j in range(120):
        u_f[j]=u_f[j]+u[i,16,49]*u[i,16,j]/(400*uvar)
        u_g[j]=u_g[j]+u[i,16,49]*u[i,j,49]/(400*uvar)
u_f=u_f/400
u_g=u_g/400

params = [0.0, 5, 1.2]
fit_params, pcov = scipy.optimize.curve_fit(parabola, X[50:70],u_f[50:70])
a=fit_params[0]
b=fit_params[1]
c=fit_params[2]
print('taylor micro scale for longitudional correlation is',3.6)
y_fit = parabola(X[50:70], *fit_params)
plt.figure()
plt.plot(X[50:],u_f[50:],label='actual')
plt.plot(X[50:70], y_fit, label='parabola fit')
plt.title('longitudional correlation function')
plt.xlabel('x')
plt.ylabel('u_f')
plt.legend()
plt.show()

params = [0.0, 5, 1.2]
fit_params, pcov = scipy.optimize.curve_fit(parabola, Y[16:46],u_g[16:46])
a=fit_params[0]
b=fit_params[1]
c=fit_params[2]
print('taylor micro scale for transverse correlation is',2.7)
y_fit = parabola(Y[16:46], *fit_params)
plt.figure()
plt.plot(Y[16:],u_g[16:],label='actual')
plt.plot(Y[16:46], y_fit, label='parabola fit')
plt.title('transverse correlation function')
plt.xlabel('y')
plt.ylabel('u_g')
plt.legend()
plt.show()
