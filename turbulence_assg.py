## AYUSH BISEN 21105025
## TURBULENCE QUESTION1

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import scipy.optimize


def parabola(x, a, b, c):
    return a * x ** 2 + b * x + c

data = np.genfromtxt('hotwire.dat')
data=data[:24000,:]
t=data[:,0]

u=data[:,1]
umean=np.ones_like(u)*np.mean(u)
u=u-umean #calculating u_fluctuation

uvar=0
for i in range(len(u)):
    uvar=uvar+u[i]*u[i]
uvar=uvar/len(u)

ucorr=[]
for i in range(len(t)):
    sum = 0
    count = 0
    for j in range(len(t)-i):
        sum=sum+u[j]*u[j+i]
        count=count+1
    avg=sum/(count)
    U=avg/uvar
    ucorr.append(U)
ucorr=np.array(ucorr)

integral_timescale=auc(t,ucorr) # calculating area under curve ucorr v/s t
print('Integral time scale is',abs(integral_timescale))

params = [0.0, 5, 1.2]
fit_params, pcov = scipy.optimize.curve_fit(parabola, t[:3], ucorr[:3])
a=fit_params[0]
b=fit_params[1]
c=fit_params[2]
taylor_micro_scale=(-b - (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
print('taylor micro scale is',taylor_micro_scale)

y_fit = parabola(t[:50], *fit_params)
plt.figure()
ax=plt.subplot()
ax.plot(t[:50], y_fit, label='parabola fit')
ax.plot(t[:100],ucorr[:100],label='actual')
plt.ylim(0,1)
plt.legend()
plt.title('CORRELATION FUNCTION R')
plt.xlabel('time(t)')
plt.ylabel('R')
plt.show()




