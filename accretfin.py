from scipy import integrate
from pylab import *
import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams

from plotparam import *

ax=plt.subplot(111)
plt.subplots_adjust(hspace=0.4)
#ax.text(.55,0.70,'Period Evolution',
       #horizontalalignment='center',transform=ax.transAxes, fontsize=22)

t=linspace(10**0,10**5,500000)
msun=2*10**33
R=1.2*10**6
B=10**15.
G=6.674*10**(-8)
M0=1.4*msun
P0=10**-3
c=3*10**10.
u=B*R**3
eta=1
#initial period

def I(M):
    return .35*M*R**2

def md(t):
    return ((.5*M0*50*t**(-5./3))**(-1)+
    (eta*10**(-3)*.5*M0*t**(.5))**(-1))**(-1)

def rm(M,t):
    return (u**(4./7)*(G*M)**(-1./7)*md(t)**(-2./7))

def rc(f,M,t):
    return (G*M/f**2)**(1./3)

def ndip(f,t):
    return -u**2*f**3/(6*c**3)

def naccret(f,M,t):
    if rm(M,t)>R:
        return (1-f/(G*M/rm(M,t)**3)**.5)*(G*M*rm(M,t))**.5*md(t)
    else:   
        return (1-f/(G*M/R**3)**.5)*(G*M*R)**.5*md(t)
        
def mdot(f,M,t):
    if rm(M,t) > rc(f,M,t): 
        return 0.0
    else:
        return md(t)

def df(f,M,t):
    return (ndip(f,t)+ naccret(f,M,t))/I(M)

def h((f,M),t):
    return (df(f,M,t),mdot(f,M,t))

yinit=(2*pi/P0,M0)
A=integrate.odeint(h,yinit,t)
#############################
def aeriv(z,t):
    return array([((50*t**(-5./3))**(-1)
    +(eta*10**(-3)*t**(.5))**(-1))**(-1)])


zinit=array([0])
z=integrate.odeint(aeriv,zinit,t)

#############################
px=np.linspace(10**(-3),1,500000)
#print A.T[1]/msun
i=0
while A.T[1][i]/msun < 3 and i<len(t)-1:
    i+=1

print t[i]

#print t[ np.abs(A.T[1]-2.5*msun) < 0.1]
plt.subplot(2,1,1)
plt.figure(1)
p1,=loglog(t,2*pi/(A.T[0]),label=r'Mass Accretion', color='k')
plt.xlabel("Time [s]", color='k')
plt.ylabel("Period [s]", color='k')
plt.twinx()
plt.ylabel("Mass Accretion [$M_{\odot}$]", color ='b')
plt.tick_params(labelcolor='b')
p3,=plt.plot(t,(A.T[1])/msun, label=r'Period, \n $\eta=.1$',color='b')
plt.legend((p1,p3),("Period", "Mass Accretion"),loc=(0.7,0.06), frameon=False)
plt.subplot(2,1,2)
plt.xlabel("Time [s]", color='k')
plt.ylabel("Fastness Parameter", color='k')
p2,=loglog(t, (rm(A.T[1],t)/rc(A.T[0],A.T[1],t))**(1.5),color='g')
plt.figure(2)
plt.xlabel('Period [s]')
plt.ylabel('Centri/Grav')
p5,=loglog(2*pi/A.T[0], R*A.T[0]**2/(G*A.T[1]/R**2), color='k')
p1=np.linspace(10**(-3),1,1000)
cg1=(R*4*pi**2/p1**2)/(G*M0/R**2)
plt.xlim(10**-3,10**-1)
p6,=plt.plot(p1,cg1, color='b')
plt.legend((p5,p6),("accretion", "no accretion"),loc=(0.7,0.75), frameon=False)
plt.figure(3)
plt.xlabel('Time [s]')
plt.ylabel('Centri/Grav')
p4,=loglog(t,(R*4*A.T[0]**2)/(G*A.T[1]/R**2))
show()
