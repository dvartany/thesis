from scipy import integrate
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from pylab import *
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='palatino')
rc('font', weight='medium')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=2)
rc("lines", linewidth=3)
rc('axes', labelsize=18) #24
rc("axes", linewidth=2) #2)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)
rc('legend', fontsize=13) #16
rc('xtick.major', pad=8) #8)
rc('ytick.major', pad=8) #8)
rc('xtick.major', size=13)
rc('ytick.major', size=13)
rc('xtick.minor', size=7)
rc('ytick.minor', size=7)
rcParams['text.latex.preamble']=[r"\usepackage{color}"]
rcParams['text.latex.preamble']=[r"\usepackage{xcolor}"]
rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

ax=plt.subplot(111)
plt.subplots_adjust(hspace=0.4)
ax.text(.5,0.90,'Period Evolution, $\eta=.1$',
        horizontalalignment='center',transform=ax.transAxes, fontsize=22)

t=linspace(10**0,10**5,500000)
R=1.2*10**6
B=10**15.
G=6.674*10**(-8)
M0=2*1.4*10**33.
c=3*10**10.
u=B*R**3
r=8.13*10**14
eta=10

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
        return 0
    else:
        return md(t)

def df(f,M,t):
    return (ndip(f,t)+ naccret(f,M,t))/I(M)

def h((f,M),t):
    return (df(f,M,t),mdot(f,M,t))

yinit=(2000*np.pi,M0)
A=integrate.odeint(h,yinit,t)

#loglog(t,2*pi/(A.T[0]))
#loglog(t, (rm(A.T[1],t)/rc(A.T[0],A.T[1],t))**(1.5))
print M
#show()