# script for CHSH figure from thesis connected to 'short_form_super_stern_gerlach_alpha.nb'

import numpy as np
import matplotlib.pyplot as plt
import math
#import scienceplots
plt.style.use(['science','bright'])

from constants_stern_gerlach import *

# simulation constants
steps   = int(1e2)
paths   = 10000
dt      = T / steps
sigma_m = np.sqrt(dt*hbar_/(mAg_*lc**2))

phi     = 1*math.pi
epsi    = 0*0.5*math.pi

parts = 20
max_angle = math.pi/2

# magnet
def rhoM(eps,ph,p1,p2,t,d1,d2):
   A1 = np.exp(-4*p1*zcl(t)/sigma0**2)
   A2 = np.exp(-4*p2*zcl(t)/sigma0**2)
   return 2*(2*(A1-1)*(A2-1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(A1*(np.cos(d1)+1)-np.cos(d1)+1)*(A2*(np.cos(d2)-1)-np.cos(d2)-1)+(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)-np.cos(d1)-1)*(A2*(np.cos(d2)+1)-np.cos(d2)+1))
def szAV1M(eps,ph,p1,p2,t,d1,d2):
	#return (np.cos(eps*0.5))**2*Rpm*sz_pm(t,z,angle,AB)+(np.sin(eps*0.5))**2*Rmp*sz_mp(t,z,angle,AB)+.5*np.sin(eps)*Am*(np.exp(1j*phi)*sz_pm(t,z,angle,AB)+np.exp(-1j*phi)*sz_mp(t,z,angle,AB)) 
    #return -(-np.exp(-j*ph)*(1+np.exp(2*j*ph))*1*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*zcl(t)/sigma0**2)+1)\
    #        *(np.exp(-4*p2*zcl(t)/sigma0**2)-1)+(np.cos(eps*0.5)**2)*(-np.cos(d1)-(np.cos(d1)+1)*np.exp(-4*p1*zcl(t)/sigma0**2)+1)\
    #        *(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*zcl(t)/sigma0**2)-1)-(np.sin(eps*0.5)**2)*(np.cos(d1)\
    #            +(np.cos(d1)-1)*np.exp(-4*p1*zcl(t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*zcl(t)/sigma0**2)+1))\
    #        /(2*(np.exp(-j*ph)*(1.+np.exp(2*j*ph))*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(np.exp(-4*p2*zcl(t)/sigma0**2)-1)+(np.cos(eps*0.5)**2)\
    #        *(-np.cos(d1)+(np.cos(d1)+1)*np.exp(-4*p1*zcl(t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*zcl(t)/sigma0**2)-1)+(np.sin(eps*0.5)**2)*(-np.cos(d1)+(np.cos(d1)-1)\
    #        *np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*zcl(t)/sigma0**2)+1)))
    A1 = np.exp(-4*p1*zcl(t)/sigma0**2)
    A2 = np.exp(-4*p2*zcl(t)/sigma0**2)
    return -((-2*(A1+1)*(A2-1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(-A1*(np.cos(d1)+1)-np.cos(d1)+1)*(A2*(np.cos(d2)-1)-np.cos(d2)-1)\
            -(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)+np.cos(d1)+1)*(A2*(np.cos(d2)+1)-np.cos(d2)+1)))/rhoM(eps,ph,p1,p2,t,d1,d2)

def szAV2M(eps,ph,p1,p2,t,d1,d2):
	#return (np.cos(eps*0.5))**2*Rpm*sz_mp(t,z,angle,AB)+(np.sin(eps*0.5))**2*Rmp*sz_pm(t,z,angle,AB)+.5*np.sin(eps)*Am*(np.exp(-1j*phi)*sz_pm(t,z,angle,AB)+np.exp(1j*phi)*sz_mp(t,z,angle,AB)) 
    #return -(-np.exp(-j*ph)*(1+np.exp(2*j*ph))*1*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(np.exp(-4*p2*zcl(t)/sigma0**2)+1)\
    #        +(np.cos(eps*0.5)**2)*(np.cos(d1)-(np.cos(d1)+1)*np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*zcl(t)/sigma0**2)+1)\
    #        -(np.sin(eps*0.5)**2)*(-np.cos(d1)+(np.cos(d1)-1)*np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*zcl(t)/sigma0**2)-1))\
    #        /(2*(np.exp(-j*ph)*(1.+np.exp(2*j*ph))*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(np.exp(-4*p2*zcl(t)/sigma0**2)-1)\
    #            +(np.cos(eps*0.5)**2)*(-np.cos(d1)+(np.cos(d1)+1)*np.exp(-4*p1*zcl(t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*zcl(t)/sigma0**2)-1)\
    #            +(np.sin(eps*0.5)**2)*(-np.cos(d1)+(np.cos(d1)-1)*np.exp(-4*p1*zcl(t)/sigma0**2)-1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*zcl(t)/sigma0**2)+1)))
    A1 = np.exp(-4*p1*zcl(t)/sigma0**2)
    A2 = np.exp(-4*p2*zcl(t)/sigma0**2)        
    return -((-2*(A1-1)*(A2+1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(-A1*(np.cos(d1)+1)+np.cos(d1)-1)*(A2*(np.cos(d2)-1)+np.cos(d2)+1)\
            -(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)-np.cos(d1)-1)*(A2*(np.cos(d2)+1)+np.cos(d2)-1)))/rhoM(eps,ph,p1,p2,t,d1,d2)
   

#after
def rhoA(eps,ph,p1,p2,t,d1,d2):
   A1 = np.exp(-4*p1*(-zM+vM*t)/sigma0**2)
   A2 = np.exp(-4*p2*(-zM+vM*t)/sigma0**2)
   return 2*(2*(A1-1)*(A2-1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(A1*(np.cos(d1)+1)-np.cos(d1)+1)*(A2*(np.cos(d2)-1)-np.cos(d2)-1)+(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)-np.cos(d1)-1)*(A2*(np.cos(d2)+1)-np.cos(d2)+1))

def szAV1A(eps,ph,p1,p2,t,d1,d2):
    #return -(-np.exp(-j*ph)*(1+np.exp(2*j*ph))*1*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*(-zM+vM*t)/sigma0**2)+1)*(np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)+np.cos(eps*0.5)*np.cos(eps*0.5)*(-np.cos(d1)-(np.cos(d1)+1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)-(np.sin(eps*0.5)**2)*(np.cos(d1)+(np.cos(d1)-1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)+1))\
    #    /(2*(np.exp(-j*ph)*(1+np.exp(2*j*ph))*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)*(np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)+(np.cos(eps*.5)**2)*(-np.cos(d1)+(np.cos(d1)+1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)+(np.sin(eps/2)**2)*(-np.cos(d1)+(np.cos(d1)-1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)+1)))
    A1 = np.exp(-4*p1*(-zM+vM*t)/sigma0**2)
    A2 = np.exp(-4*p2*(-zM+vM*t)/sigma0**2)
    return -((-2*(A1+1)*(A2-1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(-A1*(np.cos(d1)+1)-np.cos(d1)+1)*(A2*(np.cos(d2)-1)-np.cos(d2)-1)\
            -(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)+np.cos(d1)+1)*(A2*(np.cos(d2)+1)-np.cos(d2)+1)))/rhoA(eps,ph,p1,p2,t,d1,d2)

def szAV2A(eps,ph,p1,p2,t,d1,d2):
    #return -(-np.exp(-j*ph)*(1+np.exp(2*j*ph))*1*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)\
    #        *(np.exp(-4*p2*(-zM+vM*t)/sigma0**2)+1)+(np.cos(eps*0.5)**2)*(np.cos(d1)-(np.cos(d1)+1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)\
    #       *(np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)+1)-(np.sin(eps*0.5)**2)*(-np.cos(d1)\
    #            +(np.cos(d1)-1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)*(np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1))\
    #        /(2*(np.exp(-j*ph)*(1+np.exp(2*j*ph))*np.sin(d1)*np.sin(d2)*np.sin(eps)*(np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)*(np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)+(np.cos(eps*0.5)**2)\
    #        *(-np.cos(d1)+(np.cos(d1)+1)*np.exp(-4*p1*(-zM+vM*t)/sigma0**2)+1)*(-np.cos(d2)+(np.cos(d2)-1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)-1)+(np.sin(eps*0.5)**2)*(-np.cos(d1)+(np.cos(d1)-1)\
    #        *np.exp(-4*p1*(-zM+vM*t)/sigma0**2)-1)*(-np.cos(d2)+(np.cos(d2)+1)*np.exp(-4*p2*(-zM+vM*t)/sigma0**2)+1)))
    A1 = np.exp(-4*p1*(-zM+vM*t)/sigma0**2)
    A2 = np.exp(-4*p2*(-zM+vM*t)/sigma0**2)        
    return -((-2*(A1-1)*(A2+1)*np.sin(d1)*np.sin(d2)*np.sin(eps)*np.cos(ph)+(np.cos(eps*0.5)**2)*(-A1*(np.cos(d1)+1)+np.cos(d1)-1)*(A2*(np.cos(d2)-1)+np.cos(d2)+1)\
            -(np.sin(eps*0.5)**2)*(A1*(np.cos(d1)-1)-np.cos(d1)-1)*(A2*(np.cos(d2)+1)+np.cos(d2)-1)))/rhoA(eps,ph,p1,p2,t,d1,d2)
    
    

def two_particle_stern_gerlach_correlation(paths,angle1,angle2,random,ep,ph):
    pos     = np.zeros((2,paths)) # [path]
    t       = 0
    if random:
        # choose random starting orientation of spin th0 ph0
        th0_random   = np.random.uniform(0, math.pi, paths)
        theta0_1     = th0_random-angle1
        theta0_2     = th0_random-angle2

    # choose random starting position from
    pos = np.random.normal(0, np.sqrt(0.5)*sigma0, (2,paths)) # particle 1 and 2

    # count particles going up or down
    Npp     = 0.
    Npm     = 0.
    Nmp     = 0.
    Nmm     = 0.

    for i in range(1,steps):
        # generate gaussian increments
        dw = sigma_m * np.random.normal(0.0,1.0,(2,paths))
        # motion through the magnet
        if t < Tmagnet:
            # z - component
            pos[0] = pos[0] + (2*(-vcl(t)+zcl(t)/tau0)*szAV1M(ep,ph,pos[0],pos[1],t,angle1,angle2) - pos[0]/tau0) * dt + dw[0]
            pos[1] = pos[1] + (2*(-vcl(t)+zcl(t)/tau0)*szAV2M(ep,ph,pos[0],pos[1],t,angle1,angle2) - pos[1]/tau0) * dt + dw[1]
        
        # motion after the magnet
        else:
            pos[0] = pos[0] + (2*(-vM+(-zM+vM*t)/tau0)*szAV1A(ep,ph,pos[0],pos[1],t,angle1,angle2) - pos[0]/tau0) * dt + dw[0]
            pos[1] = pos[1] + (2*(-vM+(-zM+vM*t)/tau0)*szAV2A(ep,ph,pos[0],pos[1],t,angle1,angle2) - pos[1]/tau0) * dt + dw[1]
        
        t = t + dt

    for j in range(0,paths):
        if (pos[0][j]>0): # if 1:+
            if (pos[1][j]<0): # if 2:-
                Npm=Npm+1
            else: # if 2:+
                Npp=Npp+1
        if (pos[0][j]<0): # if 1:-
            if (pos[1][j]<0): # if 2:-
                Nmm=Nmm+1
            else: # if 2:+
                Nmp=Nmp+1
    
    return (Npp-Npm-Nmp+Nmm)/paths

Eab = np.zeros(parts+1)
EAb = np.zeros(parts+1)
EaB = np.zeros(parts+1)
EAB = np.zeros(parts+1)
angle = np.zeros(parts+1)

from matplotlib import cm
#color = iter(cm.rainbow(np.linspace(0,1,8)))
ibm_colorblind = iter(['#000000','#648fff','#785ef0','#dc267f','#fe6100','#ffb000'])

fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
# singlet
epsi    = 0.5*math.pi
phi     = 0
for j in range(0,parts+1):
    angle[j] = j/parts*max_angle
    print(angle[j])
    # E(a,b) ... singlet state gives a\cot b = -np.cos\delta
    Eab[j] = two_particle_stern_gerlach_correlation(paths,0.,angle[j],False,epsi,phi)
    # E(a',b) 
    EAb[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],angle[j],False,epsi,phi)
    # E(a,b') 
    EaB[j] = two_particle_stern_gerlach_correlation(paths,0.,3.*angle[j],False,epsi,phi)
    # E(a',b') 
    EAB[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],3.*angle[j],False,epsi,phi)
#c=next(ibm_colorblind)
#l2, = ax1.plot(angle,-(np.cos(angle)),'-',alpha=0.7,c=c,label=r'$E_{Q}(a,b)\ (\mathrm{singlet})$')
#l2, = ax1.plot(angle,-(np.cos(3*angle)),':',alpha=0.7,c=c,label=r'$E_{Q}(a,b)\ (\mathrm{singlet})$')
#l11, = ax1.plot(angle,Eab,'o',c=c,mfc='none',label=r'$E(a,b)$')

#l12, = ax1.plot(angle,EaB,'o',c=c,mfc='none',label=r'$E(a,B)$')
c=next(ibm_colorblind)
l3, = ax1.plot(angle,np.abs(Eab-EaB+EAb+EAB),'o', mfc='none',c=c,label=r'$|S|\ \mathrm{singlet}$')
#l1, = ax1.plot(angle,np.abs(-np.cos(angle)+np.cos(3*angle)-np.cos(angle)-np.cos(angle)),'-',c=c,label=r'$|S_Q|\ \mathrm{singlet}$')
# triplet
epsi    = 0.5*math.pi
phi     = math.pi
for j in range(0,parts+1):
    angle[j] = j/parts*max_angle
    print(angle[j])
    # E(a,b) ... singlet state gives a\cot b = -np.cos\delta
    Eab[j] = two_particle_stern_gerlach_correlation(paths,0,angle[j],False,epsi,phi)
    # E(a',b) 
    EAb[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],angle[j],False,epsi,phi)
    # E(a,b') 
    EaB[j] = two_particle_stern_gerlach_correlation(paths,0,-angle[j],False,epsi,phi)
    # E(a',b') 
    EAB[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],-angle[j],False,epsi,phi)
c=next(ibm_colorblind)
l4, = ax1.plot(angle,np.abs(Eab+EaB-EAb+EAB),'o', mfc='none',c=c,label=r'$|\tilde S|\ \mathrm{triplet}$')
#l44, = ax1.plot(angle,EAb,'o',c=c,label=r'$E_{Ab}\ \mathrm{triplet}$')
#l45, = ax1.plot(angle,-np.cos(3*angle),'-',c=c,label=r'$E^Q_{ab}\ \mathrm{triplet}$')
#l1, = ax1.plot(angle,np.abs(-np.cos(2*angle)-np.cos(-2*angle)+np.cos(2*angle)-np.cos(6*angle)),'-',c=c,label=r'$|S_Q|\ \mathrm{triplet}$')
# separable
epsi = 0
for j in range(0,parts+1):
    angle[j] = j/parts*max_angle
    print(angle[j])
    # E(a,b) ... singlet state gives a\cot b = -np.cos\delta
    Eab[j] = two_particle_stern_gerlach_correlation(paths,0,angle[j],False,epsi,phi)
    # E(a',b) 
    EAb[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],angle[j],False,epsi,phi)
    # E(a,b') 
    EaB[j] = two_particle_stern_gerlach_correlation(paths,0,3*angle[j],False,epsi,phi)
    # E(a',b') 
    EAB[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],3*angle[j],False,epsi,phi)
c=next(ibm_colorblind)
l6, = ax1.plot(angle,np.abs(Eab-EaB+EAb+EAB),'o', mfc='none',c=c,label=r'$|S|\ \mathrm{separable}$')

# mixed
epsi = math.pi*0.0625 
phi = 0
for j in range(0,parts+1):
    angle[j] = j/parts*max_angle
    print(angle[j])
    # E(a,b) ... singlet state gives a\cot b = -np.cos\delta
    Eab[j] = two_particle_stern_gerlach_correlation(paths,0,angle[j],False,epsi,phi)
    # E(a',b) 
    EAb[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],angle[j],False,epsi,phi)
    # E(a,b') 
    EaB[j] = two_particle_stern_gerlach_correlation(paths,0,3*angle[j],False,epsi,phi)
    # E(a',b') 
    EAB[j] = two_particle_stern_gerlach_correlation(paths,2.*angle[j],3*angle[j],False,epsi,phi)
c=next(ibm_colorblind)
l6, = ax1.plot(angle,np.abs(Eab-EaB+EAb+EAB),'o', mfc='none',c=c,label=r'$|S|\ (\epsilon = \frac{\pi}{16})\  \mathrm{partially\ entangled}$')

ax1.axhline(y=2, xmin=0, xmax=2,c='gray',label=r'$\mathrm{upper\ boundary}$')


ax1.set(xlabel=r'${\delta}/{\pi}$', ylabel=r'$|S|,|\tilde S|$')
ax1.set_xlim(-0.1, max_angle+0.1)
ax1.xaxis.set_ticks([0,max_angle])
ax1.set_xticklabels(['0',max_angle/math.pi])

ax1.legend(loc=0,frameon=True,handlelength=1.5)

plt.show()
fig.savefig('Eab_epsilon_0.pdf')
