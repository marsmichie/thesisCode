# script for bell inequality test connected to 'short_form_super_stern_gerlach_alpha.nb'
# for the separable problem

import numpy as np
import matplotlib.pyplot as plt
import math
#import scienceplots
plt.style.use(['science','bright'])
from constants_stern_gerlach import *


plot_zoom = True

# simultaion constants
steps   = int(1e2)
paths   = 100
dt      = T / steps
sigma_m = np.sqrt(dt*hbar_/(mAg_*lc**2))

def two_particle_stern_gerlach(paths,theta0_1,theta0_2,random):
    pos     = np.zeros((2,paths)) # [path]
    t       = 0
    if random:
        # choose random starting orientation of spin th0 ph0
        theta0_1     = np.random.uniform(0, math.pi, paths)
        theta0_2     = np.pi - theta0_1 #np.random.uniform(0, math.pi, paths)
    ph0     = np.random.uniform(0, 2*math.pi, paths)
    

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
            pos[0] = pos[0] + (2*(-vcl(t)+zcl(t)/tau0)*OphAV(t,pos[0],theta0_1) - pos[0]/tau0) * dt + dw[0]
            pos[1] = pos[1] + (2*(-vcl(t)+zcl(t)/tau0)*OphAV(t,pos[1],theta0_2) - pos[1]/tau0) * dt + dw[1]

        # motion after the magnet
        else:
            pos[0] = pos[0] + (2*(-vM+(zM+vM*t)/tau0)*OphAV_magnet(t,pos[0],theta0_1) - pos[0]/tau0) * dt + dw[0]
            pos[1] = pos[1] + (2*(-vM+(zM+vM*t)/tau0)*OphAV_magnet(t,pos[1],theta0_2) - pos[1]/tau0) * dt + dw[1]
        
        t = t + dt


    for j in range(0,paths):
        if (pos[0][j]>0): # if 1:+
            if (pos[1][j]>0): # NOTE: switched 'direction', if 2:-
                Npm=Npm+1
            else: # if 2:+
                Npp=Npp+1
        if (pos[0][j]<0): # if 1:-
            if (pos[1][j]>0): # NOTE: switched 'direction', if 2:-
                Nmm=Nmm+1
            else: # if 2:+
                Nmp=Nmp+1
    
    return [Npp/paths,Npm/paths,Nmp/paths,Nmm/paths]

parts = 5
rho_pp = np.zeros(parts+1)
rho_pm = np.zeros(parts+1)
rho_mp = np.zeros(parts+1)
rho_mm = np.zeros(parts+1)
angle = np.zeros(parts+1)
rho_mm = np.zeros(parts+1)
rho_AB = np.zeros((parts+1,4))
rho_BC = np.zeros((parts+1,4))
rho_AC = np.zeros((parts+1,4))

from matplotlib import cm
#color = iter(cm.rainbow(np.linspace(0,1,8)))
ibm_colorblind = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])#iter(['#000000','#648fff','#785ef0','#dc267f','#fe6100','#ffb000'])
for j in range(0,parts+1):
    angle[j] = j/parts*math.pi*0.54
    
    # P(A,Not(B)) ... singlet state gives 0.5*sin²(th/2)
    rho_AB[j] = two_particle_stern_gerlach(paths,0,angle[j],False)
    #rho_AB[j] = 0.5*(rho_AB[j]+two_particle_stern_gerlach(paths,np.pi,math.pi-angle[j],False))
    # P(B,Not(C)) ... singlet state gives 0.5*sin²(th/2)
    rho_BC[j] = two_particle_stern_gerlach(paths,angle[j],2*angle[j],False)
    #rho_BC[j] = 0.5*(rho_BC[j]+two_particle_stern_gerlach(paths,np.pi-angle[j],np.pi-2*angle[j],False))
    # P(A,Not(C)) ... singlet state gives 0.5*sin²(th)
    rho_AC[j] = two_particle_stern_gerlach(paths,0,2*angle[j],False)
    #rho_AC[j] = 0.5*(rho_AC[j]+two_particle_stern_gerlach(paths,np.pi,np.pi-2*angle[j],False))
    print("AB: "+str(rho_AB))
    print("BC: "+str(rho_BC))
    
np.savetxt('bell_no_entanglement.txt',(angle,rho_AB[:,0],rho_BC[:,0],rho_AC[:,0]),delimiter=',')
angle/math.pi
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(ibm_colorblind)
#l1, = ax1.plot(angle,rho_AC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar C}=\frac{N^{A\bar C}}{N}$')
l2, = ax1.plot(angle,np.sin(angle)**2,linewidth=1.5,c=c,label=r'$\rho^{A\bar C}_{Q}$')
c=next(ibm_colorblind)
#l3, = ax1.plot(angle,rho_AB[:,0]+rho_BC[:,0],'<',c=c,label=r'$\rho^{A\bar B}+\rho^{B\bar C}$')
l4, = ax1.plot(angle,2*np.sin(angle*0.5)**2,linewidth=1.5, c=c,label=r'$\rho^{A\bar B}_{Q}+\rho^{B\bar C}_{Q}$')
c=next(ibm_colorblind)
#l5, = ax1.plot(angle,rho_AB[:,0],'o',c=c,label=r'$\rho^{A\bar B}=\frac{N^{A\bar B}}{N}$')
#l6, = ax1.plot(angle,np.sin(angle*0.5)**2,ls=':',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{A\bar B}_{Q}$')
#c=next(ibm_colorblind)
#l7, = ax1.plot(angle,rho_BC[:,0],'o',c=c,mfc='none',label=r'$\rho^{B\bar C}=\frac{N^{B\bar C}}{N}$')
#l8, = ax1.plot(angle,np.sin(angle*0.5)**2,ls='--',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{B\bar C}_Q$')

ax1.set(xlabel=r'$\frac{\delta}{\pi}$', ylabel=r'$\mathrm{normalized\ probability}$')
ax1.set_ylim(-0.1, 1.2)
ax1.set_xlim(-0.2, math.pi*0.5+0.2)
ax1.xaxis.set_ticks([0,math.pi*0.5])
ax1.set_xticklabels(['0','0.5'])
# Create a legend for the first line.
#first_legend = ax1.legend(handles=[l1,l3,l5,l7], loc='upper left',title="QHE separable",frameon=True,handlelength=1.)
# Add the legend manually to the Axes.
#ax1.add_artist(first_legend)
# Create another legend for the second line.
ax1.legend(handles=[l2,l4], loc='upper center',bbox_to_anchor=(0.43,1.0),title="singlet state",frameon=True,handlelength=1.)

plt.show()
fig.savefig('bell_rho_pp_0.pdf')


# 2
ibm_colorblind = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(ibm_colorblind)
l1, = ax1.plot(angle,rho_AC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar C}=\frac{N^{A\bar C}}{N}$')
l2, = ax1.plot(angle,np.sin(angle)**2,linewidth=1.5,c=c,label=r'$\rho^{A\bar C}_{Q}$')
c=next(ibm_colorblind)
l3, = ax1.plot(angle,rho_AB[:,0]+rho_BC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar B}+\rho^{B\bar C}$')
l4, = ax1.plot(angle,2*np.sin(angle*0.5)**2,linewidth=1.5, c=c,label=r'$\rho^{A\bar B}_{Q}+\rho^{B\bar C}_{Q}$')
c=next(ibm_colorblind)
#l5, = ax1.plot(angle,rho_AB[:,0],'o',c=c,label=r'$\rho^{A\bar B}=\frac{N^{A\bar B}}{N}$')
#l6, = ax1.plot(angle,np.sin(angle*0.5)**2,ls=':',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{A\bar B}_{Q}$')
#c=next(ibm_colorblind)
#l7, = ax1.plot(angle,rho_BC[:,0],'o', mfc='none',c=c,mfc='none',label=r'$\rho^{B\bar C}=\frac{N^{B\bar C}}{N}$')
#l8, = ax1.plot(angle,np.sin(angle*0.5)**2,ls='--',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{B\bar C}_Q$')

ax1.set(xlabel=r'$\frac{\delta}{\pi}$', ylabel=r'$\mathrm{normalized\ probability}$')
ax1.set_ylim(-0.1, 1.2)
ax1.set_xlim(-0.2, math.pi*0.5+0.2)
ax1.xaxis.set_ticks([0,math.pi*0.5])
ax1.set_xticklabels(['0','0.5'])
# Create a legend for the first line.
first_legend = ax1.legend(handles=[l1,l3], loc='upper left',title="QHE separable",frameon=True,handlelength=1.)
# Add the legend manually to the Axes.
ax1.add_artist(first_legend)
# Create another legend for the second line.
ax1.legend(handles=[l2,l4], loc='upper center',bbox_to_anchor=(0.43,1.0),title="singlet state",frameon=True,handlelength=1.)

plt.show()
fig.savefig('bell_rho_pp_1.pdf')

# 3
ibm_colorblind = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(ibm_colorblind)
l1, = ax1.plot(angle,rho_AC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar C}=\frac{N^{A\bar C}}{N}$')
l2, = ax1.plot(angle,np.sin(angle)**2,linewidth=1.5,c=c,label=r'$\rho^{A\bar C}_{Q}$')
c=next(ibm_colorblind)
l3, = ax1.plot(angle,rho_AB[:,0]+rho_BC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar B}+\rho^{B\bar C}$')
l4, = ax1.plot(angle,2*np.sin(angle*0.5)**2,linewidth=1.5, c=c,label=r'$\rho^{A\bar B}_{Q}+\rho^{B\bar C}_{Q}$')
c=next(ibm_colorblind)
l5, = ax1.plot(angle,rho_AB[:,0],'o',mfc='none',c=c,label=r'$\rho^{A\bar B}=\frac{N^{A\bar B}}{N}$')
l6, = ax1.plot(angle,np.sin(angle*0.5)**2,ls=':',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{A\bar B}_{Q}$')
#c=next(ibm_colorblind)
#l7, = ax1.plot(angle,rho_BC[:,0],'o',c=c,mfc='none',label=r'$\rho^{B\bar C}=\frac{N^{B\bar C}}{N}$')
#l8, = ax1.plot(angle,np.sin(angle*0.5)**2,ls='--',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{B\bar C}_Q$')

ax1.set(xlabel=r'$\frac{\delta}{\pi}$', ylabel=r'$\mathrm{normalized\ probability}$')
ax1.set_ylim(-0.1, 1.2)
ax1.set_xlim(-0.2, math.pi*0.5+0.2)
ax1.xaxis.set_ticks([0,math.pi*0.5])
ax1.set_xticklabels(['0','0.5'])
# Create a legend for the first line.
first_legend = ax1.legend(handles=[l1,l3,l5], loc='upper left',title="QHE separable",frameon=True,handlelength=1.)
# Add the legend manually to the Axes.
ax1.add_artist(first_legend)
# Create another legend for the second line.
ax1.legend(handles=[l2,l4,l6], loc='upper center',bbox_to_anchor=(0.43,1.0),title="singlet state",frameon=True,handlelength=1.)

plt.show()
fig.savefig('bell_rho_pp_2.pdf')

# 4
ibm_colorblind = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(ibm_colorblind)
l1, = ax1.plot(angle,rho_AC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar C}=\frac{N^{A\bar C}}{N}$')
l2, = ax1.plot(angle,np.sin(angle)**2,linewidth=1.5,c=c,label=r'$\rho^{A\bar C}_{Q}$')
c=next(ibm_colorblind)
l3, = ax1.plot(angle,rho_AB[:,0]+rho_BC[:,0],'o', mfc='none',c=c,label=r'$\rho^{A\bar B}+\rho^{B\bar C}$')
l4, = ax1.plot(angle,2*np.sin(angle*0.5)**2,linewidth=1.5, c=c,label=r'$\rho^{A\bar B}_{Q}+\rho^{B\bar C}_{Q}$')
c=next(ibm_colorblind)
l5, = ax1.plot(angle,rho_AB[:,0],'o',c=c,mfc='none',label=r'$\rho^{A\bar B}=\frac{N^{A\bar B}}{N}$')
l6, = ax1.plot(angle,np.sin(angle*0.5)**2,ls=':',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{A\bar B}_{Q}$')
c=next(ibm_colorblind)
l7, = ax1.plot(angle,rho_BC[:,0],'o', c=c,mfc='none',label=r'$\rho^{B\bar C}=\frac{N^{B\bar C}}{N}$')
l8, = ax1.plot(angle,np.sin(angle*0.5)**2,ls='--',alpha=0.7,linewidth=1.5,c=c,label=r'$\rho^{B\bar C}_Q$')

ax1.set(xlabel=r'$\frac{\delta}{\pi}$', ylabel=r'$\mathrm{normalized\ probability}$')
ax1.set_ylim(-0.1, 1.2)
ax1.set_xlim(-0.2, math.pi*0.5+0.2)
ax1.xaxis.set_ticks([0,math.pi*0.5])
ax1.set_xticklabels(['0','0.5'])
# Create a legend for the first line.
first_legend = ax1.legend(handles=[l1,l3,l5,l7], loc='upper left',title="QHE separable",frameon=True,handlelength=1.)
# Add the legend manually to the Axes.
ax1.add_artist(first_legend)
# Create another legend for the second line.
ax1.legend(handles=[l2,l4,l6,l8], loc='upper center',bbox_to_anchor=(0.43,1.0),title="singlet state",frameon=True,handlelength=1.)

plt.show()
fig.savefig('bell_rho_pp_3.pdf')