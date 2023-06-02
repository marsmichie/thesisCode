# script for generating the stern gerlach figure 6.7

import numpy as np
import matplotlib.pyplot as plt
import math
#import scienceplots
plt.style.use(['science','bright'])

draw_arrows     = True
classical_path  = False
plotting_paths  = True
plot_zoom       = True
random_plot     = False
draw_quiver     = False
draw_velo       = False

# simulation constants
steps   = int(1e3)
paths   = 100
number_of_arrows = int(steps/10)


# constants from stern_gerlach_constants.nb
e_ 	= 1.602e-19
hbar_	= 1.05457e-34
mEl_	= 9.10938e-31
mAg_	= 1.80409e-25
kB_	= 1.38065e-23
Temp_	= 1500

sigma0_	= 20e-5
d1_	= 0.03
d2_	= 0.06
B0_	= 5
b_	= -1.5e3

v0_     = np.sqrt(4*kB_*Temp_/mAg_) # velocity of atoms coming from the oven in y direction
tau0_	= mAg_* sigma0_**2/hbar_
Tmagnet_= d1_/v0_ # classical time spent of atom in magnet
Tafter_	= d2_/v0_ # classical time spent of atom between magnet and screen

g	= 2
gammaEl_= g*e_/mEl_
gammaAg_= g*e_/mAg_

def vcl_(t):
	return gammaEl_*hbar_*b_*t*0.5/mAg_
def zcl_(t):
	return gammaEl_*hbar_*b_*t**2*0.25/mAg_

vM_	= vcl_(Tmagnet_)
zM_	= zcl_(Tmagnet_)
# characteristic length and time
tc	= Tmagnet_
lc  = np.abs(vM_)*tc

# dimensionless constants
v0      = tc / lc * v0_
Tmagnet	= Tmagnet_/tc
Tafter	= Tafter_/tc
T       = Tmagnet + Tafter
tau0	= tau0_/tc
vM	    = tc / lc * vM_
sigma0  = sigma0_ / lc
d1      = d1_ / lc
d2      = d2_ / lc
zM	    = zM_ / lc
dt      = T / steps
sigma_m = np.sqrt(dt*hbar_/(mAg_*lc**2))


# in magnet
def vcl(tD):
	return vcl_(tD*tc)*tc/lc
def zcl(tD):
	return zcl_(tD*tc)/lc
def OphAV(tD,zD,th0):
	return -0.5-(1+np.cos(th0))/(-1+np.exp(4*zD*zcl(tD)/sigma0**2)*(-1+np.cos(th0))-np.cos(th0))

# after magnet
def OphAV_magnet(tD,zD,th0):
	return -0.5-(1+np.cos(th0))/(-1+np.exp(4*zD*(-zM+vM*tD)/sigma0**2)*(-1+np.cos(th0))-np.cos(th0))


print('v0 \t' + str(v0))
print('s0 \t' + str(sigma0))
print('d1 \t' + str(d1))
print('d2 \t' + str(d2))
print('sig_m \t' + str(sigma_m))
print('T_mag \t' + str(Tmagnet))
print('T_aft \t' + str(Tafter))
print('zM \t' + str(zM))
print('vM \t' + str(vM))
print('dt \t' + str(dt))
print('v+z/tau\t' + str((vM+zM/tau0)))
#print('4*zM*(zM+vM*tD=1)/sigma0**2\t' + str(4*zM*(zM+vM*1)/sigma0**2))
#print('exp(4*zM*(zM+vM*tD=1)/sigma0**2)\t' + str(np.exp(4*zM*(zM+vM*1)/sigma0**2)))


pos     = np.zeros((2,int(steps),paths)) # [yz][step][path]
cl_pos  = np.zeros((2,int(steps),paths)) # [yz][step][path]
t       = np.zeros(int(steps),dtype=np.float64)
for i in range(steps-1):
    t[i+1] = t[i] + dt


def stern_gerlach(theta0,random):
    if random:
        # choose random starting orientation of spin th0 ph0
        th0     = np.random.uniform(0, math.pi, paths)
    else:
        th0     = theta0
    ph0     = np.random.uniform(0, 2*math.pi, paths)

    # values of th(z,t) and ph(z,t) for 'spin projection'
    th      = np.zeros((int(steps),paths))
    th[0]   = th0
    szAv    = np.zeros((int(steps),paths))

    # choose random starting position from
    pos[0][0] = np.random.normal(0, 0.5*sigma0, paths) # y position
    pos[1][0] = np.random.normal(0, np.sqrt(0.5)*sigma0, paths)

    cl_pos[0][0] = pos[0][0]
    cl_pos[1][0] = pos[1][0]

    szAv[0] = OphAV(0,pos[1][0],th[0])

    # generate gaussian increments
    dw = sigma_m * np.random.normal(0.0,1.0,(2,int(steps),paths))

    # count particles going up or down
    Nplus   = 0.
    Nminus  = 0.

    for i in range(1,steps):
        # y - component
        pos[0][i] = pos[0][i-1] + (v0 + (-v0 * t[i-1] -  pos[0][i-1]) / tau0) * dt + dw[0][i]
        cl_pos[0][i] = cl_pos[0][i-1] + (v0) * dt
        #
        #  motion through the magnet
        if t[i] < Tmagnet:
            # z - component
            pos[1][i] = pos[1][i-1] + (2*(-vcl(t[i-1])+zcl(t[i-1])/tau0)*OphAV(t[i-1],pos[1][i-1],th[0]) - pos[1][i-1]/tau0) * dt + dw[1][i]
            cl_pos[1][i] = cl_pos[1][i-1] + (2*(-vcl(t[i-1]))*szAv[0]) * dt
            # spin projection 
            szAv[i] = OphAV(t[i],pos[1][i],th0)

        # motion after the magnet
        else:
            # ftz = (vB*t[i] + delta_z) * pos[2][i] / (d0**2)
            pos[1][i] = pos[1][i-1] + (2*(-vM+(zM+vM*t[i-1])/tau0)*OphAV_magnet(t[i-1],pos[1][i-1],th[0]) - pos[1][i-1]/tau0) * dt + dw[1][i]
            cl_pos[1][i] = cl_pos[1][i-1] + (2*(-vM)*szAv[0]) * dt 
            # spin projection angles
            szAv[i] = OphAV_magnet(t[i],pos[1][i],th0)

    for j in range(0,paths):
        if (pos[1][-1][j]>0):
            Nplus=Nplus+1
        if (pos[1][-1][j]<0):
            Nminus=Nminus+1
    print(str(theta0/math.pi)+'$\pi$')
    print('N+/N ',Nplus/paths)
    print('N-/N ',Nminus/paths)

    

    if plotting_paths:
        tmp_paths = 8
        if paths < 8:
            tmp_paths = paths
        y1_max = 1
        y2_max = 0.5
        if draw_quiver: # spin expectation plot
            fig, ax = plt.subplots(figsize =(4, 3.2))
            arrows = 30
            xx = np.arange(0,t[-1],(t[-1]/arrows))
            yy = np.arange(-y1_max/(lc*1000),y1_max/(lc*1000),(2*y1_max/(lc*1000*arrows)))
            XX, YY = np.meshgrid(xx, yy)
            sz    = -0.5*(1.0-(2*(1+np.cos(theta0))/(1+np.exp(4*zcl(XX)*YY/sigma0**2)*(1-np.cos(theta0))+np.cos(theta0))))
            y_dir = np.sqrt(0.25-sz**2)
            z_dir = sz
            ax.quiver(XX/Tmagnet, YY*1000*lc, y_dir*1e-1, z_dir*1e-1, scale=2,color='#000000', alpha=0.9, width=0.003)
            for i in range(paths):
                c = '#56B4E9'
                ax.plot(t/Tmagnet,pos[1,:,i]*lc*1000.0, linewidth=1.5,color=c,alpha=0.6)

            ax.set(xlabel='$t/T_m$ ', ylabel='$z$ (mm)')
            ax.set_xlim(0, 2)
            ax.set_ylim(-y1_max,y1_max)
            ax.set_xticks([0,1,2])
            ax.set_yticks([-y1_max,0,y1_max])
            plt.title('')#traj: %i, steps: %i, dim: %i, neur_p_l: %i, num_l: %i'% (M, N, D, neuron_per_layer, number_of_deep_layers))
            plt.legend()
            plt.tight_layout()
            plt.show()
            fig.savefig('sg_expectation'+str(np.pi/theta0)+'.pdf')

        if draw_velo: # spin expectation plot
            y1_max = 1.5
            fig, ax = plt.subplots(figsize =(4, 3.2))
            arrows = 30*2/2.8
            xx = np.arange(0,t[-1]*2,(t[-1]/arrows))
            yy = np.arange(-y1_max/(lc*1000),y1_max/(lc*1000),(2*y1_max/(lc*1000*arrows)))
            XX, YY = np.meshgrid(xx, yy)
            vy = v0 
            sz = -0.5*(1.0-(2*(1+np.cos(theta0))/(1+np.exp(4*zcl(XX)*YY/sigma0**2)*(1-np.cos(theta0))+np.cos(theta0))))
            vz = -2*vcl(XX)*sz*19
            uy = 0
            uz = -YY-2*zcl(XX)*sz
            ax.quiver(XX/Tmagnet, YY*1000*lc, vy, vz, scale=5500,color='#000000', alpha=0.7, width=0.003)
            ax.quiver(XX/Tmagnet, YY*1000*lc, uy, uz, scale=25,color='#CC79A7', alpha=0.9, width=0.003,scale_units='x')
            for i in range(paths):
                c = '#56B4E9'
                ax.plot(t/Tmagnet,pos[1,:,i]*lc*1000.0, linewidth=1.5,color=c,alpha=0.6)

            ax.set(xlabel='$t/T_m$ ', ylabel='$z$ (mm)')
            ax.set_xlim(0, 2.7)
            ax.set_ylim(-y1_max,y1_max)
            ax.set_xticks([0,2.7],labels=[0,1])
            ax.set_yticks([-1,0,1])
            plt.title('')#traj: %i, steps: %i, dim: %i, neur_p_l: %i, num_l: %i'% (M, N, D, neuron_per_layer, number_of_deep_layers))
            plt.legend()
            plt.tight_layout()
            plt.show()
            fig.savefig('sg_velo'+str(theta0)+'.pdf')

        # new plot
        # plotting the paths
        fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(8,3.5))

        y1_max = 1
        y2_max = 0.5

        from matplotlib import cm
        color = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
        for i in range(tmp_paths):
            c = next(color)
            if True: # z(y) or z(t)
                ax2.plot(pos[0,:,i]*lc*100.0,pos[1,:,i]*lc*1000.0, linewidth=1.5,color=c)
                if classical_path:
                    ax2.plot(cl_pos[0,:,i]*lc*100.0,cl_pos[1,:,i]*lc*1000.0, ls=':',linewidth=1,color=c)
                ax2.axvline(x=d1_*100, color='grey', ls = '--', linewidth=0.5)
            else:
                ax2.plot(t/Tmagnet,pos[1,:,i]*lc*1000.0, linewidth=1.5,color=c)
                ax2.axvline(x=1, color='grey', ls = '--', linewidth=0.5)
            #ax1.text(2.0,.044, r'$B\neq 0$')
            #ax1.text(3.2,.044, r'$B=0$')
            #ax2.plot(pos[0,:,i]*lc*100.0,szAv[:,i], linewidth=2,color=c)
            
            ax1.plot(t/Tmagnet,szAv[:,i], linewidth=1.5,color=c)
            ax1.axvline(x=1, color='grey', ls = '--', linewidth=0.5)
            
            if draw_arrows:
                y_pos = pos[0,0::number_of_arrows,i]*lc*100.0
                z_pos = pos[1,0::number_of_arrows,i]*lc*1000.0
                y_dir = np.sqrt(0.25-szAv[0::number_of_arrows,i]**2)
                z_dir = szAv[0::number_of_arrows,i]
                ax2.quiver(y_pos, z_pos, y_dir*1e-1, z_dir*1e-1, scale=0.7,color=c, alpha=0.8, width=0.008)

        ax2.set(xlabel='$y$ (cm)', ylabel='$z$ (mm)')
        ax2.set_xticks([0,3,6,9])
        ax2.set_yticks([-y1_max,0,y1_max])
        '''
        y1_max = 0.1
        y2_max = 0.5
        # plotting the paths zoom
        color = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
        if plot_zoom:
            plot_small=int(d1/(d1+d2)*steps)
            arr_num_small = int(plot_small/10)
            for i in range(tmp_paths):
                c = next(color)
                ax1.plot(pos[0,:plot_small,i]*lc*100.0,pos[1,:plot_small,i]*lc*10000.0, linewidth=1.5,color=c)
                if classical_path:
                    ax1.plot(cl_pos[0,:plot_small,i]*lc*100.0,cl_pos[1,:plot_small,i]*lc*10000.0, ls=':',linewidth=1,color=c)
                if draw_arrows:
                    y_pos = pos[0,0:plot_small:arr_num_small,i]*lc*100.0
                    z_pos = pos[1,0:plot_small:arr_num_small,i]*lc*10000.0
                    y_dir = np.sqrt(0.25-szAv[0:plot_small:arr_num_small,i]**2)
                    z_dir = szAv[0:plot_small:arr_num_small,i]
                    ax1.quiver(y_pos, z_pos, y_dir*1e-1, z_dir*1e-1, scale=0.7,color=c, alpha=0.8, width=0.008)

            ax1.set_yticks([-4,0,4])
            ax1.set(xlabel='$y$ (cm)', ylabel=r'$z$ ($10^{-4}$m)')
            ax1.set_xlim(-0.1, d1_*100+0.1)
        '''
        
        ax1.set(xlabel='$t/T_m$', ylabel=r'$\bar{s}_z(t)/\hbar$')
        ax1.set_yticks([-y2_max,0,y2_max],labels=[r'$-\frac{1}{2}$',r'$0$',r'$\frac{1}{2}$'])
        ax1.set_xticks([0,1,2])
        ax1.set_xlim(0, 2)
        
        plt.tight_layout()
        plt.title('')#traj: %i, steps: %i, dim: %i, neur_p_l: %i, num_l: %i'% (M, N, D, neuron_per_layer, number_of_deep_layers))
        plt.legend()
        plt.show()
        if random:
            fig.savefig('sg_random.pdf')
        else:
            fig.savefig('sg_'+str(theta0)+'.pdf')


        color = iter(['#000000','#56B4E9','#CC79A7','#009E73','#E69F00','#F0E442','#0072B2','#D55E00'])
        fig, (ax2) = plt.subplots(1, 1, figsize=(7*0.8,4*0.8))
        # plotting the paths zoom
        for i in range(tmp_paths):
            c = next(color)
            ax2.plot(t/Tmagnet,szAv[:,i], linewidth=1.5,color=c)
        ax2.axvline(x=1, color='grey', ls = '--', linewidth=0.5)
        ax2.set(xlabel='$t/T_m$', ylabel=r'$\bar{s}_z(t)/\hbar$')
        ax2.set_yticks([-y2_max,0,y2_max],labels=[r'$-\frac{1}{2}$',r'$0$',r'$\frac{1}{2}$'])
        ax2.set_xticks([0,1,2])
        ax2.set_xlim(0, 2)

        
        plt.title('')#traj: %i, steps: %i, dim: %i, neur_p_l: %i, num_l: %i'% (M, N, D, neuron_per_layer, number_of_deep_layers))
        plt.legend()
        plt.show()
        if random:
            fig.savefig('sg_expectation_random.pdf')
    
    return [Nplus/paths,Nminus/paths]

if random_plot:
    stern_gerlach(0.5*np.pi,True)

if draw_quiver:
    stern_gerlach(0.2*np.pi,False)
if draw_quiver:
    stern_gerlach(0.5*np.pi,False)


parts = 4
rho_plus = np.zeros(parts+1)
rho_minus = np.zeros(parts+1)
angle = np.zeros(parts+1)

for j in range(2,parts+1):
    angle[j] = j/parts*math.pi
    [rho_plus[j],rho_minus[j]] = stern_gerlach(angle[j],False)
    
np.savetxt('rho_plus_minus.txt',(angle,rho_plus,rho_minus),delimiter=',')

from matplotlib import cm
#color = iter(cm.rainbow(np.linspace(0,1,8)))
color = iter(['#000000','#dc267f','#648fff','#785ef0','#fe6100','#ffb000'])
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(color)
ax1.plot(angle/math.pi,rho_plus,'o',c=c,label=r'$\rho_+=\frac{N_+}{N}$')
ax1.plot(angle/math.pi,np.cos(angle*0.5)**2,'--',linewidth=1.5,c=c,label=r'$\rho^{QM}_+$')
ax1.set(xlabel=r'$\delta/\pi$', ylabel=r'$\mathrm{normalized\ probability}$')
plt.legend()
plt.show()
color = iter(['#000000','#dc267f','#648fff','#785ef0','#fe6100','#ffb000'])
fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))
c=next(color)
ax1.plot(angle/math.pi,rho_plus,'o',c=c,label=r'$\rho_+=\frac{N_+}{N}$')
ax1.plot(angle/math.pi,np.cos(angle*0.5)**2,'--',linewidth=1.5,c=c,label=r'$\rho^{QM}_+$')
ax1.set(xlabel=r'$\delta/\pi$', ylabel=r'$\mathrm{normalized\ probability}$')
c=next(color)
ax1.plot(angle/math.pi,rho_minus,'o',c=c,label=r'$\rho_-=\frac{N_-}{N}$')
ax1.plot(angle/math.pi,np.sin(angle*0.5)**2,'--',linewidth=1.5,c=c,label=r'$\rho^{QM}_-$')
plt.legend()
plt.show()
fig.savefig('rhoplusminus.pdf')
