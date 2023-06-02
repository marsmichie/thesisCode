# script for plotting two free spins before correlation test connected to 'epr_dec22.nb'

import numpy as np
import matplotlib.pyplot as plt
import math
#import scienceplots
plt.style.use(['science','bright']) 

plot_free_spins = True

dt  =   0.6e-2
# simulation constants
h   =   1.#1.05e-34#1.
m   =   1.#1.6e-27#1.
mu, sigma = 0, (dt/m)**(1/2) # mean and standard deviation
steps   = 500
paths   = 1
time    = [i*dt for i in range(0,steps)]

epsi    = 1*math.pi*0.5
delta   = 0.*math.pi # \pi np.singlet, 0 triplet if epsi=\pi/2

cut     = 1 # number of steps not shown
th00    = 0.25*math.pi 
#
# global variables connected to possible magnetic field
#
e       = 1.#1.6e-19#1.
B       = 0.5e1 # in Tesla
q       = e    # charge of particle
omegaC  = 0.5 * q * B / m
g_factor= 1.
gamma   = g_factor * q / m # gyromagnetic factor with g-factor
length_scale= np.sqrt(h/(omegaC * m))
time_scale  = 1./omegaC
print("length scale: "+str(length_scale))
print("time scale: "+str(time_scale))

#
# spin specfic variables
#
radius_particle = 0.8e-15 / length_scale # in SI units
inertia = 1.#m * radius_particle**2 # inetria depending on the chosen length scale
mu, sigmaI = 0, (dt/inertia)**(1/2) # mean and standard deviation ?NOTE: dimension check!
print("value of diffusion constant sigmaI: "+str(sigmaI))

#
# wiener matrix
#
def h_theta(t,p):
    return np.array([np.cos(p),np.sin(p),0])
def h_phi(t,p):
    return np.array([-np.sin(p)*np.cos(t)/np.sin(t),np.cos(p)*np.cos(t)/np.sin(t),1])
def h_chi(t,p):
    return np.array([np.sin(p)/np.sin(t),-np.cos(p)/np.sin(t),0])

#
# global variables for starting values of spin
#
th_e  = 0.5*np.pi
ph_e  = 1.0*np.pi

angle1 = 0.25*math.pi
angle2 = 0.25*math.pi



szAv        = np.zeros((int(steps),2,paths))
plot_pos    = np.zeros((int(steps),2,paths))
tt       = np.zeros(int(steps))
for i in range(steps-1):
    tt[i+1] = tt[i] + dt


def R(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*(2*np.sin(th1)*np.sin(th2)*np.sin(epsi)*np.cos(delta-ph1+ph2)-np.cos(th1-th2)-np.cos(th1+th2)+2*np.cos(epsi)*(np.cos(th1)-np.cos(th2))+2)

def s1x(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./32*h*(np.sin(th2)*np.sin(epsi)*(np.cos(delta-th1+ph2)+np.cos(delta+th1+ph2)+2*np.sin(delta+ph2))+2*np.sin(th1)*(np.cos(th2)*np.cos(ph1)-np.cos(epsi)*(np.cos(th2)*np.sin(ph1)+np.cos(ph1))+np.sin(ph1)))
def s1y(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*h*(np.sin(th2)*np.sin(epsi)*(np.cos(th1)*np.sin(delta +ph2)-np.cos(delta +ph2))+np.sin(th1)*np.cos(ph1)*(np.cos(th2)*np.cos(epsi)-1)+np.sin(th1)*np.sin(ph1)*(np.cos(th2)-np.cos(epsi)))
def s1z(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*h*(np.sin(th1)*np.sin(th2)*np.sin(epsi)*np.sin(delta-ph1+ph2)+np.cos(th1)-np.cos(th2)-np.cos(th1)*np.cos(th2)*np.cos(epsi)+np.cos(epsi))

def s2x(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*h*(np.sin(th1)*np.sin(epsi)*(np.cos(th2)*np.cos(delta-ph1)-np.sin(delta-ph1))+np.sin(th2)*(np.cos(ph2)*(np.cos(th1)+np.cos(epsi))+np.cos(th1)*np.cos(epsi)*np.sin(ph2)+np.sin(ph2)))
def s2y(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*(h*np.sin(th2)*np.sin(ph2)*(np.cos(th1)+np.cos(epsi))-h*(np.sin(th1)*np.sin(epsi)*(np.cos(th2)*np.sin(delta-ph1)+np.cos(delta -ph1))+np.sin(th2)*np.cos(ph2)*(np.cos(th1)*np.cos(epsi)+1)))
def s2z(epsi,delta,angle1,angle2,th1,ph1,th2,ph2):
    return 1./16.*h*(-np.sin(th1)*np.sin(th2)*np.sin(epsi)*np.sin(delta-ph1+ph2)-np.cos(th1)+np.cos(th2)+np.cos(epsi)*(np.cos(th1)*np.cos(th2)-1))

# define the random starting points
def random_starting_points(t,p,c,n):
    tmp_rand = np.random.uniform(0,1,n)
    t[0] = 0.5*np.pi #(0.25*np.pi + 0.5*np.arccos(2*tmp_rand-1)) # 0.25 + ... in order to avoid boundaries at the start
    tmp_rand = np.random.uniform(0,1,n)
    p[0] = 0. #2*math.pi*(tmp_rand-0.5)
    #phi[0] = math.pi + math.pi*(tmp_rand-0.5)
    tmp_rand = np.random.uniform(0,1,n)
    c[0] = 0*2*math.pi*(tmp_rand-0.5)
    return

# arrays for the angles
th1     = np.zeros((steps,paths))
ph1     = np.zeros((steps,paths))
ch1     = np.zeros((steps,paths))
th2     = np.zeros((steps,paths))
ph2     = np.zeros((steps,paths))
ch2     = np.zeros((steps,paths))

sx1     = np.zeros((steps,paths))
sy1     = np.zeros((steps,paths))
sz1     = np.zeros((steps,paths))
sx2     = np.zeros((steps,paths))
sy2     = np.zeros((steps,paths))
sz2     = np.zeros((steps,paths))

def simulate(epsi,delta):
    if epsi < 0.1:
        th00 = np.random.uniform(1.5,math.pi-1.5,paths)
        th1[0] = th00
        th2[0] = math.pi-th00
        ph1[0] = 0*np.random.uniform(0.,2*math.pi-0.1,paths)
        ph2[0] = ph1[0]+math.pi*0.5
    else:
        if delta < 0.5:
            th00 = 0.2
            th1[0] = th00
            th2[0] = math.pi-th1[0]
            ph1[0] = 0*np.random.uniform(0.,2*math.pi-0.1,paths)
            ph2[0] = ph1[0]+math.pi*0.5
        else:
            th00 = np.random.uniform(0.3,math.pi-0.3,paths)
            th1[0] = th00
            th2[0] = math.pi-th00
            ph1[0] = np.random.uniform(0.,2*math.pi-0.1,paths)
            ph2[0] = ph1[0]+math.pi*0.5

    
    ch1[0] = 0*math.pi
    ch2[0] = 0*math.pi
    for j in range(paths):
        print("particle: "+str(j)+"/"+str(paths), end="\r", flush=True)
        for i in range(1,steps):
            # generate 3d random gaussian numbers
            dw  = mu + sigmaI * np.random.randn(3)

            RR  = R(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])

            sx1[i-1,j] = s1x(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR
            sx2[i-1,j] = s2x(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR
            sy1[i-1,j] = s1y(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR
            sy2[i-1,j] = s2y(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR
            sz1[i-1,j] = s1z(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR
            sz2[i-1,j] = s2z(epsi,delta,angle1,angle2,th1[i-1,j],ph1[i-1,j],th2[i-1,j],ph2[i-1,j])/RR

            # velocity updates
            th1[i,j]    = th1[i-1,j]  + (sx1[i-1,j]*np.cos(ph1[i-1,j])+sy1[i-1,j]*np.sin(ph1[i-1,j]) + 0.5 * h/inertia * np.cos(th1[i-1,j]) / np.sin(th1[i-1,j])) * dt + np.dot(h_theta(th1[i-1,j],ph1[i-1,j]),dw)
            ph1[i,j]    = ph1[i-1,j]    + sz1[i-1,j] * dt + np.dot(h_phi(th1[i-1,j],ph1[i-1,j]),dw)
            ch1[i,j]    = ch1[i-1,j]    + (sz1[i-1,j]*np.cos(th1[i-1,j])+np.sin(th1[i-1,j])*(sx1[i-1,j]*np.sin(ph1[i-1,j])-sy1[i-1,j]*np.cos(ph1[i-1,j])))* dt + np.dot(h_chi(th1[i-1,j],ph1[i-1,j]),dw)

            # generate 3d random gaussian numbers
            dw  = mu + sigmaI * np.random.randn(3)

            th2[i,j]    = th2[i-1,j]  + (sx2[i-1,j]*np.cos(ph2[i-1,j])+sy2[i-1,j]*np.sin(ph2[i-1,j]) + 0.5 * h/inertia * np.cos(th2[i-1,j]) / np.sin(th2[i-1,j])) * dt + np.dot(h_theta(th2[i-1,j],ph2[i-1,j]),dw)
            ph2[i,j]    = ph2[i-1,j]    + sz2[i-1,j] * dt + np.dot(h_phi(th2[i-1,j],ph2[i-1,j]),dw)
            ch2[i,j]    = ch2[i-1,j]    + (sz2[i-1,j]*np.cos(th2[i-1,j])+np.sin(th2[i-1,j])*(sx2[i-1,j]*np.sin(ph2[i-1,j])-sy2[i-1,j]*np.cos(ph2[i-1,j])))* dt + np.dot(h_chi(th2[i-1,j],ph2[i-1,j]),dw)

            #
            # restriction for theta to stay between [0,pi]
            #
            if th1[i,j] > 0.97*math.pi :
                #print("high prev "+str(theta[i-1,j])+" curr "+str(theta[i,j]))
                th1[i,j] = math.pi-0.05#th1[i-1,j]
            elif th1[i,j] < 0.05 :
                #print("low prev "+str(th1[i-1,j])+" curr "+str(th1[i,j]))
                th1[i,j] = 0.05#th1[i-1,j]
            #
            if th2[i,j] > 0.97*math.pi :
                #print("high prev "+str(theta[i-1,j])+" curr "+str(theta[i,j]))
                th2[i,j] = math.pi-0.05#th1[i-1,j]
            elif th2[i,j] < 0.05 :
                #print("low prev "+str(th1[i-1,j])+" curr "+str(th2[i,j]))
                th2[i,j] = 0.05#th1[i-1,j]
            #
            # restriction for phi to stay between [-pi,pi]
            #
            if ph1[i,j] > math.pi :
                ph1[i,j] = ph1[i,j] % (math.pi) - math.pi
            elif ph1[i,j] < -math.pi :
                ph1[i,j] = ph1[i,j] % (math.pi) + 0* math.pi
            if ch1[i,j] > math.pi :
                ch1[i,j] = ch1[i,j] % (math.pi) - math.pi
            elif ch1[i,j] < -math.pi:
                ch1[i,j] = ch1[i,j] % (math.pi) + 0*math.pi
            if ph2[i,j] > math.pi :
                ph2[i,j] = ph2[i,j] % (math.pi) - math.pi
            elif ph2[i,j] < -math.pi :
                ph2[i,j] = ph2[i,j] % (math.pi) + 0* math.pi
            if ch2[i,j] > math.pi :
                ch2[i,j] = ch2[i,j] % (math.pi) - math.pi
            elif ch2[i,j] < -math.pi:
                ch2[i,j] = ch2[i,j] % (math.pi) + 0*math.pi

    if plot_free_spins:
        fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(5*1.5,2*1.5))
        lw  = 1.5 # linewidth plots

        tmp_paths = 1
        if paths < 5:
            tmp_paths = paths
        from matplotlib import cm
        color = iter(['#000000','#648fff','#785ef0','#dc267f','#fe6100','#ffb000'])


        if paths > 5:
            color = iter(cm.rainbow(np.linspace(0, 1, tmp_paths)))
        for i in range(tmp_paths):
            c = next(color)
            ax2.plot(tt[cut:-1],sx1[cut:-1,i],'-', linewidth=lw,color=c,label='$s_{1,x}$')
            ax1.plot(-tt[cut:-1],sx2[cut:-1,i],'-', linewidth=lw,color=c,label='$s_{2,x}$')
            c = next(color)
            ax2.plot(tt[cut:-1],sy1[cut:-1,i],'-', linewidth=lw,color=c,label='$s_{1,y}$')
            ax1.plot(-tt[cut:-1],sy2[cut:-1,i],'-', linewidth=lw,color=c,label='$s_{2,y}$')
            c = next(color)
            ax2.plot(tt[cut:-1],sz1[cut:-1,i], linewidth=lw,color=c,label='$s_{1,z}$')
            ax1.plot(-tt[cut:-1],sz2[cut:-1,i], linewidth=lw,color=c,label='$s_{2,z}$')
            """
            ax4.plot(tt[cut:-1],th1[cut:-1,i], linewidth=lw,color=c,label=r'$\vartheta_{1}$')
            ax4.plot(tt[cut:-1],ph1[cut:-1,i],'--', linewidth=lw,color=c,label=r'$\varphi_{1}$')
            ax3.plot(-tt[cut:-1],th2[cut:-1,i], linewidth=lw,color=c,label=r'$\vartheta_{2}$')
            ax3.plot(-tt[cut:-1],ph2[cut:-1,i],'--', linewidth=lw,color=c,label=r'$\varphi_{2}$')
            """
        #plot_lim = 10
        #ax1.set_ylim(-plot_lim,plot_lim)
        #ax2.set_ylim(-plot_lim,plot_lim)


        ax1.set(xlabel='$t$ (a.u.)', ylabel=r'${s}_{2}/\hbar$')
        ax2.set(xlabel='$t$ (a.u.)', ylabel=r'${s}_{1}/\hbar$')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.set_ticks_position("right")
        ax1.set_xticks([-3,-2,-1,0],labels=['3','2','1','0'])
        ax2.set_xticks([0,1,2,3])
        ax1.legend()
        ax2.legend()
        #ax3.legend()
        #ax4.legend()
        plt.tight_layout()
        plt.legend()
        plt.show()
        fig.savefig('free_spin_0'+str(epsi/math.pi)+'_delta_'+str(delta)+'phd.pdf')
    return 


fig, (ax1) = plt.subplots(1, 1, figsize=(5*1,4*1))
lw  = 1.5 # linewidth plots

tmp_paths = 5
if paths < 5:
    tmp_paths = paths
tmp_paths = 1
from matplotlib import cm
#color = iter(cm.rainbow(np.linspace(0,1,8)))
color = iter(['#000000','#648fff','#785ef0','#dc267f','#fe6100','#ffb000'])

#if paths > 5:
#    color = iter(cm.rainbow(np.linspace(0, 1, tmp_paths)))
simulate(0, 0)
c = next(color)
s1xav=np.average(sx1[cut:-1],axis=1)
s2xav=np.average(sx2[cut:-1],axis=1)
s1yav=np.average(sy1[cut:-1],axis=1)
s2yav=np.average(sy2[cut:-1],axis=1)
s1zav=np.average(sz1[cut:-1],axis=1)
s2zav=np.average(sz2[cut:-1],axis=1)
ax1.plot(tt[cut:-1],(np.average(sx1[cut:-1]*sx2[cut:-1]+sy1[cut:-1]*sy2[cut:-1]+sz1[cut:-1]*sz2[cut:-1],axis=1)-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt(np.average(sx1[cut:-1]**2+sy1[cut:-1]**2+sz1[cut:-1]**2,axis=1)-s1xav**2-s1yav**2-s1zav**2)
*np.sqrt(np.average(sx2[cut:-1]**2+sy2[cut:-1]**2+sz2[cut:-1]**2,axis=1)-s2xav**2-s2yav**2-s2zav**2)),'-', linewidth=lw,color=c,label=r'separable')
for i in range(0):
    ax1.plot(tt[cut:-1],((sx1[cut:-1,i]*sx2[cut:-1,i]+sy1[cut:-1,i]*sy2[cut:-1,i]+sz1[cut:-1,i]*sz2[cut:-1,i])-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt((((sx1[cut:-1,i]-s1xav)**2+(sy1[cut:-1,i]-s1yav)**2+(sz1[cut:-1,i]-s1yav)**2)*((sx2[cut:-1,i]-s2xav)**2+(sy2[cut:-1,i]-s2yav)**2+(sz2[cut:-1,i]-s2zav)**2)))),'--', linewidth=lw,color=c,alpha=0.5,label=r'$s_1\cdot s_2$ separable')
#ax1.plot(tt[cut:-1],np.average(sx1[cut:-1]*sx2[cut:-1],axis=1),'--', linewidth=lw,color=c,label='$s_{1,x}s_{2,x}$')
#ax1.plot(tt[cut:-1],np.average(sy1[cut:-1]*sy2[cut:-1],axis=1),':', linewidth=lw,color=c,label='$s_{1,y}s_{2,y}$')
#ax1.plot(tt[cut:-1],np.average(sz1[cut:-1]*sz2[cut:-1],axis=1),'-.', linewidth=lw,color=c,label='$s_{1,z}s_{2,z}$')
simulate(math.pi*0.5, 0)
c = next(color)
s1xav=np.average(sx1[cut:-1],axis=1)
s2xav=np.average(sx2[cut:-1],axis=1)
s1yav=np.average(sy1[cut:-1],axis=1)
s2yav=np.average(sy2[cut:-1],axis=1)
s1zav=np.average(sz1[cut:-1],axis=1)
s2zav=np.average(sz2[cut:-1],axis=1)
ax1.plot(tt[cut:-1],(np.average(sx1[cut:-1]*sx2[cut:-1]+sy1[cut:-1]*sy2[cut:-1]+sz1[cut:-1]*sz2[cut:-1],axis=1)-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt(np.average(sx1[cut:-1]**2+sy1[cut:-1]**2+sz1[cut:-1]**2,axis=1)-s1xav**2-s1yav**2-s1zav**2)
*np.sqrt(np.average(sx2[cut:-1]**2+sy2[cut:-1]**2+sz2[cut:-1]**2,axis=1)-s2xav**2-s2yav**2-s2zav**2)),'-', linewidth=lw,color=c,label=r'triplet')
for i in range(0):
    ax1.plot(tt[cut:-1],((sx1[cut:-1,i]*sx2[cut:-1,i]+sy1[cut:-1,i]*sy2[cut:-1,i]+sz1[cut:-1,i]*sz2[cut:-1,i])-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt((((sx1[cut:-1,i]-s1xav)**2+(sy1[cut:-1,i]-s1yav)**2+(sz1[cut:-1,i]-s1yav)**2)*((sx2[cut:-1,i]-s2xav)**2+(sy2[cut:-1,i]-s2yav)**2+(sz2[cut:-1,i]-0*s2zav)**2)))),'--', linewidth=lw,color=c,alpha=0.5,label=r'$s_1\cdot s_2$ triplet')

#ax1.plot(tt[cut:-1],np.average(sx1[cut:-1]*sx2[cut:-1],axis=1),'--', linewidth=lw,color=c,label='$s_{1,x}s_{2,x}$')
#ax1.plot(tt[cut:-1],np.average(sy1[cut:-1]*sy2[cut:-1],axis=1),':', linewidth=lw,color=c,label='$s_{1,y}s_{2,y}$')
#ax1.plot(tt[cut:-1],np.average(sz1[cut:-1]*sz2[cut:-1],axis=1),'-.', linewidth=lw,color=c,label='$s_{1,z}s_{2,z}$')
simulate(math.pi*0.5, math.pi)
c = next(color)
s1xav=np.average(sx1[cut:-1],axis=1)
s2xav=np.average(sx2[cut:-1],axis=1)
s1yav=np.average(sy1[cut:-1],axis=1)
s2yav=np.average(sy2[cut:-1],axis=1)
s1zav=np.average(sz1[cut:-1],axis=1)
s2zav=np.average(sz2[cut:-1],axis=1)
ax1.plot(tt[cut:-1],(np.average(sx1[cut:-1]*sx2[cut:-1]+sy1[cut:-1]*sy2[cut:-1]+sz1[cut:-1]*sz2[cut:-1],axis=1)-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt(np.average(sx1[cut:-1]**2+sy1[cut:-1]**2+sz1[cut:-1]**2,axis=1)-s1xav**2-s1yav**2-s1zav**2)
*np.sqrt(np.average(sx2[cut:-1]**2+sy2[cut:-1]**2+sz2[cut:-1]**2,axis=1)-s2xav**2-s2yav**2-s2zav**2)),'-', linewidth=lw,color=c,label=r'singlet')
for i in range(0):
    ax1.plot(tt[cut:-1],((sx1[cut:-1,i]*sx2[cut:-1,i]+sy1[cut:-1,i]*sy2[cut:-1,i]+sz1[cut:-1,i]*sz2[cut:-1,i])-(s1xav*s2xav+s1yav*s2yav+s1zav*s2zav))/(np.sqrt((((sx1[cut:-1,i]-s1xav)**2+(sy1[cut:-1,i]-s1yav)**2+(sz1[cut:-1,i]-s1yav)**2)*((sx2[cut:-1,i]-s2xav)**2+(sy2[cut:-1,i]-s2yav)**2+(sz2[cut:-1,i]-s2zav)**2)))),'--', linewidth=lw,color=c,alpha=0.5,label=r'$s_1\cdot s_2$ singlet')


ax1.set(xlabel=r'$t/\tau_{rot}$', ylabel=r'$C(s_1(t),s_2(t))$')
#ax1.set(xlabel=r'$t/\tau_{rot}$', ylabel=r'$\frac{s_1(t)\cdot s_2(t)}{|s_1(t)||s_2(t)|}$')
ax1.legend()
plt.tight_layout()
plt.legend()
plt.show()
fig.savefig('free_spin_corr_'+str(epsi/math.pi)+'_delta_'+str(delta)+'th0_'+str(th00[0])+'.pdf')