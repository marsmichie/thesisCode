import numpy as np
import math

# constants from stern_gerlach_constants.nb
e_ 	= 1.602e-19
hbar_	= 1.05457e-34
mEl_	= 9.10938e-31
mAg_	= 1.80409e-25
kB_	= 1.38065e-23
Temp_	= 1500

sigma0_	= 20e-5
d1_	= 0.03
d2_	= 0.03
B0_	= 5
b_	= -1.5e3

v0_     = np.sqrt(4*kB_*Temp_/mAg_) # velocity of atoms coming from the oven in y direction
tau0_	= mAg_* sigma0_**2/hbar_
Tmagnet_= d1_/v0_ # classical time spent of atom in magnet
Tafter_	= d2_/v0_ # classical time spent of atom between magnet and screen

g	= 2
gammaEl_= g*e_/mEl_
gammaAg_= g*e_/mAg_

omegaB_	= gammaEl_*B0_

def vcl_(t):
	return gammaEl_*hbar_*b_*t*0.5/mAg_
def zcl_(t):
	return gammaEl_*hbar_*b_*t**2*0.25/mAg_

vM_	= vcl_(Tmagnet_)
zM_	= zcl_(Tmagnet_)
# characteristic length and time
tc	= Tmagnet_
lc  = np.abs(vcl_(Tmagnet_))*tc#np.abs(vM_)*tc

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
omegaB	= omegaB_ * tc
h_m 	= hbar_ / mAg_ * tc / lc**2


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
#print('sig_m \t' + str(sigma_m))
print('T_mag \t' + str(Tmagnet))
print('T_aft \t' + str(Tafter))
print('zM \t' + str(zM))
print('vM \t' + str(vM))
#print('dt \t' + str(dt))
print('v+z/tau\t' + str((vM+zM/tau0)))
#print('4*zM*(zM+vM*tD=1)/sigma0**2\t' + str(4*zM*(zM+vM*1)/sigma0**2))
#print('exp(4*zM*(zM+vM*tD=1)/sigma0**2)\t' + str(np.exp(4*zM*(zM+vM*1)/sigma0**2))) 
