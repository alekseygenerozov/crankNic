import tangoPlot as tp
import numpy as np
import matplotlib.pyplot as plt

alpha = 0.1    # viscosity param
M7    = 1.0    # primary mass (in 10^7 solar masses)
Md7   = 1E-4   # mass inflow  (normalized for eddington limit for ^^ )

params = tp.readParams('params.in')
l = tp.genGrid(params)
r = l*l
isco = r[0]
f = 1 - np.sqrt(isco/r)

sigma = 36.5*(M7*Md7**3/alpha**4*f**3/r**3)**(.2)
nu    = 8.97E-6*(M7**4*Md7**2*alpha**4*r**3*f**2)**(1.0/5.0)
T     = 488*(Md7**2/(alpha*M7*r**(9.0/2.0))*f**2)**(1.0/5.0)
H     = .003*(Md7*f)**(1.0/5.0)*(alpha*M7)**(-.1)*r**(21.0/20.0)
beta  = 1.0/( 1 + .584*M7*H*T**3/sigma )

def onePlot(ax,y,ylabel):
	ax.plot(r,y,'k-',linewidth=2)
	plt.xlabel('gravitational radii')
	plt.ylabel(ylabel)

onePlot(plt.subplot(2,2,1),sigma,'sigma (10^5 g/cm^2)')
onePlot(plt.subplot(2,2,2),nu,'nu/(GM/c)')
onePlot(plt.subplot(2,2,3),T,'T (10^5 K)')
onePlot(plt.subplot(2,2,4),H,'H/(GM/c^2)')

plt.savefig('sigma.ps')

plt.clf()
plt.semilogy(r,beta,'k-',linewidth=2)
plt.xlabel('gravitational radii')
plt.ylabel('log( beta )')
plt.savefig('beta')
