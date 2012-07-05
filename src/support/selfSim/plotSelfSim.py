import matplotlib
import matplotlib.pyplot as plt
import numpy
from numpy import log, exp, sqrt, pi
from tangoPlot import *
from scipy.special import erf

matplotlib.rc('lines',markersize=4,linewidth=1);

def plotData(fName,t,params,stroke):
	
	# grab data
	data = readDataFile(fName)
	l = data[:,0]
	FJ = data[:,1]

	FJa = analytic(l,t,params)

	# plot it against analytic
	ax.loglog(l,FJ,stroke,linewidth=1,
	           label='t = '+str(t)+' $\\times\; l_{\\rm in}D_J^{-1}$')
	ax.loglog(l,FJa,'k-',linewidth=.5,label=":)")

	return data

def plotResids(fName,t,params,stroke):
	data = readDataFile(fName)
	l = data[:,0]
	FJ = data[:,1]
	FJa = analytic(l,t,params)
	resids = FJ - FJa
	plt.semilogx(l,resids,stroke)

#
#	ANALYTIC
#		
#		Self-similar solution from Rafikov 2012 
#		paper
#
def analytic(l,t,params):
	
	Mdot = 1.0
	lin = l[0]
	lmlin = l - lin

	chi = params['inner_bndry_value']
	D0 = params['D0']

	FJ = sqrt(4.0*D0*t/pi)*exp(-lmlin*lmlin/(4.0*D0*t))
	FJ += lin + lmlin*erf(lmlin/(2.0*sqrt(D0*t)))
	FJ *= (1.0-chi)*Mdot
	FJ += chi*Mdot*l

	return FJ


###====================================================

#
#   Some IO prep ....
#
fBase = 'savedOutputs/chi_p'
paramBase = 'savedOutputs/params_chi_p'

T = [0,1,10,100]
color = ['b','g','r','c']

params_0 = readParams(paramBase+'0.dat')
params_5 = readParams(paramBase+'5.dat')


#
# Plot Evolution
#

ax = plt.subplot(1,1,1)

for t,c in zip(T,color):
	plotData(fBase + '0_' + strFN(t) + '.dat',t,params_0,c+'-')
	plotData(fBase + '5_' + strFN(t) + '.dat',t,params_5,c+'--')

ax.axis((0.9,35.0,0.9,35.0))
plt.xlabel('$l/l_{\\rm in}$')
plt.ylabel('$F_J/(\\dot{M}_{\\infty}l_{\\rm in})$')
plt.title('Self Similar Solution')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::4],labels[::4],'lower right')
plt.savefig("numerical")
plt.savefig("numerical.ps")

#
# Plot Residuals
#

plt.clf()
for t,c in zip(T,color):
	plotResids(fBase + '0_' + strFN(t) + '.dat',t,params_0,c+'-')
	plotResids(fBase + '5_' + strFN(t) + '.dat',t,params_5,c+'--')

plt.xlim((0.9,35.0))
plt.xlabel('$l / l_{\\rm in}$')
plt.ylabel('residuals')
plt.title('Residuals for Self-similar solution')

plt.savefig("resids")
plt.savefig("resids.ps")
