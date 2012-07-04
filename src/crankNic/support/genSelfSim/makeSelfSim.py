import numpy as np
from numpy import log, exp, sqrt, pi
from tangoPlot import *
from scipy.special import erf

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

	FJ = np.zeros(l.size)

	if( t != 0.0 ):
		FJ += sqrt(4.0*D0*t/pi)*exp(-lmlin*lmlin/(4.0*D0*t))
		FJ += lin + lmlin*erf(lmlin/(2.0*sqrt(D0*t)))
	else:
		FJ += l
	FJ *= (1.0-chi)*Mdot
	FJ += chi*Mdot*l

	return FJ

params = readParams("params.in")
l = genGrid("params.in")
FJ = analytic(l,params['tStart'],params)

for ll, ff in zip(l,FJ):
	print str(ll) + "\t" + str(ff)
