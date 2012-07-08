#
#	MAKE ANALYTIC
#
#	This program constructs the analytic solution to the progression
#	of a delta function of surface density in an 'alpha disc' as given
#	by Shakura and Sunyaev (height integrated 1-D disk). 
#		It produces a plot of this surface density at 4 representative 
#	times, that should exactly match the plot on page 83 of Frank, King
#	& Raine. The data is also output to text files to compare with
#	simulations
#
#	*** Numpy, scipy and matplotlib are all required
#
#		created by:		Munier Salem
#		on:						October 19, 2011
#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import pow
from scipy import arange, pi, sqrt, exp
from scipy.special import iv
from matplotlib import rc
from math import isnan
from tangoPlot import * 

params = readParams("params.in")
l0 = params['l0']
M  = params['M']
pi = 3.14159
Lambda = params['lambda']
N = int(params['N'])

lMin = params['lMin']
lMax = params['lMax']
Lambda = 1.03
nu = params['D0']*4.0/3.0	# conversion assumes M = 1

print "l0 = " + str(l0)
print "lMin = " + str(lMin)
print "lMax = " + str(lMax)
print "N = " + str(N)
print "lambda = " + str(Lambda)
print "nu = " + str(nu)

l = genGrid("params.in")
dl = l[1] - l[0]

#
# ---- Our analytic expression
#
def Fj(x,t):
		floor = 1E-4
		tmp = 3.0*M*M*M*nu/(l0*l0*l0*t)
		tmp = tmp*x**(.25)*exp(-(1+x*x)/t)*iv(.25,2*x/t) + floor

		for i in range(len(tmp)):
			if( isnan(tmp[i]) ):
				tmp[i] = floor
		return tmp


for i in range(64):	# FIXME
	t = .008 + .008*i
	x = (l/l0)**2
	F_j = Fj(x,t)

#	plt.plot(l,f)

	tmp = ""
	fNum = i
	if( fNum < 10 ):
		tmp = "00" + str(fNum)
	elif( fNum < 100 ):
		tmp = "0" + str(fNum)
	else:
		tmp = str(fNum)

	f = open('analytic/T' + tmp + '.dat','w')
	for i in range(N):
		f.write(str(l[i]) + '\t' + str(F_j[i]) + '\n')
	f.close()
