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

R0 = 1.0 		# initial radius of delta fcn
M  = 1.0		# intiial mass
pi = 3.14159

N = 1000
rMin = 0.1
rMax = 2.0
Lambda = 1.0

# we read in from param file ...
f = open('params.in','r')
params = f.read()
f.close()
params = params.split('\n');
Np = len(params)

for i in range(Np):
	s = params[i]
	s = s.split('=')
	var = s[0].strip()
	if( len(s) > 1 ):
		num = s[1].split('//')
		num = float(num[0].strip())
	
		if( var == 'r0' ): 	R0 = num
		if( var == 'rMin' ): rMin = num
		if( var == 'rMax' ): rMax = num
		if( var == 'N' ): N = int(num)
		if( var == 'lambda' ): Lambda = num

print "R0 = " + str(R0)
print "rMin = " + str(rMin)
print "rMax = " + str(rMax)
print "N = " + str(N)
print "lambda = " + str(Lambda)

#
# ---- Our analytic expression
#
def S(x,t):
		floor = 1E-4
		tmp = M/(R0*R0*pi*t)
		tmp = tmp*x**(-.25)*exp(-(1+x*x)/t)*iv(.25,2*x/t) + floor
		return tmp

dr = 0
if( Lambda == 1 ):
	dr = (rMax - rMin)/(N-1.0)
else:
	dr = (rMax - rMin)*(Lambda-1.0)/(pow(Lambda,N-1)-1.0)

r = np.zeros(N)

for i in range(N):
	if( Lambda == 1 ):
		r[i] = rMin + i*dr
	else:
		if( i == 0 ):
			r[i] = rMin
		else:
			r[i] = rMin + dr*( pow(Lambda,i) - 1.0 )/(Lambda-1.0)

for i in range(64):	# FIXME
	t = .008 + .008*i
	sigma = S(r/R0,t)

	plt.plot(r,sigma*pi*R0**2/M)

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
		f.write(str(r[i]) + '\t' + str(sigma[i]) + '\n')
	f.close()

#
#plt.xlabel('$r/R_0$')
#plt.ylabel('$\pi \Sigma R_0^2/m$')
#plt.title('Surface Density')
#
#plt.legend( ('$\\tau = .008$', '$\\tau = .032$', '$\\tau = .128$', '$\\tau = .512$'),
#		loc='upper right')
#
#plt.savefig('analyticAlphaDisk')
