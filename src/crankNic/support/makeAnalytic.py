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

l0 = 1.0 		# initial radius of delta fcn
M  = 1.0		# intiial mass
pi = 3.14159

N = 500
lMin = 0.1
lMax = 2.0
Lambda = 1.03
nu = 1.0/12.0

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
	
		if( var == 'l0' ): 	l0 = num
		if( var == 'lMin' ): lMin = num
		if( var == 'lMax' ): lMax = num
		if( var == 'N' ): N = int(num)
		if( var == 'lambda' ): Lambda = num
		if( var == 'nu0' ): nu = num

print "l0 = " + str(l0)
print "lMin = " + str(lMin)
print "lMax = " + str(lMax)
print "N = " + str(N)
print "lambda = " + str(Lambda)
print "nu = " + str(nu)

#
# ---- Our analytic expression
#
def Fj(x,t):
		floor = 1E-4
		tmp = 3.0*M*M*M*nu/(l0*l0*l0*t)
		tmp = tmp*x**(.25)*exp(-(1+x*x)/t)*iv(.25,2*x/t) + floor
		return tmp

dl = 0
if( Lambda == 1 ):
	dl = (lMax - lMin)/(N-1.0)
else:
	dl = (lMax - lMin)*(Lambda-1.0)/(pow(Lambda,N-1)-1.0)

l = np.zeros(N)

for i in range(N):
	if( Lambda == 1 ):
		l[i] = lMin + i*dl
	else:
		if( i == 0 ):
			l[i] = lMin
		else:
			l[i] = lMin + dl*( pow(Lambda,i) - 1.0 )/(Lambda-1.0)

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

#
#plt.xlabel('$r/R_0$')
#plt.ylabel('$\pi \Sigma R_0^2/m$')
#plt.title('Surface Density')
#
#plt.legend( ('$\\tau = .008$', '$\\tau = .032$', '$\\tau = .128$', '$\\tau = .512$'),
#		loc='upper right')
#
#plt.savefig('analyticAlphaDisk')
