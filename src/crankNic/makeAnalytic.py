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

import matplotlib.pyplot as plt
from scipy import arange, pi, sqrt, exp
from scipy.special import iv
from matplotlib import rc

R0 = 1.0 		# initial radius of delta fcn
M  = 1.0		# intiial mass
pi = 3.14159

N = 500
rMin = 0.2
rMax = 2.0

# we read in from param file ...
f = open('params.dat','r')
params = f.read()
f.close()
params = params.split('\n');
N = len(params)

for i in range(N):
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

print "R0 = " + str(R0)
print "rMin = " + str(rMin)
print "rMax = " + str(rMax)
print "N = " + str(N)


#
# ---- Our analytic expression
#
def S(x,t):
    tmp = M/(R0*R0*pi*t)
    return tmp*x**(-.25)*exp(-(1+x*x)/t)*iv(.25,2*x/t)

dr = (rMax - rMin)/(N+1.0)
r = arange(rMin,rMax,dr)

for i in range(4):	# FIXME
	t = .008*4**i
	sigma = S(r/R0,t)
	plt.plot(r,sigma*pi*R0**2/M)
	
	f = open('analytic_T' + str(i) + '.dat','w')
	for j in range(N):
		f.write(str(r[j]) + '\t' + str(sigma[j]) + '\n')
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
