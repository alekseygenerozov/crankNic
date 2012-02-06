import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi, sqrt, log

def XI( x , alpha , beta, gamma, dr ):

	A = 4.0*alpha*sin(x/2)*sin(x/2) + beta*dr*dr*(1.5+gamma)
	B = 0.5*dr*(3.0*alpha-2.0*beta)*sin(x)

	ap1 = A+1.0;
	denom = ap1*ap1 + B*B

	R = 2.0*ap1/denom-1.0
	I = 2.0*B/denom

	return sqrt( R*R + I*I )

def makePlot( a , b, c, dr ):
	N  = 500
	xMin = 0;
	xMax = 2*pi;
	dx = (xMax - xMin)/(N - 1.0);

	x = np.zeros(N)
	s = np.zeros(N)

	for i in range(N):
		x[i] = i * dx + xMin
		s[i] = XI(x[i],a,b,c,dr)

	plt.plot( x , s )


makePlot( 0.1 , 1000.0 , 0.001 , .01 )
makePlot( 0.1 , 1000.0 , 0.01 , .01 )
makePlot( 0.1 , 1000.0 , 0.10 , .01 )
makePlot( 0.1 , 1000.0 , 1.00 , .01 )
makePlot( 0.1 , 1000.0 , 10.00 , .01 )
makePlot( 0.1 , 1000.0 , 100.00 , .01 )
makePlot( 0.1 , 1000.0 , 1000.00 , .01 )
makePlot( 0.1 , 1000.0 , 10000.00 , .01 )

plt.axis([-0.1,7.0,0,1.1]);
plt.xlabel('$k \\delta r$')
plt.ylabel('$| \\xi | $')

plt.legend([
		'$\\gamma = 0.0010$',
		'$\\gamma = 0.010$',
		'$\\gamma = 0.10$',
		'$\\gamma = 1.00$',
		'$\\gamma = 10.00$',
		'$\\gamma = 100.00$',
		'$\\gamma = 1000.00$',
		'$\\gamma = 10000.00$',
		],'upper right')


plt.savefig("xi")
