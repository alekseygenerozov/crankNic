import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi, sqrt, log

def XI( x , alpha , dr ):

	A = 4.0*alpha*sin(x/2)*sin(x/2)
	B = 3.0/2.0*dr*alpha*sin(x)

	ap1 = A+1.0;
	denom = ap1*ap1 + B*B

	R = 2.0*ap1/denom-1.0
	I = 2.0*B/denom

	return sqrt( R*R - I*I )
#sqrt(R*R - I*I)

def makePlot( a , dr ):
	N  = 200
	dx = 4*pi/( N - 1.0)

	x = np.zeros(N)
	s = np.zeros(N)

	for i in range(N):
		x[i] = i * dx
		s[i] = XI(x[i],a,dr)

	plt.plot( x , s )


makePlot( 0.1 , .5 )
makePlot( 0.1 , .3 )
makePlot( 0.1 , .1 )
makePlot( 0.1 , .03 )
makePlot( 0.1 , .01 )

plt.xlabel('$k \\delta r$')
plt.ylabel('$| \\xi | $')

plt.legend([
		'$\\delta r = 0.5$',
		'$\\delta r = 0.3$',
		'$\\delta r = 0.1$',
		'$\\delta r = 0.03$',
		'$\\delta r = 0.01$',
#		'$\\delta r = 10.0$',
		],'upper right')


plt.savefig("xi")
