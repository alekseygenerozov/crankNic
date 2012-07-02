import numpy as nmp
import matplotlib.pyplot as plt
from tangoPlot import *

def scaledPlot(d,params,chi):
	tau = params['tEnd']
	np  = params['np']
	nd  = params['nd']

	nn = -1.0/( nd + np - 2.0)
	denom = pow(tau,nn)

	x = d[:,0]/denom
	f = d[:,1]/denom

	plt.plot(x,f)
	genSelfSim(x,f,params,chi)

def genSelfSim(x,f,params,chi):

	nd   = params['nd']
	np   = params['np']
	nn = -1.0/(nd+np-2.0)

	N = x.size

	dx   = x[1]-x[0]
	dx2 = dx*dx

	# initialize FJ
	F = nmp.zeros(N)

	# copy into it
	for i in range(len(f)):
		F[i] = f[i]

	dt = .1*dx2
	tol = 1.0E-5
	maxDiff = 10.0*tol

	min_iter = 10
	count = 0
	w = 1.5

	while( maxDiff > tol or count < min_iter ):
		count = count + 1
		maxDiff = 0.0

		for j in range(1,int(N)-2):

			FOld = F[j]

			F_ll = (F[j+1] - 2.0*F[j] + F[j-1])/dx2
			F_l = (F[j+1]-F[j-1])/(2.0*dx)

			Fd = pow(F[j],nd)
			xp = pow(x[j],np)

			F[j] += w*dt*( F_ll*Fd*xp - nn*(1.0-np)*(F[j]-x[j]*F_l))

			diff = abs(FOld - F[j])
			if( diff > maxDiff ):
				maxDiff = diff

		F[0] = F[1] - chi*dx    # F'(xMin) = chi
		F[N-1] = F[N-2] + dx          # free ( F' = 0 )

	plt.plot(x[::6],F[::6],'k.',markersize=6)





params = readParams()
scaledPlot(readDataFile("savedOutputs/chi_p0.dat"),params,.0)
scaledPlot(readDataFile("savedOutputs/chi_p2.dat"),params,.2)
scaledPlot(readDataFile("savedOutputs/chi_p5.dat"),params,.5)
scaledPlot(readDataFile("savedOutputs/chi_1p.dat"),params,1.)

plt.axis((0.0,3.0,0.0,3.0))
plt.savefig('scaled')
