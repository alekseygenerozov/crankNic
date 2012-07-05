import matplotlib as mpl
import numpy as nmp
import matplotlib.pyplot as plt
from tangoPlot import *

mpl.rcParams['figure.subplot.wspace'] = 0.0

def scaledPlot(d,params,chi,ax,analytic=True):
	tau = params['tEnd']
	np  = params['np']
	nd  = params['nd']
	N = int(params['N'])
	nn = -1.0/( nd + np - 2.0)
	denom = pow(tau,nn)

	x = d[:,0]/denom
	f = d[:,1]/denom

	ax.plot(x,f)
	ax.text(0.1,f[N/6],"$\\chi = $ " + str(chi),fontsize=11)
	if(analytic):
		genSelfSim(x,f,params,chi,ax)

def genSelfSim(x,f,params,chi,ax):

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

			Fd = pow(abs(F[j]),nd)
			xp = pow(x[j],np)

			F[j] += w*dt*( F_ll*Fd*xp - nn*(1.0-np)*(F[j]-x[j]*F_l))

			diff = abs(FOld - F[j])
			if( diff > maxDiff ):
				maxDiff = diff

		F[0] = F[1] - chi*dx    # F'(xMin) = chi
		F[N-1] = F[N-2] + dx          # free ( F' = 0 )

	ax.plot(x[::8],F[::8],'k.',markersize=5)


ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)

params = readParams("savedOutputs/params_2.out")
scaledPlot(readDataFile("savedOutputs/chi_p0_2.dat"),params,.0,ax2)
scaledPlot(readDataFile("savedOutputs/chi_p2_2.dat"),params,.2,ax2)
scaledPlot(readDataFile("savedOutputs/chi_p5_2.dat"),params,.5,ax2)
scaledPlot(readDataFile("savedOutputs/chi_1p_2.dat"),params,1.,ax2)

params = readParams("savedOutputs/params_1.out")
scaledPlot(readDataFile("savedOutputs/chi_p0_1.dat"),params,.0,ax1, False)
scaledPlot(readDataFile("savedOutputs/chi_p2_1.dat"),params,.2,ax1, False)
scaledPlot(readDataFile("savedOutputs/chi_p5_1.dat"),params,.5,ax1, False)
scaledPlot(readDataFile("savedOutputs/chi_1p_1.dat"),params,1.,ax1, False)

plt.subplot(1,2,1)
ax1.axis((0.0,3.0,0.0,3.0))
ax1.text(1.5,0.5,"$n_d = 2/5$\n$n_p = -6/5$",fontsize=16)
plt.xlabel('$\\xi = 1/\\tau^{5/14}$')
plt.ylabel('$f(\\xi)$')

plt.subplot(1,2,2)
ax2.axis((0.0,3.0,0.0,3.0))
plt.xlabel('$\\xi = 1/\\tau^{2/5}$')
ax2.text(1.5,0.5,"$n_d = 3/10$\n$n_p = -4/5$",fontsize=16)
plt.setp( ax2.get_yticklabels(), visible=False)

plt.savefig('images/numerical')
plt.savefig('images/numerical.ps')
