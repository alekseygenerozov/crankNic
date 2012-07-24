import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy
from numpy import log
from tangoPlot import *

matplotlib.rc('lines',markersize=4,linewidth=1);

ax1 = plt.subplot(2,1,1)
plt.axis((0.0,2.0,0.0,1.0))
plt.ylabel('$F_J$')
plt.title('$F_J(l)$ for Consant $\\nu$ Disk')
ax1.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$','analytic'],'upper right')

ax2 = plt.subplot(2,1,2)
plt.xlabel('$l/l_0$')
plt.axis((0.0,2.0,-.00005,.00005))

def onePlot(n):
	sim = readDataFile(n)
	ax1.plot(sim[:,0],sim[:,1])

	fName = "analytic/T0"
	if( n > 10 ):
		fName += str(n) + ".dat"
	else:
		fName += "0" + str(n) + ".dat"

	ana = readDataFile(fName,0)
	ax1.plot(ana[:,0],ana[:,1],'k.',markersize=2)

	resid = sim[:,1] - ana[:,1]
	ax2.plot(sim[:,0],resid)

onePlot(0)
onePlot(3)
onePlot(15)
onePlot(63)

plt.savefig("numerical")
plt.clf()

"""
#
# Residuals:
#
def resid(n):
	fName = "analytic/T0"
	if( n > 10 ):
		fName += str(n) + ".dat"
	else:
		fName += "0" + str(n) + ".dat"

	sim = readDataFile(n)
	ana = readDataFile(fName)
	resid = sim[:,1] - ana[:,1]
	r = sim[:,0]

	plt.plot(r,fakeLog(resid))

resid(0)
resid(3)
resid(15)
resid(63)

plt.clf()
plt.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$'],'lower right')
plt.xlabel('$l/l_0$')
plt.ylabel('$residuals$')
plt.title('Residuals')
plt.axis((0,2.0,-.0001,.0001))
plt.savefig("resids")
"""
