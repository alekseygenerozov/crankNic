
import matplotlib
import matplotlib.pyplot as plt
import numpy
from numpy import log
from tangoPlot import *

matplotlib.rc('lines',markersize=4,linewidth=1);


def fakeLog(x):
	return x #log(1+x)

#
# ------ We Read In Data
#

def grabData(fName,pType):
	
	hdrLines = 0
	if( pType == 0 ):
		hdrLines = 3

	data = readDataFile(fName,hdrLines)

	if( type(data) is bool and not data ):
		return False

	if(pType == 0):
		plt.plot( data[:,0] , data[:,1]*3.14159) 
	else:
		plt.plot( data[:,0] , data[:,1]*3.14159 , 'k.', markersize=4)

	return data

sim1 = grabData("outputFiles/T000.dat",0)
sim2 = grabData("outputFiles/T003.dat",0)
sim3 = grabData("outputFiles/T015.dat",0)
sim4 = grabData("outputFiles/T063.dat",0)

exact1 = grabData("analytic/T000.dat",1)
exact2 = grabData("analytic/T003.dat",1)
exact3 = grabData("analytic/T015.dat",1)
exact4 = grabData("analytic/T063.dat",1)


plt.axis((0.0,2.0,0.0,4.0))
plt.xlabel('$l/l_0$')
plt.ylabel('$F_J$')
plt.title('$F_J(l)$ for Consant $\\nu$ Disk')
plt.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$','analytic'],'upper right')
plt.savefig("numerical")


#
# Residuals:
#
r = sim1[ : , 0 ]
resid1 = (sim1[ : , 1 ] - exact1[ : , 1 ])
resid2 = (sim2[ : , 1 ] - exact2[ : , 1 ])
resid3 = (sim3[ : , 1 ] - exact3[ : , 1 ])
resid4 = (sim4[ : , 1 ] - exact4[ : , 1 ])


plt.clf()
plt.plot( r , fakeLog(resid1))
plt.plot( r , fakeLog(resid2))
plt.plot( r , fakeLog(resid3))
plt.plot( r , fakeLog(resid4))
plt.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$'],'lower right')
plt.xlabel('$l/l_0$')
plt.ylabel('$residuals$')
plt.title('Residuals')
plt.axis((0,2.0,-.0001,.0001))

plt.savefig("resids")

"""
#
#  mDot:
# 
plt.clf()
plt.plot(r,sim1[ :  , 3 ])
plt.plot(r,sim2[ :  , 3 ])
plt.plot(r,sim3[ :  , 3 ])
plt.plot(r,sim4[ :  , 3 ])
plt.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$'],'lower right')
plt.xlabel('$r/R_0$')
plt.ylabel('$\\dot{M} = 3\\pi\\nu\\Sigma$')
plt.title('Mass Flow')
plt.savefig("mDot")
"""
