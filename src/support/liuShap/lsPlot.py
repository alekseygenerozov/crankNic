from tangoPlot import *
import matplotlib as mpl

plotAllStd()
n = 0
params = readParams()
ana = readDataFile("analytic.dat",0)
while(plotSigmaStd(params,n,'png',ana)):
	n += 1

# comparison plot

sim = readDataFile(n-1)
(r ,s ) = convert(params,sim[:,0],sim[:,1])
(ra,sa) = convert(params,ana[:,0],ana[:,1])

import matplotlib.pyplot as plt

plt.loglog(ra,sa,'k--',linewidth=1)
plt.loglog(r ,s ,'k-',linewidth=2)
plt.axis((r.min(),r.max(),1E-6,1))
plt.savefig('comparison')
