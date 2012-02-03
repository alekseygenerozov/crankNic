#
#     ANALYZE OUTPUT LEVEL INFORMATION
#
#         Takes the Output Level Info file from
#         enzo and plots the coverage per refined
#         level as a function of time using 
#	  matplotlib
#

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('lines',markersize=4,linewidth=1);


#
# ------ We Read In Data
#

def grabData(fileName,type):

	exists = 1
	try:
		f = open(fileName,'r')
	except IOError as e:
		print "WARNING ...  file '" + fileName + "' does not exist."
		exists = 0

	if( exists == 1 ):
		dataString = f.read()
		f.close()
		dataLines = dataString.split('\n');
		N = len(dataLines) - 1
		data = numpy.zeros((N,2))
	
		for i in range(N):
		  line = dataLines[i].split('\t')
		  data[i][0] = float( line[0] )
		  data[i][1] = float( line[1] )
	
		if(type == 0):
			plt.plot( data[:,0] , data[:,1]*3.14159 ) 
		else:
			plt.plot( data[::12,0] , data[::12,1]*3.14159 , 'k*')

grabData("T008.dat",0)
grabData("T032.dat",0)
grabData("T128.dat",0)
grabData("T512.dat",0)

grabData("analytic_T0.dat",1)
grabData("analytic_T1.dat",1)
grabData("analytic_T2.dat",1)
grabData("analytic_T3.dat",1)


#plt.axis((0.0,2.0,0.0,4.0))
plt.xlabel('$r/R_0$')
plt.ylabel('$\pi \Sigma R_0^2/m$')
plt.title('$\Sigma(r)$ for S&S Disk')
plt.legend(['$\\tau = .008$','$\\tau = .032$','$\\tau = .128$','$\\tau = .512$','analytic'],'upper right')
plt.savefig("numerical")
