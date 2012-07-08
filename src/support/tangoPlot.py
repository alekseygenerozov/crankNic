import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from math import pow

mpl.rcParams['figure.subplot.hspace'] = 0.0
mpl.rcParams['figure.subplot.wspace'] = 0.0

#
#		NORMALIZE
#			Assumes 1D np array
#
def normalize(x):
	return x/np.max(np.abs(x))

#
#		strFN
#			Converts # to string with prefix of zeros
#				Assumes we're going from 000 to 999
#
def strFN(n):

	n = int(n)
	fstr = ""
	if( n < 10 ):
		fstr += "00" + str(n)
	elif( n < 100 ):
		fstr += "0" + str(n)
	else:
		fstr += str(n)
	return fstr

#
#		N2FNAME
#			Converts number to filename
#				Assumes data resides in folder outputFiles
#
def n2fName(n):
	if( n < 0 ):
		return "ERROR.dat"
	return 'outputFiles/T' + strFN(n) + '.dat'

def n2IName(n):
	if( n < 0 ):
		return "ERROR"
	return "images/T" + strFN(n)

#
#		FAKE LOG
#			Allows you to take logarithm of negative
#			numbers, which is totally mathematically
#			rigorous, right!?!
#
#			Assumes 1D np array, uses log base-10
#
def fakeLog(x):
	return np.sign(x)*np.log10(np.abs(x))

#
#		GEN GRID
#			Makes the logarithmic mesh, given param file
#
def genGrid(gridSpecs="params.out"):

	# grab parameters
	lMin=0.0
	lMax=0.0
	Lambda=0.0
	N=0

	inType = type(gridSpecs)

	if( inType is str ):
		params = readParams(gridSpecs)
		lMin = params['lMin']
		lMax = params['lMax']
		Lambda = params['lambda']
		N = int(params['N'])
	elif( inType is list ):
		lMin = gridSpecs[0]
		lMax = gridSpecs[1]
		Lambda = gridSpecs[2]
		N = int(gridSpecs[3])
	elif( inType is dict ):
		lMin = gridSpecs['lMin']
		lMax = gridSpecs['lMax']
		Lambda = gridSpecs['lambda']
		N = int(gridSpecs['N'])
	else:
		print "ERROR IN genGrid --- Improper input"
		return False

	dl = 0
	if( Lambda == 1.0 ):
		dl = (lMax - lMin)/(N-1.0)
	else:
		dl = (lMax - lMin)*(Lambda-1.0)/(pow(Lambda,N-1)-1.0)

	l = np.zeros(N)

	for i in range(N):
		if( Lambda == 1.0 ):
			l[i] = lMin + i*dl
		else:
			if( i == 0 ):
				l[i] = lMin
			else:
				l[i] = lMin + dl*( pow(Lambda,i) - 1.0 )/(Lambda-1.0)
	return l

#
#		READ PARAMS
#			Munier Salem, July 2012
#
#			Opens file params.out and grabs all the
#			content, passing it back in a dictionary
#			whose keys are strings and whose values
#			are all FLOATS
#			
#					Assumes format: 
#					" var = #			// comment "
#
def readParams(fName='params.out'):
	f = open(fName,'r')
	fparams = f.read()
	f.close()
	fparams = fparams.split('\n');
	Np = len(fparams)
	
	
	params = {}
	for i in range(Np):
		s = fparams[i]
		s = s.split('=')
		var = s[0].strip()
		if( len(s) > 1 ):
			num = s[1].split('//')
			num = float(num[0].strip())
	
			params[var] = num

	return params

#
#		READ DATA FILE
#			Munier Salem, July 2012
#
#			Opens data file and reads in all rows
#			and columns into a multi-dim numpy array,
#			passing that back to the user.
#
#				Assumes #s separated by SINGLE TABs
#
def readDataFile(n):

	hdrLines = 3

	try:
		fName = n2fName(n)
		f = open(fName,'r')
	except IOError as e:
		print "ERROR IN GRAB DATA ... File " + fName + " does not exist"
		return False

	dataString = f.read()
	f.close()

	dataLines = dataString.split('\n');
	nRows = len(dataLines) - 1
	nCols = len(dataLines[hdrLines].split('\t'))
	data = np.zeros((nRows,nCols))

	for i in range(hdrLines,nRows):
		line = dataLines[i].split('\t')
		for j in range(nCols):
			data[i][j] = float( line[j] )

	return data


#
#		PLOT DATA STD
#
#			Plots triple-paned view of F_J and Tidal 
#			Torque Profile, takes param dict and 
#			either data file name, data file # or data 
#			itself.
#
def plotDataStd(params,n,imType="png"):

	if( (type(n) is np.ndarray) ):
		data = n
	else:
		data = readDataFile(n)

	if( type(data) is bool and not data ):
		return False

	plt.clf()

	l   = data[:,0]
	FJ  = data[:,1]
	trk = data[:,2]

	# UPPER LEFT: F_J Plotted over whole region
	ax1 = plt.subplot2grid((4,4), (0,0), rowspan=3,colspan=3)
	ax1.loglog(l,FJ,'b-')
	plt.ylabel('$F_J$')
	ax1.axis((l.min(),l.max(),1E-1,1E4))
	plt.setp( ax1.get_xticklabels(), visible=False)

	# LOWER LEFT: Torque profile over whole region
	ax2 = plt.subplot2grid((4,4), (3,0), colspan=3,sharex=ax1)
	ax2.semilogx(l,trk,'b-')
	plt.ylabel('$\\Lambda$')
	plt.xlabel('$l$')
	ax2.axis((l.min(),l.max(),trk.min(),trk.max()))

	# RIGHT: Zoom in on region of secondary, both plotted
	subMin = 8.0
	subMax = 12.0
	ax3 = plt.subplot2grid((4,4),(0,3),rowspan=4)
	ax3.plot(l,normalize(trk),'k.-',linewidth=.5,markersize=2)
	ax3.axis((subMin,subMax,-1.0,1.0))
	
	lsub = l[ l > subMin]
	lsub = lsub[ lsub < subMax]
	Fsub = FJ[ l > subMin ]
	Fsub = Fsub[ lsub < subMax ] 
	Fsub = Fsub/100 - 0.8
	ax3.plot(lsub,Fsub,'b-',linewidth=1.5)

	plt.xlabel('l')
	plt.ylabel('$\\Lambda$')
	plt.setp( ax3.get_yticklabels(), visible=False)

	fName = n2IName(n)
	if( imType == "ps" ):
		plt.savefig(fName + '.ps')
	else:
		plt.savefig(fName)

	return True

#
#			PLOT ALL STD
#
#			A one-liner call that reads in all output
#			files (with an optional param for going 
#			by twos, threes, tens, etc ... ) that pulls
#			sim info from params.out (or another file)
#			and produces a std 3-pane plot of F_J and trk
#
def plotAllStd(skip=1,pFile="params.out"):

	params = readParams(pFile)
	
	n = 0
	while(plotDataStd(params,n)): 
		n += skip

#
#		OPT RES
#
#			Given region of interest, available pts, and target spot
#			will determine a logarithmic stretch factor to optimize
#			resolution at that spot
#
def optRes(lMin,lMax,l_a,N):

  dLam = .00001
  Lambda = np.arange(1.0+dLam,1.03,dLam)

  def dla(ll):
    return (ll-1.0)*(l_a-lMin+(lMax-lMin)/(pow(ll,N)-1.0))

  lambdaMin = Lambda[0]
  dlMin = dla(lambdaMin)
  for ll in Lambda:
    dl = dla(ll)
    if( dl < dlMin ):
      lambdaMin = ll
      dlMin = dl

  return ( lambdaMin , dlMin )

