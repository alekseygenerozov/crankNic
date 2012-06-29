import numpy as np

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

	try:
		fName = ""
		if( type(n) is int ):
			fName += n2fName(n)
		else:
			fName = n
		f = open(fName,'r')
	except IOError as e:
		print "ERROR IN GRAB DATA ... File " + fileName + " does not exist"
		return False

	dataString = f.read()
	f.close()

	dataLines = dataString.split('\n');
	nRows = len(dataLines) - 1
	nCols = len(dataLines[0].split('\t'))
	data = np.zeros((nRows,nCols))

	for i in range(nRows):
		line = dataLines[i].split('\t')
		for j in range(nCols):
			data[i][j] = float( line[j] )

	return data
