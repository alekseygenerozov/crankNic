import numpy as np

#
#		NORMALIZE
#			Assumes 1D np array
def normalize(x):
	return x/np.max(np.abs(x))	

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
def readParams():
	f = open('params.out','r')
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
def readDataFile(fileName):

	try:
		f = open(fileName,'r')
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
