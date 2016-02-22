import os, sys, shutil, subprocess, getopt, math, copy
import numpy as np
from cStringIO import StringIO
import scipy
import scipy.linalg
import scipy.optimize

# file_path_c='test_c.csv'
# file_path_r='test_r.csv'
# file_path_output='out_Z.csv'


def main(argv):
	file_path_c = None
	file_path_r = None
	file_path_output = None
	try:
		opts, args = getopt.getopt(argv, "c:r:o:", ["ColumnDegrees=", "RowDegrees=", "Output="])
	except getopt.GetoptError:
		print 'Call using: python MaxEntMatrix.py -c <ColumnDegrees.csv> -r <RowDegrees.csv> -o <Output.csv>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Call using: python MaxEntMatrix.py -c <ColumnDegrees.csv> -r <RowDegrees.csv> -o <Output.csv>'
			sys.exit(2)
		elif opt in ("-c", "--ColumnDegrees"):
			file_path_c = arg
		elif opt in ("-r", "--RowDegrees"):
			file_path_r = arg
		elif opt in ("-o", "--Output"):
			file_path_output = arg
	#Run the main program
	calc_max_ent(file_path_c, file_path_r, file_path_output)

def G(x,r,c): #This is the G function on page 3 of Barvinok 2009
	m = len(r)
	n = len(c)
	s = x[0:m]
	t = x[m:]
	res = -np.sum(r*s) - np.sum(c*t) + np.sum(np.log(1+np.exp(t[:,np.newaxis]+s))) #-\sum_i r_i*s_i - \sum_i c_i*t_i + sum_{i,j} \log(1+e^{s_i+t_j})
	return res


def calc_max_ent(file_path_c, file_path_r, file_path_output):
	assert isinstance(file_path_c,basestring), file_path_c
	assert isinstance(file_path_r,basestring), file_path_r
	assert isinstance(file_path_output,basestring), file_path_output
	
	#Read in files
	c_degrees = np.genfromtxt(file_path_c,delimiter=',')
	r_degrees = np.genfromtxt(file_path_r,delimiter=',')
	
	m = len(r_degrees)
	n = len(c_degrees)
	x0 = np.concatenate((r_degrees/np.sum(r_degrees), c_degrees/np.sum(c_degrees)))

	#BFGS
	res = scipy.optimize.minimize(G, x0, args=(r_degrees,c_degrees))
	res = res.x
	
	x = np.exp(res[0:m])
	y = np.exp(res[m:])
	Z = np.zeros((m,n))
	for i in range(m):
		for j in range(n):
			Z[i,j] = x[i]*y[j]/(1+x[i]*y[j])
	
	
	np.savetxt(file_path_output, Z, delimiter=',', fmt='%.10f')


	


if __name__ == "__main__":
	main(sys.argv[1:])
