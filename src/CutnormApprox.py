#!/usr/bin/env python

import os, sys, shutil, subprocess, getopt, math, copy
import numpy as np
from cStringIO import StringIO
import scipy
import scipy.linalg

# file_path_max_ent='MedBilayerG1_Z.csv'
# file_path_sample='MedBilayerG1_warswap.out_E_out.csv'
# SDPA_exec_path = "/home/david/Desktop/SDPA/sdpa-7.3.8/./sdpa"
# file_path_output='test.txt'


def main(argv):
	file_path_max_ent = None
	file_path_sample = None
	file_path_output = None
	SDPA_exec_path = "sdpa"
	num_threads = 1
	try:
		opts, args = getopt.getopt(argv, "m:s:o:e:t:", ["MaxEntropyMatrix=", "SampleAveMatrix=", "Output=","SDPAExecutable=","NumThreads="])
	except getopt.GetoptError:
		print 'Call using: python CutnormApprox.py -m <MaxEntropyMatrix.csv> -s <SampleAveMatrix.csv> -o <Output.txt> -e <ExecutableForSDPA> -t <NumThreads>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Call using: python CutnormApprox.py -m <MaxEntropyMatrix.csv> -s <SampleAveMatrix.csv> -o <Output.txt> -e <ExecutableForSDPA> -t <NumThreads>'
			sys.exit(2)
		elif opt in ("-m", "--MaxEntropyMatrix"):
			file_path_max_ent = arg
		elif opt in ("-s", "--SampleAveMatrix"):
			file_path_sample = arg
		elif opt in ("-o", "--Output"):
			file_path_output = arg
		elif opt in ("-e", "--SDPAExecutable"):
			SDPA_exec_path = arg
		elif opt in ("-t", "--NumThreads"):
			num_threads = arg
	#Run the main program
	calc_cutnorm(file_path_max_ent, file_path_sample, file_path_output,SDPA_exec_path,num_threads)

def calc_cutnorm(file_path_max_ent, file_path_sample, file_path_output,SDPA_exec_path,num_threads):
	assert isinstance(file_path_max_ent,basestring), file_path_max_ent
	assert isinstance(file_path_sample,basestring), file_path_sample
	assert isinstance(file_path_output,basestring), file_path_output
	
	#Read in files
	Z = np.genfromtxt(file_path_max_ent,delimiter=',')
	S = np.genfromtxt(file_path_sample,delimiter=',')
	num_rows = np.shape(Z)[0]
	num_columns = np.shape(Z)[1]
	num_var = num_rows*num_columns
	if np.shape(Z) != np.shape(S):
		print("Error: shape of matrices are not the same")
		sys.exit(2)
	#Make input file for SDPA
	base_path = os.path.dirname(os.path.abspath(file_path_output))
	SDPA_input_file = os.path.join(base_path,os.path.splitext(os.path.basename(file_path_output))[0]+'_SDPAinput.dat-s')
	SDPA_output_file = os.path.join(base_path,os.path.splitext(os.path.basename(file_path_output))[0]+'_SDPAoutput.txt')
	fid = open(SDPA_input_file,'w')
	D = Z-S
	make_SDPA_input(fid,D)
	fid.close()
	
	#Run SDPA
	cmd = SDPA_exec_path + " " + SDPA_input_file + " " + SDPA_output_file + " " + "-numThreads " + str(num_threads)
	res = subprocess.check_output(cmd, shell = True)
	
	#Parse SDPA output
	fid = open(SDPA_output_file,"r")
	SDPA_output = fid.read()
	fid.close()
	y_mat_text = find_between(SDPA_output,"yMat = \n{\n{ ", "   }\n}\n    main").replace("{","").replace("}","").replace(" ","").replace(",\n","\n") 	#Yes, the SDPA output format (if you can call it that) is a BEAR to parse!
	Y = np.genfromtxt(StringIO(y_mat_text),delimiter=',')
	#eigen_values,eigen_vectors = np.linalg.eig(Y) #could probably use eigh
	#eigen_values,eigen_vectors = scipy.linalg.eigh(Y,eigvals=(num_var-2,num_var-1)) #assumes Y symmetric, use this and a loop if it's too inneficient to compute ALL the eigenvalues, just loop through getting the top k until the are <=0
	eigen_values,eigen_vectors = scipy.linalg.eigh(Y)
	idx = eigen_values.argsort()[::-1]
	eigen_values = eigen_values[idx]
	all_eigen_values = eigen_values[:]
	eigen_vectors = eigen_vectors[:,idx]
	for positive_index in range(len(eigen_values)):
		if eigen_values[positive_index]<=0:
			break

	eigen_values = eigen_values[0:positive_index]
	eigen_vectors = eigen_vectors[:,0:positive_index]
	soln = np.zeros(np.shape(eigen_vectors))
	for i in range(len(eigen_values)):
		soln[:,i] = np.sqrt(eigen_values[i])*eigen_vectors[:,i] #coords/variables are row-wise. eg: [u1;u2;v1;v2]

	#Perform Alon & Noar rounding
	approx_opt = 0
	uis_opt = list()
	vjs_opt = list()
	for dummy in range(1,1000): #I'll use the sample average, can make this a parameter later
		G = np.zeros(positive_index)
		for i in range(len(G)):
			G[i] = np.random.normal()
			
		rounded_soln = np.sign(soln.dot(G))
		uis = rounded_soln[0:num_rows]
		vjs = rounded_soln[num_rows:num_rows+num_columns]
		#Calculate approximation
		approx = 0
		for i in range(num_rows):
			for j in range(num_columns):
				approx = approx + D[i,j]*uis[i]*vjs[j]

		approx = np.abs(approx)
		if approx > approx_opt:
			approx_opt = copy.copy(approx)
			uis_opt = copy.copy(uis)
			vjs_opt = copy.copy(vjs)

	#Save output
	fid = open(file_path_output,'w')
	fid.write("#Approximation\n")
	fid.write("%f\n" % approx_opt)
	fid.write("#Ymat\n")
	for i in range(num_var):
		for j in range(num_var):
			if j == num_var-1:
				fid.write("%f" % Y[i,j])
			else:
				fid.write("%f," % Y[i,j])
		fid.write("\n")

	fid.write("#ui's (rows)\n")
	for i in range(len(uis_opt)):
		if i == len(uis_opt)-1:
			fid.write("%f\n" % uis_opt[i])
		else:
			fid.write("%f," % uis_opt[i])

	fid.write("#vj's (rows)\n")
	for i in range(len(vjs_opt)):
		if i == len(vjs_opt)-1:
			fid.write("%f\n" % vjs_opt[i])
		else:
			fid.write("%f," % vjs_opt[i])

	fid.write("#Eigenvalues\n")
	for i in range(len(all_eigen_values)):
		if i == len(all_eigen_values)-1:
			fid.write("%f\n" % all_eigen_values[i])
		else:
			fid.write("%f," % all_eigen_values[i])


	fid.close()



def make_SDPA_input(fid,D):
	n_rows = np.shape(D)[0]
	n_columns = np.shape(D)[1]
	num_var = n_columns*n_rows
	fid.write("\"Example\"\n")
	fid.write("%d = mDIM\n" % num_var)
	fid.write("1 = nBLOCK\n")
	fid.write("%d = bLOCKsTRUCT\n" % num_var)
	for it in range(n_rows*n_columns):
		fid.write("1 ")
	fid.write("\n")
	Is = list()
	Js = list()
	for i in range(n_rows):
		for j in range(n_columns):
			fid.write("0 1 %d %d %f\n" % (i+1, j+n_rows+1, D[i,j]))
			Is.append(i+1)
			Js.append(j+n_rows+1)
	Is = list(set(Is))
	Js = list(set(Js))
	for i in Is:
		fid.write("%d 1 %d %d 1\n" % (i,i,i))
	for j in Js:
		fid.write("%d 1 %d %d 1\n" % (j,j,j))


def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""




if __name__ == "__main__":
	main(sys.argv[1:])
