# Cut Norm Approximation #

## What is the cut norm approximation? ##
The code contained herein uses SDPA to calculate a semidefinite approximation to the cut norm and then applied a rounding technique of [1]. We use the so-called infinity-to-one norm to approximate the cut norm, which is given by the quadratic integer program contained in equation (1) in [1]. The semidefinite relaxation (SDR) is given in equation (2) of [1]. We then apply the algorithm in Section 5.1 of [1] to round the SDR to a feasible solution of equation (1).

As we want to use this as a norm: ||A-B||, two inputs are required. If you want the cut norm of just a single matrix, A, make B=0 the zero matrix.

## Installation ##
The following will allow you to install SDPA on a fresh install of Ubuntu (tested on Ubuntu 14.04.3 LTS).
```bash
#dependencies
sudo apt-get update
sudo apt-get install g++ patch
sudo apt-get install gfortran
sudo apt-get install liblapack-dev liblapack-doc-man liblapack-doc liblapack-pic liblapack3 liblapack-test liblapack3gf liblapacke liblapacke-dev
sudo apt-get install libatlas-base-dev

#SDPA
cd /Desktop
mkdir SDPA
cd SDPA
wget -O sdpa_7.3.8.tar.gz http://sourceforge.net/projects/sdpa/files/sdpa/sdpa_7.3.8.tar.gz/download
tar xvfz sdpa_7.3.8.tar.gz
cd sdpa_7.3.8
./configure
make
sudo cp sdpa /usr/local/bin

#numpy and scipy
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
```

The script(s) contained in ``src`` are written in python and should run upon installation.

## Usage ##
The script ``CutnormApprox.py`` is run via
```
./CutnormApprox.py -m <MaxEntropyMatrix.csv> -s <SampleAveMatrix.csv> -o <Output.txt> -e <ExecutableForSDPA>'
```

For example, using the included test data, upon running
```
./CutnormApprox.py -m ../test/A.csv -s ../test/B.csv -o ../test/out.txt
```

The file ``out.txt`` should match the included file ``/test/AB_out.txt``.

## Output format ##
The output format consists of three fields: the approximation value (``#Approximation``), the pertinent part of the output of SDPA (``#Ymat``) and the rounded variable values (``#ui's (rows)`` and ``#vj's (columns)``).

The ``uis`` and ``vjs`` correspond to the ``x_i`` and ``y_i`` in equation (1) of [1].


## Citations ##
[1] Alon, Noga, and Assaf Naor. "Approximating the cut-norm via Grothendieck's inequality." SIAM Journal on Computing 35.4 (2006): 787-803.