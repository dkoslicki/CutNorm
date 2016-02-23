# Cut Norm Approximation #

## What is the cut norm approximation? ##
The code contained herein uses CSDP to calculate a semidefinite approximation to the cut norm and then applied a rounding technique of [1]. We use the so-called infinity-to-one norm to approximate the cut norm, which is given by the quadratic integer program contained in equation (1) in [1]. The semidefinite relaxation (SDR) is given in equation (2) of [1]. We then apply the algorithm in Section 5.1 of [1] to round the SDR to a feasible solution of equation (1).

As we want to use this as a norm: ||A-B||, two inputs are required. If you want the cut norm of just a single matrix, A, make B=0 the zero matrix.

## Installation ##
The following will allow you to install CSDP on a fresh install of Ubuntu (tested on Ubuntu 14.04.3 LTS).
```bash
#Get CSDP binary
wget http://www.coin-or.org/download/binary/Csdp/csdp6.1.0linuxp4.tgz
tar -xzvf csdp6.1.0linuxp4.tgz

#The binary executable is now located at: csdp6.1.0linuxp4/bin/./csdp

#numpy and scipy
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
```

The script(s) contained in ``src`` are written in python and should run upon installation.

## Usage ##
The script ``CutnormApprox.py`` is run via
```
python CutnormApprox.py -m <MaxEntropyMatrix.csv> -s <SampleAveMatrix.csv> -o <Output.txt> -e <ExecutableForCSDP>'
```

For example, using the included test data, upon running
```
python CutnormApprox.py -m ../test/A.csv -s ../test/B.csv -o ../test/out.txt
```

The file ``out.txt`` should match the included file ``/test/AB_out.txt``.

The script ``MaxEntMatrix.py`` is run via
```bash
python MaxEntMatrix.py -c <ColumnDegrees.csv> -r <RowDegrees.csv> -o <Output.csv>
```

The file ``ColumnDegrees.csv`` is a csv file with a single line representing the degrees of the columns.

The file ``RowDegrees.csv`` is a csv file with a single line representing the degrees of the rows.

The file ``Output.csv`` is a csv file containing the maximum entropy matrix Z from Barvinok 2009.

## Output format ##
The output format consists of five fields:
1. The approximation value of the infinity to one norm (``#Approximation of infinity-to-one norm``)
2. The ``#Interval of cut norm approximation``
3. The rounded solutions for the ``#uis``
4. The rounded solutions for the ``#vjs``
5. The ``#Eigenvalues``. Note that if there are only two distinct eigenvalues, then the lower bound of ``#Interval of cut norm approximation`` is exactly equal to the cut norm.

The ``uis`` and ``vjs`` correspond to the ``x_i`` and ``y_i`` in equation (1) of [1].

## CSDP details ##
The CSDP algorithm [2] uses a predictor–corrector variant of the primal–dual method of Helmberg, Rendl, Vanderbei, and Wolkowicz [3] and Kojima, Shindoh, and Hara [4]. The sparse input format is utilized, and input/output files can be found with extensions ``_CSDPinput.dat-s`` and ``_CSDPoutput.txt``.


## Citations ##
[1] Alon, N., and Naor, A. "Approximating the cut-norm via Grothendieck's inequality." SIAM Journal on Computing 35.4 (2006): 787-803.

[2] Borchers, Brian. "CSDP, AC library for semidefinite programming." Optimization methods and Software 11, no. 1-4 (1999): 613-623.

[3] Helmberg, Christoph, Franz Rendl, Robert J. Vanderbei, and Henry Wolkowicz. "An interior-point method for semidefinite programming." SIAM Journal on Optimization 6, no. 2 (1996): 342-361.

[4] Kojima, Masakazu, Susumu Shindoh, and Shinji Hara. "Interior-point methods for the monotone semidefinite linear complementarity problem in symmetric matrices." SIAM Journal on Optimization 7, no. 1 (1997): 86-125.