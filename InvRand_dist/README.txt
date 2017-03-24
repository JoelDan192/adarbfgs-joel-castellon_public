InvRand
---------------------------------------------------------------------------

1. Introduction
===============

This is a suite of randomized methods for inverting positive definite matrices implemented in MATLAB, based on methods described in

[1]	Robert M. Gower and Peter Richtarik,
    Randomized Quasi-Newton Updates are Linearly Convergent Matrix Inversion Algorithms, 
    arXiv:1602.01768

The methods implemented are

AdaRBFGS        (Adaptive Randomized Block BFGS)
MR              (Minimal Residual method with self-conditioning)
Newton-Schulz   (The Newton-Schulz method)              

2. Installation and Setup
=========================

Start Matlab and make sure that the working directory is set to the
main directory of the present package.  At the MATLAB prompt, run

  >> setup_InvRand

To test if the installation and setup have been 
completed successfully please run in the MATLAB prompt:

  >> demo_InvRand
  >> demo_LIBSVM
  >> demo_UF

3. Repeat tests in paper [1]
============================

WARNING: The following experiments are CPU and data intensive!

To run the tests carried out in the paper [1] do the following.
To repeat all random matrix tests, at the MATLAB prompt run

  >> test_all_User() 

To repeat all the UF sparse matrix tests, at the MATLAB prompt run

  >> problems =  [1401, 440, 760, 40, 868, 937, 1288];
  >> test_all_UF(problems)                     % Tests a sequence of UF Matrices [3]

NOTE: this will download approx. a GB of data to your local hard drive!
To repeat all the LIBSVM tests, run the commands

  >> get_LIBSVM_data
  >> problems =  {'aloi', 'protein','gisette_scale', 'real-sim' }
  >> prefix = './tests/LIBSVM/data/';
  >> test_all_LIBSVM(problems,prefix)           % Tests a sequence of ridge regression matrices using data from LIBSVM [2]

NOTE: the script 'get_LIBSVM_data' will download a total of 380MB to your local hard drive.

4. License
==========

 InvRand Copyright (C) 2016, Robert Gower 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.


4. The Authors
==============

If you have any bug reports or comments, please feel free to email me

  Robert Gower <r.m.gower@sms.ed.ac.uk>
  Robert Gower <gowerrobert@gmail.com>

Robert Gower
17 Febuary 2016


5. References
=============

[1]	Robert M. Gower and Peter Richtarik,
    Randomized Quasi-Newton Updates are Linearly Convergent Matrix Inversion Algorithms, 
    arXiv:1602.01768

[2] Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines.
    ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011. 
    Software available at http://www.csie.ntu.edu.tw/~cjlin/libsvm 

[3] Timothy A. Davis and Yifan Hu. 2011. The university of Florida sparse matrix collection. 
    ACM Trans. Math. Softw. 38, 1, Article 1 (December 2011), 25 pages.
    https://www.cise.ufl.edu/research/sparse/matrices/
