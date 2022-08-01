# Randomized-Tensor-Train-research
Authors: Chris Camaño and Roel Van Beeumen 
San Francisco State University and Lawrence Berkeley National Laboratory 

This work was supported in part by the U.S. department of Energy, Office of Science, Office of Workforce Development for Teachers and Students (WDTS) under the Science Undergraduate Laboratory Internship (SULI) Program. 

The following project is active research on randomizing the The density-matrix renormalization group (DMRG) algorithm also known as the approximate linerized scheme (ALS) or modfied approximate linearized scheme (MALS) 

Input data is a randomly generate Hamiltonian Heisenberg matrix with dimensions 2^n x 2^n formed in tensor train format. 

Randomization is included throughout the original algorithm with observable improvements over conventional qr decomposition approaches to right to left orthogonalization of tensor cores. 


Additional work includes an unoptimized implementation of  the sketched rayleigh-ritz form the following paper:
Yuji Nakatsukasa and Joel Tropp in Feburary 2022 preprint:
Fast & Accurate Randomized Algorithms For Linear Systems and Eigenvalue Problems. 

