# Randomized-Tensor-Train-research
Authors: Chris Camaño and Roel Van Beeumen 
San Francisco State University and Lawrence Berkeley National Laboratory 



The following project is split into two parts of active research
1. Randomizing the The density-matrix renormalization group (DMRG) algorithm also known as the approximate linerized scheme (ALS) or modfied approximate linearized scheme (MALS) 

2. Implementing the Sketched Rayleigh Ritz method for the computation of eigenpairs in sparse matrices. 




## 1. TT-DMRG 

Input data is a randomly generate Hamiltonian Heisenberg matrix with dimensions 2^n x 2^n formed in tensor train format. 

Randomization is included throughout the original algorithm with observable improvements over conventional qr decomposition approaches to right to left orthogonalization of tensor cores. 

## 2. Sketched Rayleigh-Ritz 
Additional work includes an unoptimized implementation of  the sketched rayleigh-ritz form the following paper:
Yuji Nakatsukasa and Joel Tropp in Feburary 2022 preprint:Fast & Accurate Randomized Algorithms For Linear Systems and Eigenvalue Problems. 
This implementation outperforms eigs() in time on Heisenberg Hamiltonian matrices of size $2^13 \times 2^13$ or greater. 



This work was supported in part by the U.S. department of Energy, Office of Science, Office of Workforce Development for Teachers and Students (WDTS) under the Science Undergraduate Laboratory Internship (SULI) Program. 
