clc;
clear;

m =25.0 ;
M_ = 700.0;
k = 3e3;
L = 1.0;
g = 9.8;

K = [(2*m + M_)/L + k -k; -k (2*m + M_)/L + k]; % Stiffness Matrix
M = [2*m+M_ 0;0 2*m+M_]; % Mass Matrix
TOL = 1e-6; % Convergence tolerence.


% Task 7:

x1 = [1;0]; % initial guess for Task 7.

[lambda, phi, n ]= forward_iter(K,M,x1,TOL);
lambda
phi
n