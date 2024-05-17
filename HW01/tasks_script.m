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

% Froward Iteration
[lambda, phi, n ]= forward_iter(K,M,x1,TOL);
lambda % 9 -> (largest eigen value)
phi % [0.7071; -0.7071]
n % 5

% Inverse Iteration
[lambda, phi, n ]= inverse_iter(K,M,x1,TOL);
lambda % 1 -> (smallest eigen value)
phi % [0.7071; 0.7071]
n % 5

% Task 8:
x2 = [1;1]; % initial guess for Task 8.
[lambda, phi, n ]= forward_iter(K,M,x2,TOL);
lambda % 1 -> (smallest eigen value)
phi % [0.7071; 0.7071]
n % 2