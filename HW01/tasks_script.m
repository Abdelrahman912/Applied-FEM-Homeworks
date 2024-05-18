clc;
clear;

m =25.0 ;
M_ = 700.0;
k = 3e3;
L = 1.0;
g = 9.8;

K = [((2*m + M_)*g)/L + k -k; -k ((2*m + M_)*g)/L + k]; % Stiffness Matrix
fprintf('Stiffness Matrix: \n');
disp(K);
M = [2*m+M_ 0;0 2*m+M_]; % Mass Matrix
fprintf('Mass Matrix: \n');
disp(M);
TOL = 1e-6; % Convergence tolerence.


% Task 7:

x1 = [1;0]; % initial guess for Task 7.

% Froward Iteration
[lambda, phi, n ]= forward_iter(K,M,x1,TOL);
fprintf('Forward Iteration: \n');
fprintf('eigen value: %f \n',lambda) % 17.8 -> (largest eigen value)
fprintf('eigen vector: [%f, %f] \n',phi(1),phi(2)) % [0.7071; -0.7071]
fprintf('number of iterations: %d \n \n',n) % 12

% Inverse Iteration
[lambda, phi, n ]= inverse_iter(K,M,x1,TOL);
fprintf('Inverse Iteration: \n');
fprintf('eigen value: %f \n',lambda) % 9.8 -> (smallest eigen value)
fprintf('eigen vector: [%f, %f] \n',phi(1),phi(2)) % [0.7071; 0.7071]
fprintf('number of iterations: %d \n \n',n) % 13


% Task 8:
x2 = [1;1]; % initial guess for Task 8.
[lambda, phi, n ]= forward_iter(K,M,x2,TOL);
fprintf('Forward Iteration with initial guess [%f, %f]: \n',x2(1),x2(2));
fprintf('eigen value: %f \n',lambda) % 9.8 -> (smallest eigen value)
fprintf('eigen vector: [%f, %f] \n',phi(1),phi(2)) % [0.7071; 0.7071]
fprintf('number of iterations: %d \n',n) % 2
