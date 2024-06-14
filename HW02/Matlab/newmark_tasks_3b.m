clc;
clear;

fprintf('Homework 2: Generalized-Alpha method\n');
fprintf('Name: Abdelrahman Fathy Abdelhaleem Aly Abdelrahman \n');
fprintf('Matr.-Nr.: 108023251500 \n \n');
fprintf('Start Task 3b (and Task 4): \n');

% Data Inputs
L = 1.0 ;
pA = 1.0 ;
EA = 500;
g = 9.8;
dt = 0.01;
t_f = 5;

% Linear Momentum
M = [(pA*L^3)/3 0; 0 (pA*L)/3];
fprintf('Mass Matrix: \n');
disp(M);
K = [(pA*g*L^2)/2 0; 0 EA/L];
fprintf('Stiffness Matrix: \n');
disp(K);

% Initial Conditions
u0 = [0; -L/5];
fprintf('Initial Displacement: [%f, %f] \n',u0(1),u0(2));
v0 = [sqrt(g/6*L); 0];
fprintf('Initial Velocity: [%f, %f] \n',v0(1),v0(2));

% Task 3a:
alpha_1 = 0;
alpha_2 = 0;
p_inf = 0.1;

fprintf('Task 3b (and Task 4): \n');
[u, v , a  , e_abs, eta, e_cum, t, t_steps] = newmark(M,K,alpha_1,alpha_2,p_inf,u0,v0,t_f,dt); 



% Create a figure for each plot and save it
% Define a directory to save the plots
outputDir = 'task3b_plots';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Plot u
figure;
plot(t, u(1, :), 'r', t, u(2, :), 'b');
xlabel('Time');
ylabel('u');
title('u vs Time (3b)');
legend('\theta', 'u');
saveas(gcf, fullfile(outputDir, 'u_vs_time.png'));

% Plot v
figure;
plot(t, v(1, :), 'r', t, v(2, :), 'b');
xlabel('Time');
ylabel('v');
title('v vs Time (3b)');
legend('$\dot{\theta}$', '$\dot{u}$','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 'v_vs_time.png'));

% Plot a
figure;
plot(t, a(1, :), 'r', t, a(2, :), 'b');
xlabel('Time');
ylabel('a');
title('a vs Time (3b)');
legend('$\ddot{\theta}$', '$\ddot{u}$','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 'a_vs_time.png'));

% Plot e_abs
figure;
plot(t, e_abs, 'r');
xlabel('Time');
ylabel('e\_abs','Interpreter','latex');
title('e\_abs vs Time (3b)','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 'e_abs_vs_time.png'));

% Plot eta
figure;
plot(t, eta, 'r');
xlabel('Time');
ylabel('$\eta$','Interpreter','latex');
title('$\eta$ vs Time (3b)','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 'eta_vs_time.png'));

% Plot e_cum
figure;
plot(t, e_cum, 'r');
xlabel('Time');
ylabel('e\_cum','Interpreter','latex');
title('e\_cum vs Time (3b)','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 'e_cum_vs_time.png'));

% Plot t_steps
figure;
plot(t, t_steps, 'r');
xlabel('Time');
ylabel('t\_steps','Interpreter','latex');
title('t\_steps vs Time (3b)','Interpreter','latex');
saveas(gcf, fullfile(outputDir, 't_steps_vs_time.png'));

% Close all figures
close all;

fprintf('Task 3b (and Task 4): plots is saved in [%s] folder \n', outputDir);
