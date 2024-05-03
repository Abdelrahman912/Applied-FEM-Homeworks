clc;
clear;

dt = 0.01; %sec
L = 1.0; %m
rhoA = 1; %kg/m
EA = 500; %N
g = 9.8; %m/s2

alpha_1 = 0;
alpha_2 = 0; 
p_inf = 0.1;

M = [(rhoA * L^3)/3 0 ; 0 (rhoA * L) /3];
K = [(rhoA * g * L^2)/2 0; 0 EA/L];
C = alpha_1 * M + alpha_2 * K;
R = [0 ; 0];

u0 = [-L/5 ; 0]; %initial displacement
v0 = [0 ; sqrt(g/(6*L))]; %initial velocity
a0 = M \ (R - C*v0 - K*u0); % initial acceleration

%Chung and Hulbert (1993)
alpha_m = (2 * p_inf - 1)/(p_inf + 1);
alpha_f = p_inf / (p_inf + 1);
beta = 0.25 * (1 - alpha_m + alpha_f)^2;
gamma = 0.5 - alpha_m + alpha_f;

K_eff = M * ((1-alpha_m)/(beta*dt^2)) + C * (gamma *(1-alpha_f)/(beta*dt)) + K *(1 - alpha_f); %slides pg. 63

t_f = 5.0;
t = 0:dt:t_f;



N = length(t);
u = zeros(2,N);
v = zeros(2,N);
a = zeros(2,N);

% arrays for the response
u(:,1) = u0;
v(:,1) = v0;
a(:,1) = a0;

% arrays for the error
e = zeros(1,N); %absolute error
eta = zeros(1,N); %relative error
com_err = zeros(1,N);


for i = 1:length(t)-1
    r_eff = -K * alpha_f * u(:,i) ... 
          + C * ((gamma*(1-alpha_f)/(beta * dt))* u(:,i) + ((gamma - gamma * alpha_f - beta)/(beta))*v(:,i) + (((gamma-2*beta)*(1 - alpha_f))/(2*beta))*dt*a(:,i)) ... 
          + M * (((1-alpha_m)/(beta * dt^2)) * u(:,i) + ((1-alpha_m)/(beta*dt)) * v(:,i) + ((1-alpha_m-2*beta)/(2*beta)) * a(:,i)) ;
    u(:,i+1) = K_eff\r_eff;
    % update v and a -> slides pg. 60
    v(:,i+1) = (gamma/(beta * dt))*(u(:,i+1) - u(:,i)) - ((gamma - beta)/beta)*v(:,i) - ((gamma-2*beta)/(2*beta))*dt*a(:,i);
    a(:,i+1) = (1/(beta*dt^2))*(u(:,i+1) - u(:,i)) - (1/(beta*dt)) * v(:,i) - ((1-2*beta)/(2*beta))*a(:,i);
    
    %calculate our errors
    e(i+1) = norm(((6 * beta - 1)/6) * (a(:,i+1) - a(:,i))*dt^2);
    eta(i+1) = (e(i+1))/norm(u(i+1) - u(i));
    com_err(i+1) = sum(e);
end

plot(t, u,'black')
hold on
plot(t,v,'blue')
hold on
%plot(t,a)




