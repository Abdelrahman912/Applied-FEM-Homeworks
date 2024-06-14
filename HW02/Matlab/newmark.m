function [u, v , a  , e_abs, eta, e_cum, t, t_steps] = newmark(M,K,alpha_1,alpha_2,p_inf,u0,v0,t_f,dt)

% Rayleigh damping
C = alpha_1 * M + alpha_2 * K;

% initial acceleration
a0 = M \ (- C*v0 - K*u0); 

% Chung and Hulbert (1993)
alpha_m = (2 * p_inf - 1)/(p_inf + 1);
alpha_f = p_inf / (p_inf + 1);
beta = 0.25 * (1 - alpha_m + alpha_f)^2;
gamma = 0.5 - alpha_m + alpha_f;


% time stepping
% t_f : final time
% dt : time step
% t : time vector   
% t_steps : vector of time steps (for plotting)
t = 0:dt:t_f;
N = length(t);  
t_steps = dt * ones(1,N);

% arrays for the response
u = zeros(2,N);
v = zeros(2,N);
a = zeros(2,N);

% initial conditions
u(:,1) = u0;
v(:,1) = v0;
a(:,1) = a0;

% define the container for the errors
e_abs = zeros(1,N); %absolute error
eta = zeros(1,N); %relative error
e_cum = zeros(1,N); %cumulative error


% effective stiffness matrix (slides pg. 63)
K_eff = M * ((1-alpha_m)/(beta*dt^2)) + C * (gamma *(1-alpha_f)/(beta*dt)) + K *(1 - alpha_f); 

for i = 1:length(t)-1
    r_eff = -K * alpha_f * u(:,i) ... 
          + C * ((gamma*(1-alpha_f)/(beta * dt))* u(:,i) + ((gamma - gamma * alpha_f - beta)/(beta))*v(:,i) + (((gamma-2*beta)*(1 - alpha_f))/(2*beta))*dt*a(:,i)) ... 
          + M * (((1-alpha_m)/(beta * dt^2)) * u(:,i) + ((1-alpha_m)/(beta*dt)) * v(:,i) + ((1-alpha_m-2*beta)/(2*beta)) * a(:,i)) ;
   
    % solve for u at the next time step
    u(:,i+1) = K_eff\r_eff;

    % update v and a -> slides pg. 60
    v(:,i+1) = (gamma/(beta * dt))*(u(:,i+1) - u(:,i)) - ((gamma - beta)/beta)*v(:,i) - ((gamma-2*beta)/(2*beta))*dt*a(:,i);
    a(:,i+1) = (1/(beta*dt^2))*(u(:,i+1) - u(:,i)) - (1/(beta*dt)) * v(:,i) - ((1-2*beta)/(2*beta))*a(:,i);
    
    %calculate our errors
    e_abs(i+1) = norm(((6 * beta - 1)/6) * (a(:,i+1) - a(:,i))*dt^2);
    eta(i+1) = (e_abs(i+1))/norm(u(:,i+1) - u(:,i));
    e_cum(i+1) = sum(e_abs);
end


end







