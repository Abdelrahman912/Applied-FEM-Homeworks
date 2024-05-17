function [lambda, phi, n ] =  forward_iter(K,M,x,TOL)

    % K: stifnness matrix
    % M: mass matrix
    % x: initial guess
    % TOL: convergence tolerence

    y = K*x;
    % any number above the tolerence in order not to converge
    % from the first step.
    err = TOL*2 ;
    rho_new = 0;
    rho_old = 0;
    n = 0; % counter to count number of loops before convergence.
    
    while err >= TOL
        n = n+1;
        xbar = M\y;
        ybar = K * xbar;
        rho_old = rho_new;
        rho_new = (xbar'*ybar)/(xbar'*y);
        err = abs(rho_new - rho_old)/(rho_new);
        y = ybar/sqrt(xbar'*y);
    end
    lambda = rho_new; % eigenvalue
    phi = K\y; 
    phi = phi/norm(phi); %normalized eigenvector

end