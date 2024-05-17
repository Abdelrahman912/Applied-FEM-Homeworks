function [lambda,phi,n] = inverse_iter(K,M,x,TOL)


    % K: stifnness matrix
    % M: mass matrix
    % x: initial guess
    % TOL: convergence tolerence

    y = M*x;

    % any number above the tolerence in order not to converge
    % from the first step.
    err = TOL*2 ;
    rho_new = 0;
    rho_old = 0;
    n = 0; % counter to count number of loops before convergence.

    while err >= TOL
        n = n+1;
        xbar = K\y;
        ybar = M*xbar;
        rho_old = rho_new;
        rho_new = (xbar'*y)/(xbar'*ybar);
        err = abs(rho_new - rho_old) / rho_new;
        y = ybar/sqrt(xbar' * ybar);
    end

    lambda = rho_new;
    phi = M\y;
    phi = phi/norm(phi);

end