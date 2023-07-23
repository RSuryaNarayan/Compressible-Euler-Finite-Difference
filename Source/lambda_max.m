function l_max = lambda_max(u,v, p, rho, xi_x, xi_y, eta_x, eta_y)
    gamma = 1.4;
    a = sqrt(gamma*p/rho);
    if (~isreal(a))
        a
        rho
        p
        return;
    end
    l_xi = abs(xi_x*u+xi_y*v) + a*sqrt(xi_x^2 + xi_y^2);
    l_eta = abs(eta_x*u+eta_y*v) + a*sqrt(eta_x^2 + eta_y^2);
    l_max = [l_xi, l_eta];
end