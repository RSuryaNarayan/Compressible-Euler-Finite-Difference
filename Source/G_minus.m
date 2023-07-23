function gm = G_minus(u, v, p, rho, eta_x, eta_y,J)
    %returns the 3 column F^+ vector for the advance step 
    %all inputs are scalars hence values at a point and not vectors
    r_eta = sqrt(eta_x^2+eta_y^2);
    gamma = 1.4;
    n_x = eta_x/r_eta;
    n_y = eta_y/r_eta;
    u_d = u*n_x + v*n_y;
    q = (u^2+v^2)/2;
    a = sqrt(gamma*p/rho);
    h0 = q + a^2/(gamma-1);
    e = 1e-5;
    %construct lambda_plus matrix
    L_m = zeros(4,4);
    L_m(1,1) = r_eta*0.5*(u_d-a-sqrt((u_d-a)^2+e^2));
    L_m(2,2) = r_eta*0.5*(u_d-sqrt(u_d^2+e^2));
    L_m(3,3) = r_eta*0.5*(u_d+a-sqrt((u_d+a)^2+e^2));
    L_m(4,4) = r_eta*0.5*(u_d-sqrt(u_d^2+e^2));

    %construct left eigen vector
    R = zeros(4,4);
    R(1,1) = 1;
    R(1,2) = 1;
    R(1,3) = 1;
    
    R(2,1) = u-a*n_x;
    R(2,2) = u;
    R(2,3) = u+a*n_x;
    R(2,4) = n_y;

    R(3,1) = v-a*n_y;
    R(3,2) = v;
    R(3,3) = v+a*n_y;
    R(3,4) = -n_x;

    R(4,1) = h0 - a*u_d;
    R(4,2) = q;
    R(4,3) = h0 + a*u_d;
    R(4,4) = u*n_y - v*n_x;

    %construct left-egien vector
    R_inv = zeros(4,4);

    R_inv(1,1) = ((gamma-1)*q+a*u_d)/(2*a*a);
    R_inv(1,2) = ((1-gamma)*u-a*n_x)/(2*a*a);
    R_inv(1,3) = ((1-gamma)*v-a*n_y)/(2*a*a);
    R_inv(1,4) = ((gamma-1))/(2*a*a);

    R_inv(2,1) = (a^2-(gamma-1)*q)/(a*a);
    R_inv(2,2) = ((gamma-1)*u)/(a*a);
    R_inv(2,3) = ((gamma-1)*v)/(a*a);
    R_inv(2,4) = ((1-gamma))/(a*a);

    R_inv(3,1) = ((gamma-1)*q-a*u_d)/(2*a*a);
    R_inv(3,2) = ((1-gamma)*u+a*n_x)/(2*a*a);
    R_inv(3,3) = ((1-gamma)*v+a*n_y)/(2*a*a);
    R_inv(3,4) = ((gamma-1))/(2*a*a);

    R_inv(4,1) = v*n_x-u*n_y;
    R_inv(4,2) = n_y;
    R_inv(4,3) = -n_x;
    R_inv(4,4) = 0;

    %construct state vector:
    U = [rho;rho*u;rho*v;rho*q+p/(gamma-1)];
    U_Tilde = U/J;

    gm = R*L_m*R_inv*U_Tilde;
end