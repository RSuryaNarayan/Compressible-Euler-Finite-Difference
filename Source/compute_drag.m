M = 0.4;
rho_inf = 1;
U_inf = M;
F_D =0;
for i=1:N-1
    M = [xi_x(i,1) xi_y(i,1);eta_x(i,1) eta_y(i,1)];
    if (i==1)
        D = [0; p_n(i,1)*sqrt((x(i+1,1)-x(N-1,1))^2+(y(i+1,1)-y(N-1,1))^2)];
    else
        D = [0; p_n(i,1)*sqrt((x(i+1,1)-x(i-1,1))^2+(y(i+1,1)-y(i-1,1))^2)];
    end
    A = M\D;
    F_x_i = A(1);
    F_D = F_D + F_x_i;
end

L=0;
for i=1:N-1
    L = L+sqrt((x(i+1,1)-x(i,1))^2+(y(i+1,1)-y(i,1))^2);
end

C_D = F_D/0.5/rho_inf/U_inf^2;