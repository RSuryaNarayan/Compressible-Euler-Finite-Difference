function dt = compute_dt(u,v,rho,p,xi_x, xi_y, eta_x, eta_y, dxi, deta)
 N = length(u);
 dt = 1e2;
 for i = 1:N
        for j=1:N
            lmax = lambda_max(u(i,j),v(i,j),rho(i,j),p(i,j),xi_x(i,j),xi_y(i,j),eta_x(i,j),eta_y(i,j));
            dt_min = min(dxi/lmax(1),deta/lmax(2));
            if(dt_min<dt)
                dt = dt_min;
            end
        end
 end
end