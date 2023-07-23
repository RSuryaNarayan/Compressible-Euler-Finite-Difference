%---------------------Author and course info----------------------------%
% Name: Suryanarayan Ramachandran
% Course: AEM5253 Computational Fluid Mechanics
% Date: Dec 26th, 2022
% E-Mail: ramac106@umn.edu
% Course instructor: Dr.Pramod K Subbareddy

%---------------------Preamble------------------------------------------%
clear all; close all; clc; %#ok<CLALL> 

%----------------------Mesh generation----------------------------------%
%# of points along xi/eta
N=100;

%generate xi, eta
xi = linspace(0,2*pi,N);
eta = linspace(atanh(0.25),6.2,N);

%meshgrid to get structured mesh
[eta, xi] = meshgrid(eta, xi);

%get physical space coords
x = cosh(eta).*cos(xi);
y = sinh(eta).*sin(xi);

%---------------------Mesh Metrics--------------------------------------%
%get Delta(xi), Delta(eta)
dxi = xi(2,1)-xi(1,1);
deta = eta(1,2)-eta(1,1);

%Declare arrays for mesh-metrics and Jacobians
xi_x = zeros(N,N);  
xi_y = zeros(N,N);
eta_x = zeros(N,N);
eta_y = zeros(N,N);
J = zeros(N,N);

%Compute mesh-metrics and Jacobians
for i=1:N
    for j=1:N
        if(i>=2 && i<=N-1)
            x_xi = (x(i+1,j)-x(i-1,j))/2/dxi;
            y_xi = (y(i+1,j)-y(i-1,j))/2/dxi;
        elseif(i==1)
            x_xi = (x(i+1,j)-x(i,j))/dxi;
            y_xi = (y(i+1,j)-y(i,j))/dxi;
        elseif(i==N)
            x_xi = (x(i,j)-x(i-1,j))/dxi;
            y_xi = (y(i,j)-y(i-1,j))/dxi;
   
        
        end

        if(j>=2 && j<=N-1)
            x_eta = (x(i,j+1)-x(i,j-1))/2/deta;
            y_eta = (y(i,j+1)-y(i,j-1))/2/deta;
        elseif(j==1)
            x_eta = (x(i,j+1)-x(i,j))/deta;
            y_eta = (y(i,j+1)-y(i,j))/deta;
        elseif(j==N)
            x_eta = (x(i,j)-x(i,j-1))/deta;
            y_eta = (y(i,j)-y(i,j-1))/deta;
        end

        J(i,j) = 1/(x_xi*y_eta - x_eta*y_xi);
        xi_x(i,j) = y_eta*J(i,j); 
        xi_y(i,j) = -y_xi*J(i,j);
        eta_x(i,j)= -x_eta*J(i,j);
        eta_y(i,j)= x_xi*J(i,j);
    end
end

%----------------------Declare matrices for the solver--------------------%
%Time-level: n arrays
rho = zeros(N,N);
u = zeros(N,N);
v = zeros(N,N);
p = zeros(N,N);
%Time-level:n+1 arrays
rho_n = zeros(N,N);
u_n = zeros(N,N);
v_n = zeros(N,N);
p_n = zeros(N,N);

%matrices for storing norms
u_norm = [];
v_norm = [];
rho_norm = [];
p_norm = [];
time = [];

%--------------------Initial conditions-------------------------------%
%freestream conditions
M=0.40;
gamma = 1.4;
p_inflow = 1/gamma;
rho_inflow=1;
u_inflow = M*sqrt(gamma*p_inflow/rho_inflow);
cfl=0.5; 

%Initialize
for i=1:N
    for j=1:N
        u(i,j)= u_inflow;
        p(i,j)=p_inflow;
        rho(i,j)=rho_inflow;
        v(i,j)=0;
    end
end

%------------------Main Solution Loop--------------------------------%
%set time:
t=0;
%set max time of simulation:
t_max=10;

%Solve:
%Index notation followed: i for xi direction, j for eta direction
while(t<t_max)
    %compute minimum dt:
    dt_min = cfl*compute_dt(u,v,p,rho, xi_x, xi_y, eta_x, eta_y, dxi, deta);

    %advance i=1 and i=N, in the interior j's 
    % For i=1, i+1=2, i-1 = N-1
    for j=2:N-1
            %Compute RHS
            R_ij = (F_plus(u(1,j),v(1,j),p(1,j),rho(1,j),xi_x(1,j),xi_y(1,j),J(1,j))...
                -F_plus(u(N-1,j),v(N-1,j),p(N-1,j),rho(N-1,j),xi_x(N-1,j),xi_y(N-1,j),J(N-1,j)))/dxi...
                +(F_minus(u(2,j),v(2,j),p(2,j),rho(2,j),xi_x(2,j),xi_y(2,j),J(2,j))...
                -F_minus(u(1,j),v(1,j),p(1,j),rho(1,j),xi_x(1,j),xi_y(1,j),J(1,j)))/dxi...
                +(F_plus(u(1,j),v(1,j),p(1,j),rho(1,j),eta_x(1,j),eta_y(1,j),J(1,j))...
                -F_plus(u(1,j-1),v(1,j-1),p(1,j-1),rho(1,j-1),eta_x(1,j-1),eta_y(1,j-1),J(1,j-1)))/deta...
                +(F_minus(u(1,j+1),v(1,j+1),p(1,j+1),rho(1,j+1),eta_x(1,j+1),eta_y(1,j+1),J(1,j+1))...
                -F_minus(u(1,j),v(1,j),p(1,j),rho(1,j),eta_x(1,j),eta_y(1,j),J(1,j)))/deta;

            %Locally construct conserved state vector
            U_old = [rho(1,j);rho(1,j)*u(1,j);rho(1,j)*v(1,j);rho(1,j)*(0.5*u(1,j)^2+0.5*v(1,j)^2)+p(1,j)/(gamma-1)];

            %advance
            U_new = U_old -dt_min*J(1,j)*R_ij;

            %get back new primitive vars
            rho_new = U_new(1);
            u_new = U_new(2)/rho_new;
            v_new = U_new(3)/rho_new;
            q_new = 0.5*(u_new^2+v_new^2);
            p_new = (gamma-1)*(U_new(4)-rho_new*q_new);

            %assign to new arrays
            % since xi-direction is periodic, we note i=1 and i=N are the same 
            rho_n(1,j) = rho_new;
            rho_n(N,j) = rho_new;
            u_n(1,j) = u_new;
            u_n(N,j) = u_new;
            v_n(1,j) = v_new;
            v_n(N,j) = v_new;
            p_n(1,j) = p_new;
            p_n(N,j) = p_new;
    end

    %advance all other interior j's and i's i.e. 2:N-1
    for j=2:N-1 %iterate through the interior ellipses-outer loop, we will fill j=1, N later
        for i=2:N-1 %iterate along the interior ellipses-inner loop, we will fill i=1,N later
                %compute RHS
                R_ij = (F_plus(u(i,j),v(i,j),p(i,j),rho(i,j),xi_x(i,j),xi_y(i,j),J(i,j))...
                       -F_plus(u(i-1,j),v(i-1,j),p(i-1,j),rho(i-1,j),xi_x(i-1,j),xi_y(i-1,j),J(i-1,j)))/dxi...
                      +(F_minus(u(i+1,j),v(i+1,j),p(i+1,j),rho(i+1,j),xi_x(i+1,j),xi_y(i+1,j),J(i+1,j))...
                       -F_minus(u(i,j),v(i,j),p(i,j),rho(i,j),xi_x(i,j),xi_y(i,j),J(i,j)))/dxi...
                      +(F_plus(u(i,j),v(i,j),p(i,j),rho(i,j),eta_x(i,j),eta_y(i,j),J(i,j))...
                       -F_plus(u(i,j-1),v(i,j-1),p(i,j-1),rho(i,j-1),eta_x(i,j-1),eta_y(i,j-1),J(i,j-1)))/deta...
                      +(F_minus(u(i,j+1),v(i,j+1),p(i,j+1),rho(i,j+1),eta_x(i,j+1),eta_y(i,j+1),J(i,j+1))...
                       -F_minus(u(i,j),v(i,j),p(i,j),rho(i,j),eta_x(i,j),eta_y(i,j),J(i,j)))/deta;
                
                %construct U_old
                U_old = [rho(i,j);rho(i,j)*u(i,j);rho(i,j)*v(i,j);0.5*rho(i,j)*(u(i,j)^2+v(i,j)^2)+p(i,j)/(gamma-1)];
                
                %advance
                U_new = U_old -dt_min*J(i,j)*R_ij;

                %extract fields from solution
                rho_new = U_new(1);
                u_new = U_new(2)/rho_new;
                v_new = U_new(3)/rho_new;
                q_new = 0.5*(u_new^2+v_new^2);
                p_new = (gamma-1)*(U_new(4)-rho_new*q_new);

                %assign to new arrays: 
                rho_n(i,j) = rho_new;
                u_n(i,j) = u_new;
                v_n(i,j) = v_new;
                p_n(i,j) = p_new;
        end
    end

    %advance j=1 and j=N
    for i=1:N
       %j=1 i.e. airfoil surface
       M = [xi_x(i,1) xi_y(i,1);eta_x(i,1) eta_y(i,1)];
       D = [xi_x(i,2)*u_n(i,2)+xi_y(i,2)*v_n(i,2); 0];
       A=M\D;
       %get velocities
       u_n(i,1)= A(1); v_n(i,1)=A(2);
       %get density
       rho_n(i,1) = rho_n(i,2);
       %get pressure
       q_1 = 0.5*(u_n(i,1)^2 + v_n(i,1)^2);
       q_2 = 0.5*(u_n(i,2)^2 + v_n(i,2)^2);
       p_n(i,1) = p_n(i,2) + (gamma-1)*(rho_n(i,2)*q_2-rho_n(i,1)*q_1);

       %j=N, i.e. outer boundary
        if(i>N/3 && i<=3*N/4) 
            %one half of the domain is an inlet
            u_n(i,N)=u_inflow;
            v_n(i,N)=0;
            p_n(i,N)=p_inflow;
            rho_n(i,N)=rho_inflow;
        else 
            %other half is an outlet
            u_n(i,N)=u_n(i,N-1);
            v_n(i,N)=v_n(i,N-1);
            p_n(i,N)=p_n(i,N-1);
            rho_n(i,N)=rho_n(i,N-1);
        end       
    end
    
    %impose periodicity along i once again
    u_n(N,:)=u_n(1,:);
    v_n(N,:)=v_n(1,:);
    p_n(N,:)=p_n(1,:);
    rho_n(N,:)=rho_n(1,:);

    %compute norms
    L2_u = norm(u-u_n,2);
    L2_v = norm(v-v_n,2);
    L2_p = norm(p-p_n,2);
    L2_rho = norm(rho-rho_n,2);

    %update norm vectors
    u_norm = [L2_u, u_norm];
    v_norm = [L2_v, v_norm];
    p_norm = [L2_p, p_norm];
    rho_norm = [L2_rho, rho_norm];

    %update solution arrays
    u = u_n;
    v = v_n;
    p = p_n;
    rho = rho_n;

    %advance time
    t = t+dt_min;
    %display time
    t
    %store time
    time = [t,time];
end

%plot residuals of primitive variables vs time
loglog(time,p_norm,time,rho_norm,time,u_norm,time,v_norm,'LineWidth',2);
xlabel("time","interpreter","latex");
ylabel("$L_2$-norm","interpreter","latex");
legend("$p$","$\rho$","u","v","interpreter","latex");
title("Residuals vs. time","interpreter","latex");
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/norm","epsc");