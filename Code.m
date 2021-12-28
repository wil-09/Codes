function Code
% PSEUDOSPECTRAL METHOD FOR THE ADVECTION-CHEMOTAXIS MODEL IN 2D
clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
i = sqrt(-1);
% we allocate values to model parameters 
D1 = 4*10^(-6); D2 = 1.69*10^(-9); D3 = 8.9*10^(-6);
chi0 = 6.49.*10^(-5); tau0 = 1.6*10^(-3); 
beta =  3.5*10^(-8); sigma = 2.4231.*10.^(-3);
delta0x = 1.5*10^(-5); delta0y = 0.45*10^(-5);
% D1 = D3.*tau0./(tau0 - chi0) - 10.^(-7);
r = beta.*tau0./(tau0 - chi0) + 1.62*10^(-9);
beta1 = r.*sigma./10 + beta.*(D3./(tau0 - chi0) - D1./tau0)./10;

t0 = 0; T = 200;                                            % e-6;%1; % initial and final times simulation
dt = 0.01;                                                   % 1e-11;  % time step
nmax = round(T/dt)+1;                                        % Number of time iterations
t = 0:dt:T;                                                  % time vector
Nx = 128*2; Ny = 128*2; 
Lx = 50; Ly = 50; nplot = (T - t0)/(4*dt);                   % period 2*pi*L
% Gridspace and grid vector used in the Fourier method
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)*Lx;                            % x coordinate; S11 = length(x); 
X1 = x(2) - x(1)                                            % measuring the space step
kxm = [0:Nx/2-1 0 -Nx/2+1:-1]/Lx;                            % wave vector    1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)*Ly;                            % y coordinate
ky = [0:Ny/2-1 0 -Ny/2+1:-1]/Ly;                             % wave vector   1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;
[xx,yy] = meshgrid(x,y);
[k2xm,k2ym] = meshgrid(kxm.^2,ky.^2);
[kxm,kym] = meshgrid(kxm,ky);

% Parameters of the Analytical solutions used as initial conditions
P2 = 4;
sigma1 = beta1 - r.*sigma + beta.*(D3.*tau0 + D1.*(chi0 - tau0))./(tau0.*(chi0 - tau0));
sigma2 = r + beta*tau0./(chi0 - tau0); 

% traveling solitary pulses as initial conditions
K = sqrt(-sigma1./D2)./(4.*P2); l = sqrt(K./10);
beta3 = -2.*sigma2.^2./(15.*sigma1);
omega = delta0x.*sqrt(K - l.^2) + delta0y.*l;
V = @(xx,yy,t) -15.*sigma1.*(sech(sqrt(P2).*(sqrt(K - l.^2).*xx + l.*yy - omega.*t))).^2./(2.*sigma2); n0 = V(xx,yy,0.0); Tabn = [n0]; n_in = n0;
U = @(xx,yy,t) -D1./tau0 - D3./(chi0 - tau0) + 15.*sigma1.*tau0.*(sech(sqrt(P2).*(sqrt(K - l.^2).*xx + l.*yy - omega.*t))).^2./(2.*sigma2.*(chi0 - tau0)); c0 = U(xx,yy,0.0); Tabc = [c0]; c_in = c0;  Tabx = []; Taby = [];

% these are the noise strength we use to perturbe the initial conditions for both variables respectively
epsn1 = max(n0(:)); epsc1 = max(c0(:)); 
epsn2 = min(n0(:)); epsc2 = min(c0(:)); 

% Amplitudes of analytical solutions at time t = 0
epsn = epsn1 - epsn2;
epsc = epsc1 - epsc2;

% We perturb initial conditions with uniformly distributed random noises of strength a proportional to their amplitudes at time t = 0.
v = fftn(n0) + 0.1*epsn*randn(Nx,Ny); 
c = fftn(c0) + 0.1*epsc*randn(Nx,Ny);

% graphical output Configurations for all the code
figure
subplot(1,2,1)
mesh(x,y,n_in); view(-144,58)
xlabel 'x', ylabel 'y', zlabel 'n(x,y,0)'

subplot(1,2,2)
mesh(x,y,c_in); view(-144,58) 
xlabel 'x', ylabel 'y', zlabel 'c(x,y,0)'

% parameters of the system in the Fourier space
c_1n = r.*sigma - D1.*(k2xm + k2ym) - D2.*(k2xm + k2ym).^2 - 1i.*kxm.*delta0x - 1i.*delta0y.*kym;
c_2 = -i.*kxm.*tau0; c_3 = -i.*kym.*tau0; c_4 = -i.*kxm.*chi0; c_5 = -i.*kym.*chi0;
c_1c = beta1 - D3.*(k2xm + k2ym) - 1i.*kxm.*delta0x - 1i.*delta0y.*kym; 
c_2c = -beta3; c_4c = -beta;

% Right hand side of the first equation
RHS1n = @(v) c_1n.*v;
RHS2n = @(v) c_2.*fftn(ifftn(v).*ifftn(1i.*kxm.*v));
RHS3n = @(v) c_3.*fftn(ifftn(v).*ifftn(1i.*kym.*v));
RHS4n = @(v,c) c_4.*fftn(ifftn(v).*ifftn(1i.*kxm.*c));
RHS5n = @(v,c) c_5.*fftn(ifftn(v).*ifftn(1i.*kym.*c));
RHS6n = @(v) -r.*fftn((ifftn(v)).^2);
RHSn = @(v,c) RHS1n(v) + RHS2n(v) + RHS3n(v) + RHS4n(v,c) + RHS5n(v,c) + RHS6n(v);

% the Right hand side of the second equation
RHS1c = @(c) c_1c.*c;
RHS2c = @(v,c) c_2.*fftn(ifftn(c).*ifftn(i.*kxm.*v));
RHS3c = @(v,c) c_3.*fftn(ifftn(c).*ifftn(i.*kym.*v));
RHS4c = @(v) c_2c.*fftn(ifftn(v).^3);
RHS6c = @(v,c) c_4c.*fftn(ifftn(c).*ifftn(v));
RHSc = @(v,c) RHS1c(c) + RHS2c(v,c) + RHS3c(v,c) + RHS4c(v) + RHS6c(v,c);

% loop in time and integrate the differential equation in the Fourier space
for n = 1:nmax
    k1n = dt*RHSn(v,c);
    k1c = dt*RHSc(v,c);
    k2n = dt*RHSn(v + k1n/2, c + k1c/2);
    k2c = dt*RHSc(v + k1n/2, c + k1c/2);
    k3n = dt*RHSn(v + k2n/2, c + k2c/2);
    k3c = dt*RHSc(v + k2n/2, c + k2c/2);
    k4n = dt*RHSn(v + k3n, c + k3c);
    k4c = dt*RHSc(v + k3n, c + k3c);
    v = v + (k1n + 2*k2n + 2*k3n + k4n)/6;
    c = c + (k1c + 2*k2c + 2*k3c + k4c)/6;
    
%     return to the real space
    n_fin = real(ifftn(v));
    c_fin = real(ifftn(c));
    if(mod(n,nplot)==0)
        figure
        subplot(1,2,1)
        mesh(x,y,n_fin); view(-144,58)
        xlabel 'x'; ylabel 'y'; zlabel 'n(x,y,t)'
        
        subplot(1,2,2)
        mesh(x,y,c_fin); view(-144,58)
        xlabel 'x'; ylabel 'y'; zlabel 'c(x,y,t)'       
    end
    t = t + dt;
end
