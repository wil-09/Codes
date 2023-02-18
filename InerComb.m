function InerComb
% this code numerically appromaximates the solutions of the hyperbolic-parabolic chemotaxis model
% I wrote the code on the 14th February 2023 at 12:55AM

tic
clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
format long g;

% Gridspace
xlim = 15; x0 = -xlim; xf = xlim; dx = 0.01; t0 = 0; tf = 1;
dt = dx*(tf - t0)/(xf - x0); x = x0:dx:xf; t = t0:dt:tf; 
N = numel(x); M = numel(t);

% Parameters of the model
alpha0 = 0.25; beta = 100; beta0 = 0.5; gamma = 100;
chi0 = 80; u = 10; r = 0.1; D2 = 1; tau = 10^3;

% Parameters characterizing the medium porosity
alpha1 = 0.05; alpha2 = 0.02; m = 5; m0 = 10; 
% alpha1 = alpha1.^m; alpha2 = alpha2.^m0;

% Approximated perturbed solutions
K = 0.05; eps = 0.01;

for m0 = 0
% The steady state of the system reads
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);

omega1 = (1./6).*( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (6.*((1./3).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau) - (1./9).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2))./( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (1./3).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u);
Vg1 = (-(2.*K.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega1.^2 - (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + i.*u.*K).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u) + (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u)./tau).*omega1 - (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + i.*u.*K)./tau + 2.*chi0.*n0.*beta0.*K./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega1.^2 + (r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*u.*K) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*u.*K)./tau + (2.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*K)).*omega1);
omegar = real(omega1); omegai = imag(omega1); Vg = real(Vg1); % Vgi = imag(Vg); 

% % % % % n_xxxx = (n(j - 2) +6*n(j) -4*n(j + 1) + n(j + 2) - 4*n(j - 1))/h^4;
nn = zeros(1,N); JJ = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for k = 1:M
    for j = 1:N
        nn(k,j) = n0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        cc(k,j) = c0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        Er0(k,j) = 0*nn(k,j);
    end
end

for k = 1:M-1
    for j = 1:N
        JJ(k,j) = (nn(k + 1,j) - nn(k,j))./dt; 
    end
end

% initialisation of the intermediate variables
J(1:N) = JJ(1,1:N); n(1:N) = nn(1,1:N); n00(1:N) = nn(1,1:N); 
c(1:N) = cc(1,1:N); Er(1:N) = Er0(1,1:N);
k1n = zeros(1,N); k2n = zeros(1,N); k3n = zeros(1,N); k4n = zeros(1,N);
k1c = zeros(1,N); k2c = zeros(1,N); k3c = zeros(1,N); k4c = zeros(1,N);
k1j = zeros(1,N); k2j = zeros(1,N); k3j = zeros(1,N); k4j = zeros(1,N);
n1 = zeros(1,N); n2 = zeros(1,N); n3 = zeros(1,N);
c1 = zeros(1,N); c2 = zeros(1,N); c3 = zeros(1,N);
J1 = zeros(1,N); J2 = zeros(1,N); J3 = zeros(1,N);

Nf = []; Cf = []; Jf = []; Err = [];
for k = 1:M
    for j = 1:N
        if(j == 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == 2)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N - 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        else
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*(n(j - 2) + n(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./beta0 - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j - 1)).*(c(j + 1) - c(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 1) - 2.*c(j) + c(j - 1))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c(j - 1) - dt.*(1 + 2./dx.^2).*c(j) - dt.*alpha2.*(c(j + 1).^(1 + m0) - 2.*c(j).^(1 + m0) + c(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*(n1(j - 2) + n1(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n1(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./beta0 - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j - 1)).*(c1(j + 1) - c1(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 1) - 2.*c1(j) + c1(j - 1))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c1(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c1(j - 1) - dt.*(1 + 2./dx.^2).*c1(j) - dt.*alpha2.*(c1(j + 1).^(1 + m0) - 2.*c1(j).^(1 + m0) + c1(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*(n2(j - 2) + n2(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n2(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./beta0 - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j - 1)).*(c2(j + 1) - c2(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 1) - 2.*c2(j) + c2(j - 1))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c2(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c2(j - 1) - dt.*(1 + 2./dx.^2).*c2(j) - dt.*alpha2.*(c2(j + 1).^(1 + m0) - 2.*c2(j).^(1 + m0) + c2(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*(n3(j - 2) + n3(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n3(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./beta0 - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j - 1)).*(c3(j + 1) - c3(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 1) - 2.*c3(j) + c3(j - 1))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c3(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c3(j - 1) - dt.*(1 + 2./dx.^2).*c3(j) - dt.*alpha2.*(c3(j + 1).^(1 + m0) - 2.*c3(j).^(1 + m0) + c3(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));
            
        end
        
        n(j) = n(j) + dt.*(k1n(j) + 2.*k2n(j) + 2.*k3n(j) + k4n(j))./6;
        J(j) = J(j) + dt.*(k1j(j) + 2.*k2j(j) + 2.*k3j(j) + k4j(j))./6;
        c(j) = c(j) + dt.*(k1c(j) + 2.*k2c(j) + 2.*k3c(j) + k4c(j))./6;
        
        Nf(k,j) = n(j); Jf(k,j) = J(j); Cf(k,j) = c(j);
        Err(k,j) = abs(Nf(k,j) - nn(k,j));
%         Err(k,j) = Er(j);
        Err1(k,j) = Err(k,j);
    end
end
plot(x, Nf(round(3.*M./4),:), 'linewidth', 3); hold on
end

for m0 = 300
% The steady state of the system reads
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);

omega1 = (1./6).*( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (6.*((1./3).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau) - (1./9).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2))./( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (1./3).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u);
Vg1 = (-(2.*K.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega1.^2 - (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + i.*u.*K).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u) + (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u)./tau).*omega1 - (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + i.*u.*K)./tau + 2.*chi0.*n0.*beta0.*K./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega1.^2 + (r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*u.*K) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*u.*K)./tau + (2.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*K)).*omega1);
omegar = real(omega1); omegai = imag(omega1); Vg = real(Vg1); % Vgi = imag(Vg); 

% % % % % n_xxxx = (n(j - 2) +6*n(j) -4*n(j + 1) + n(j + 2) - 4*n(j - 1))/h^4;
nn = zeros(1,N); JJ = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for k = 1:M
    for j = 1:N
        nn(k,j) = n0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        cc(k,j) = c0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        Er0(k,j) = 0*nn(k,j);
    end
end

for k = 1:M-1
    for j = 1:N
        JJ(k,j) = (nn(k + 1,j) - nn(k,j))./dt; 
    end
end

% initialisation of the intermediate variables
J(1:N) = JJ(1,1:N); n(1:N) = nn(1,1:N); n00(1:N) = nn(1,1:N); 
c(1:N) = cc(1,1:N); Er(1:N) = Er0(1,1:N);
k1n = zeros(1,N); k2n = zeros(1,N); k3n = zeros(1,N); k4n = zeros(1,N);
k1c = zeros(1,N); k2c = zeros(1,N); k3c = zeros(1,N); k4c = zeros(1,N);
k1j = zeros(1,N); k2j = zeros(1,N); k3j = zeros(1,N); k4j = zeros(1,N);
n1 = zeros(1,N); n2 = zeros(1,N); n3 = zeros(1,N);
c1 = zeros(1,N); c2 = zeros(1,N); c3 = zeros(1,N);
J1 = zeros(1,N); J2 = zeros(1,N); J3 = zeros(1,N);

Nf = []; Cf = []; Jf = []; Err = [];
for k = 1:M
    for j = 1:N
        if(j == 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == 2)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N - 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        else
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*(n(j - 2) + n(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./beta0 - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j - 1)).*(c(j + 1) - c(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 1) - 2.*c(j) + c(j - 1))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c(j - 1) - dt.*(1 + 2./dx.^2).*c(j) - dt.*alpha2.*(c(j + 1).^(1 + m0) - 2.*c(j).^(1 + m0) + c(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*(n1(j - 2) + n1(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n1(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./beta0 - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j - 1)).*(c1(j + 1) - c1(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 1) - 2.*c1(j) + c1(j - 1))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c1(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c1(j - 1) - dt.*(1 + 2./dx.^2).*c1(j) - dt.*alpha2.*(c1(j + 1).^(1 + m0) - 2.*c1(j).^(1 + m0) + c1(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*(n2(j - 2) + n2(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n2(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./beta0 - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j - 1)).*(c2(j + 1) - c2(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 1) - 2.*c2(j) + c2(j - 1))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c2(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c2(j - 1) - dt.*(1 + 2./dx.^2).*c2(j) - dt.*alpha2.*(c2(j + 1).^(1 + m0) - 2.*c2(j).^(1 + m0) + c2(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*(n3(j - 2) + n3(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n3(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./beta0 - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j - 1)).*(c3(j + 1) - c3(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 1) - 2.*c3(j) + c3(j - 1))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c3(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c3(j - 1) - dt.*(1 + 2./dx.^2).*c3(j) - dt.*alpha2.*(c3(j + 1).^(1 + m0) - 2.*c3(j).^(1 + m0) + c3(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));
            
        end
        
        n(j) = n(j) + dt.*(k1n(j) + 2.*k2n(j) + 2.*k3n(j) + k4n(j))./6;
        J(j) = J(j) + dt.*(k1j(j) + 2.*k2j(j) + 2.*k3j(j) + k4j(j))./6;
        c(j) = c(j) + dt.*(k1c(j) + 2.*k2c(j) + 2.*k3c(j) + k4c(j))./6;
        
        Nf(k,j) = n(j); Jf(k,j) = J(j); Cf(k,j) = c(j);
        Err(k,j) = abs(Nf(k,j) - nn(k,j));
%         Err(k,j) = Er(j);
        Err2(k,j) = Err(k,j);
    end
end

plot(x, Nf(round(3.*M./4),:), 'linewidth', 3); hold on
end

for m0 = 10^3
% The steady state of the system reads
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);

omega1 = (1./6).*( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (6.*((1./3).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau) - (1./9).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2))./( - 8.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3 + 36.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u) - 108.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*sqrt(12.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^3.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) - 3.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).^2.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^2 - 54.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u).*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)) + 12.*((r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*K.*u) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u)./tau).^3 + 81.*((r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*K.*u).*(1 + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u)./tau - chi0.*n0.*beta0.*K.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2)).^2)).^(1./3) - (1./3).*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*K.*u);
Vg1 = (-(2.*K.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega1.^2 - (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + i.*u.*K).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*K - 2.*alpha2.*K.*c0.^m0 + i.*u) + (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u)./tau).*omega1 - (2.*K.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 2.*K.^3.*D2 + 1i.*u).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + i.*u.*K)./tau + 2.*chi0.*n0.*beta0.*K./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega1.^2 + (r + 1./tau).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0 + 1i.*u.*K) + (r + K.^2.*(alpha0 - alpha1.*n0.^m + D2.*K.^2) + 1i.*u.*K)./tau + (2.*(1 + r + 1./tau + K.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*K)).*omega1);
omegar = real(omega1); omegai = imag(omega1); Vg = real(Vg1); % Vgi = imag(Vg); 

% % % % % n_xxxx = (n(j - 2) +6*n(j) -4*n(j + 1) + n(j + 2) - 4*n(j - 1))/h^4;
nn = zeros(1,N); JJ = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for k = 1:M
    for j = 1:N
        nn(k,j) = n0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        cc(k,j) = c0 + eps.*(sech(x(j) - Vg.*t(k))).*cos(K.*x(j) - omegai.*t(k)).*exp(omegar.*t(k));
        Er0(k,j) = 0*nn(k,j);
    end
end

for k = 1:M-1
    for j = 1:N
        JJ(k,j) = (nn(k + 1,j) - nn(k,j))./dt; 
    end
end

% initialisation of the intermediate variables
J(1:N) = JJ(1,1:N); n(1:N) = nn(1,1:N); n00(1:N) = nn(1,1:N); 
c(1:N) = cc(1,1:N); Er(1:N) = Er0(1,1:N);
k1n = zeros(1,N); k2n = zeros(1,N); k3n = zeros(1,N); k4n = zeros(1,N);
k1c = zeros(1,N); k2c = zeros(1,N); k3c = zeros(1,N); k4c = zeros(1,N);
k1j = zeros(1,N); k2j = zeros(1,N); k3j = zeros(1,N); k4j = zeros(1,N);
n1 = zeros(1,N); n2 = zeros(1,N); n3 = zeros(1,N);
c1 = zeros(1,N); c2 = zeros(1,N); c3 = zeros(1,N);
J1 = zeros(1,N); J2 = zeros(1,N); J3 = zeros(1,N);

Nf = []; Cf = []; Jf = []; Err = [];
for k = 1:M
    for j = 1:N
        if(j == 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == 2)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j)).*(c(j + 1) - c(j)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 2) - 2.*c(j + 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j + 2)).^(1 + m) - 2.*(n(j + 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c(j) - dt.*alpha2.*(c(j + 2).^(1 + m0) - 2.*c(j + 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j)).*(c1(j + 1) - c1(j)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 2) - 2.*c1(j + 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j + 2)).^(1 + m) - 2.*(n1(j + 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c1(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c1(j) - dt.*alpha2.*(c1(j + 2).^(1 + m0) - 2.*c1(j + 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j)).*(c2(j + 1) - c2(j)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 2) - 2.*c2(j + 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j + 2)).^(1 + m) - 2.*(n2(j + 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c2(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c2(j) - dt.*alpha2.*(c2(j + 2).^(1 + m0) - 2.*c2(j + 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j + 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j + 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j + 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) - u./(tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau + u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j)).*(c3(j + 1) - c3(j)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 2) - 2.*c3(j + 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j + 2)).^(1 + m) - 2.*(n3(j + 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j + 2)./dx.^2 - dt.*(2./dx.^2 + u./dx).*c3(j + 1) + dt.*(1./dx.^2 - 1 + u./dx).*c3(j) - dt.*alpha2.*(c3(j + 2).^(1 + m0) - 2.*c3(j + 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N - 1)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        elseif(j == N)
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*n(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./(beta0) - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j) - n(j - 1)).*(c(j) - c(j - 1)))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j - 2) - 2.*c(j - 1) + c(j))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + 2.*dt.*chi0.*beta.*n(j).*(c(j) - c(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*c(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c(j) - dt.*alpha2.*(c(j - 2).^(1 + m0) - 2.*c(j - 1).^(1 + m0) + c(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*n1(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n1(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n1(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./(beta0) - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j) - n1(j - 1)).*(c1(j) - c1(j - 1)))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j - 2) - 2.*c1(j - 1) + c1(j))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + 2.*dt.*chi0.*beta.*n1(j).*(c1(j) - c1(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*c1(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c1(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c1(j) - dt.*alpha2.*(c1(j - 2).^(1 + m0) - 2.*c1(j - 1).^(1 + m0) + c1(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*n2(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n2(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n2(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./(beta0) - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j) - n2(j - 1)).*(c2(j) - c2(j - 1)))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j - 2) - 2.*c2(j - 1) + c2(j))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + 2.*dt.*chi0.*beta.*n2(j).*(c2(j) - c2(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*c2(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c2(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c2(j) - dt.*alpha2.*(c2(j - 2).^(1 + m0) - 2.*c2(j - 1).^(1 + m0) + c2(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*n3(j - 4)./(tau.*dx.^4) + 4.*dt.*D2.*n3(j - 3)./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j - 2) + dt.*(4.*D2./(tau.*dx.^4) - 2.*alpha0./(tau.*dx.^2) + u./(tau.*dx)).*n3(j - 1) + dt.*(alpha0./(tau.*dx.^2) - D2./(tau.*dx.^4) + r./tau - u./(tau.*dx)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./(beta0) - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j) - n3(j - 1)).*(c3(j) - c3(j - 1)))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j - 2) - 2.*c3(j - 1) + c3(j))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + 2.*dt.*chi0.*beta.*n3(j).*(c3(j) - c3(j - 1)).^2./(tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*c3(j - 2)./dx.^2 + dt.*(u./dx - 2./dx.^2).*c3(j - 1) + dt.*(1./dx.^2 - 1 - u./dx).*c3(j) - dt.*alpha2.*(c3(j - 2).^(1 + m0) - 2.*c3(j - 1).^(1 + m0) + c3(j).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));

        else
            k1n(j) = dt.*J(j);
            k1j(j) = -dt.*D2.*(n(j - 2) + n(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n(j) + dt.*(r - 1./tau).*J(j) - 2.*dt.*r.*n(j).*J(j)./beta0 - dt.*r.*n(j).^2./(tau.*beta0) - dt.*chi0.*((n(j + 1) - n(j - 1)).*(c(j + 1) - c(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c(j)).^2) - dt.*chi0.*n(j).*(c(j + 1) - 2.*c(j) + c(j - 1))./(tau.*dx.^2.*(1 + beta.*c(j)).^2) + dt.*chi0.*beta.*n(j).*(c(j + 1) - c(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c(j)).^3) - dt.*alpha1.*((n(j - 2)).^(1 + m) - 2.*(n(j - 1)).^(1 + m) + (n(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k1c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c(j - 1) - dt.*(1 + 2./dx.^2).*c(j) - dt.*alpha2.*(c(j + 1).^(1 + m0) - 2.*c(j).^(1 + m0) + c(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n(j)./(1 + gamma.*n(j));
            
            n1(j) = n(j) + dt.*k1n(j)./2;
            J1(j) = J(j) + dt.*k1j(j)./2;
            c1(j) = c(j) + dt.*k1c(j)./2;
            k2n(j) = dt.*J1(j);
            k2j(j) = -dt.*D2.*(n1(j - 2) + n1(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n1(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n1(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n1(j) + dt.*(r - 1./tau).*J1(j) - 2.*dt.*r.*n1(j).*J1(j)./beta0 - dt.*r.*n1(j).^2./(tau.*beta0) - dt.*chi0.*((n1(j + 1) - n1(j - 1)).*(c1(j + 1) - c1(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c1(j)).^2) - dt.*chi0.*n1(j).*(c1(j + 1) - 2.*c1(j) + c1(j - 1))./(tau.*dx.^2.*(1 + beta.*c1(j)).^2) + dt.*chi0.*beta.*n1(j).*(c1(j + 1) - c1(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c1(j)).^3) - dt.*alpha1.*((n1(j - 2)).^(1 + m) - 2.*(n1(j - 1)).^(1 + m) + (n1(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k2c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c1(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c1(j - 1) - dt.*(1 + 2./dx.^2).*c1(j) - dt.*alpha2.*(c1(j + 1).^(1 + m0) - 2.*c1(j).^(1 + m0) + c1(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n1(j)./(1 + gamma.*n1(j));
            
            n2(j) = n(j) + dt.*k2n(j)./2;
            J2(j) = J(j) + dt.*k2j(j)./2;
            c2(j) = c(j) + dt.*k2c(j)./2;
            k3n(j) = dt.*J2(j);
            k3j(j) = -dt.*D2.*(n2(j - 2) + n2(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n2(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n2(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n2(j) + dt.*(r - 1./tau).*J2(j) - 2.*dt.*r.*n2(j).*J2(j)./beta0 - dt.*r.*n2(j).^2./(tau.*beta0) - dt.*chi0.*((n2(j + 1) - n2(j - 1)).*(c2(j + 1) - c2(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c2(j)).^2) - dt.*chi0.*n2(j).*(c2(j + 1) - 2.*c2(j) + c2(j - 1))./(tau.*dx.^2.*(1 + beta.*c2(j)).^2) + dt.*chi0.*beta.*n2(j).*(c2(j + 1) - c2(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c2(j)).^3) - dt.*alpha1.*((n2(j - 2)).^(1 + m) - 2.*(n2(j - 1)).^(1 + m) + (n2(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k3c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c2(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c2(j - 1) - dt.*(1 + 2./dx.^2).*c2(j) - dt.*alpha2.*(c2(j + 1).^(1 + m0) - 2.*c2(j).^(1 + m0) + c2(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n2(j)./(1 + gamma.*n2(j));
            
            n3(j) = n(j) + dt.*k3n(j);
            J3(j) = J(j) + dt.*k3j(j);
            c3(j) = c(j) + dt.*k3c(j);
            k4n(j) = dt.*J3(j);
            k4j(j) = -dt.*D2.*(n3(j - 2) + n3(j + 2))./(tau.*dx.^4) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) - u./(2.*tau.*dx)).*n3(j + 1) + dt.*(alpha0./(tau.*dx.^2) + 4.*D2./(tau.*dx.^4) + u./(2.*tau.*dx)).*n3(j - 1) + dt.*(r./tau - 2.*alpha0./(tau.*dx.^2) - 6.*D2./(tau.*dx.^4)).*n3(j) + dt.*(r - 1./tau).*J3(j) - 2.*dt.*r.*n3(j).*J3(j)./beta0 - dt.*r.*n3(j).^2./(tau.*beta0) - dt.*chi0.*((n3(j + 1) - n3(j - 1)).*(c3(j + 1) - c3(j - 1)))./(4.*tau.*dx.^2.*(1 + beta.*c3(j)).^2) - dt.*chi0.*n3(j).*(c3(j + 1) - 2.*c3(j) + c3(j - 1))./(tau.*dx.^2.*(1 + beta.*c3(j)).^2) + dt.*chi0.*beta.*n3(j).*(c3(j + 1) - c3(j - 1)).^2./(2.*tau.*dx.^2.*(1 + beta.*c3(j)).^3) - dt.*alpha1.*((n3(j - 2)).^(1 + m) - 2.*(n3(j - 1)).^(1 + m) + (n3(j)).^(1 + m))./(tau.*dx.^2.*(1 + m));
            k4c(j) = dt.*(1./dx.^2 - 0.5.*u./dx).*c3(j + 1) + dt.*(1./dx.^2 + 0.5.*u./dx).*c3(j - 1) - dt.*(1 + 2./dx.^2).*c3(j) - dt.*alpha2.*(c3(j + 1).^(1 + m0) - 2.*c3(j).^(1 + m0) + c3(j - 1).^(1 + m0))./(dx.^2.*(1 + m0)) + dt.*beta0.*n3(j)./(1 + gamma.*n3(j));
            
        end
        
        n(j) = n(j) + dt.*(k1n(j) + 2.*k2n(j) + 2.*k3n(j) + k4n(j))./6;
        J(j) = J(j) + dt.*(k1j(j) + 2.*k2j(j) + 2.*k3j(j) + k4j(j))./6;
        c(j) = c(j) + dt.*(k1c(j) + 2.*k2c(j) + 2.*k3c(j) + k4c(j))./6;
        
        Nf(k,j) = n(j); Jf(k,j) = J(j); Cf(k,j) = c(j);
        Err(k,j) = abs(Nf(k,j) - nn(k,j));
%         Err(k,j) = Er(j);
        Err3(k,j) = Err(k,j);
    end
end

plot(x, Nf(round(3.*M./4),:), 'linewidth', 3);
xlabel 'x'; ylabel 'n(x,t = 0.75)'
% legend('\tau = 1', '\tau = 10^2', '\tau = 10^3')
% legend('D_2 = 0', 'D_2 = 10^2', 'D_2 = 10^3')
% legend('\alpha_1 = 0', '\alpha_1 = 0.005', '\alpha_1 = 1')
legend('m = 0', 'm = 10', 'm = 10^2')
% legend('u = 0', 'u = 10', 'u = 10^3')
% legend('r = 0', 'r = 0.01', 'r = 1')
end

figure
plot(x, Err1(round(3.*M./4),:), x, Err2(round(3.*M./4),:), x, Err3(round(3.*M./4),:), 'linewidth', 3); hold on
xlabel 'x'; ylabel 'Error(x,t = 0.75)'
% legend('\tau = 1', '\tau = 10^3', '\tau = 10^5')
% legend('D_2 = 0', 'D_2 = 10^2', 'D_2 = 10^3')
% legend('\alpha_1 = 0', '\alpha_1 = 0.005', '\alpha_1 = 1')
legend('m = 0', 'm = 10', 'm = 10^2')
% legend('u = 0', 'u = 10', 'u = 10^3')
% legend('r = 0', 'r = 0.01', 'r = 1')
toc
