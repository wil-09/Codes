function RK4_LR2SS
% this code uses the Rkk4 method to integrate the long-range diffusion and
% long-range haptotaxis model in (1 + 1) system

clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
i = sqrt(-1); 

% les paramètres du model
D1 = 10^(-3); D2 = D1; alpha1 = D1; alpha2 = D1; gamma = D1; 
lambda = 0.12; s = 140; tau00 = 1.65; tau0 = tau00./s; k = -0.005;

% parameters of the solution derived
h1 = 0.09; h2 = 1; mu = 0.9; h001 = 0; h002 = 020; epsilon = 0.01;

% grids space and time
xlim = 20000.0; t0 = 0; tf = 1000; dt = 0.01; 
x0 = -xlim; xf = xlim; dx = 100; x = x0:dx:xf; t = t0:dt:tf;
N = numel(x); M = numel(t); h = dx;

n0 = 1; rho0 = 1;
K0 = 1 + lambda*n0^2; K10 = (1 + lambda*n0^2).^2; K30 = 1 - lambda.^2.*n0.^4;
iomega = ((k.^2.*K10.*(D1 + D2.*k.^2) + tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2)) - sqrt((k.^2.*K10.*(D1 + D2.*k.^2) + tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2)).^2 - 4.*K0.*(tau0.*k.^2.*(D1 + D2.*k.^2).*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2) - 4.*rho0.*tau0.*lambda.^2.*n0.^4.*k.^2.*(alpha1 + alpha2.*k.^2))))./(2.*K0);
Gamma1 = 4.*rho0.*tau0.*lambda.^2.*n0.^3./(tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2) - iomega.*K10);
Vg = 2.*i.*k.*(tau0.*K30.*(D1 + D2.*k.^2) + 2.*D2.*tau0.*k.^2.*(1 - gamma.*k.^2).*(K30 - n0.*k.^2.*K0) + 4.*rho0.*tau0.*lambda.^2.*n0.^4.*(alpha1 + alpha2.*k.^2) - iomega.*(tau0.*gamma.*K30 + tau0.*n0.*K0.*(1 - gamma.*k.^2) + K10.*(D1 + 2.*D2.*k.^2)))./(tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2) + k.^2.*K10.*(D1 + 2.*D2.*k.^2) - 2.*iomega.*K0);
a1 = 2.*k.^2.*(D1 + 2.*D2.*k.^2) - iomega; b1 = -n0.*k.^2.*(alpha1 + alpha2.*k.^2); c1 = (alpha1 + alpha2.*k.^2).*Gamma1;
a2 = -4.*lambda.^2.*rho0.*tau0.*n0.^3; b2 = tau0.*(1 - 4.*gamma.*k.^2).*(K30 - 4.*n0.*K0.*k.^2) - 2.*iomega.*k.^2;
c2 = 6.*iomega.*lambda.*n0.*Gamma1.*(1 + lambda.*n0.^2).^2 + 2.*rho0.*n0.*tau0.*lambda.*(3.*lambda.*n0 - k.^2.*(3 - lambda.*n0.^2)) + tau0.*Gamma1.*(1 - gamma.*k.^2).*(k.^2.*(1 + lambda.*n0.^2).*(3 - lambda.*n0.^2) + 4.*n0.^2.*(k.^2.*(1 + lambda.*n0.^2) + lambda.*n0));
Gamma21 = Gamma1.*(Vg.*K10 - 2.*i.*k.*tau0.*(K30.*gamma + n0.*K0.*(1 - 2.*gamma.*k.^2)))./(tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2) - iomega.*K0);
B2 = (b2.*c1 - b1.*c2)./(a1.*b2 - a2.*b1); Gamma2 = (a1.*c2 - a2.*c1)./(a1.*b2 - a2.*b1);
B = i.*k.*(alpha1 + alpha2.*k.^2).*(conj(Gamma1) - Gamma1)./Vg;

Pro = Gamma21.*(K10.*Vg - 2.*i.*k.*tau0.*(gamma.*K30 + n0.*K0.*(1 - 2.*gamma.*k.^2))) - Gamma1.*tau0.*(gamma.*K30 + n0.*K0.*(1 - gamma.*k.^2));
Qro = Gamma1.*B.*(tau0.*(1 - gamma.*k.^2).*(4.*lambda.^2.*n0.^3 + k.^2.*(1 + lambda.*n0.^2).*(1 + 5.*lambda.*n0.^2)) + 6.*iomega.*lambda.*n0.*(1 + lambda.*n0.^2).^2) + Gamma1.*(6.*iomega.*lambda.*n0.*(1 + lambda.*n0.^2).^2 + 4.*tau0.*n0.*(1 - gamma.^2).*(3.*lambda.^2.*n0.^3 + 2.*lambda.*k.^2.*(3 + 2.*lambda.*n0.^2))) + conj(Gamma1).*(2.*lambda.^2.*n0.^2.*tau0.*(1 - gamma.*k.^2).*(3 + 10.*n0.*k.^2) - 6.*iomega.*lambda.*n0.*(1 + lambda.*n0.^2).*(1 + 5.*lambda.*n0.^2) )+ B2.*(16.*iomega.*lambda.*n0.*(1 + lambda.*n0.^2).^2 + 4.*lambda.*rho0.*n0.*tau0.*(3.*lambda.*n0 + 2.*k.^2.*(3 - lambda.*n0.^2))) + B2.*conj(Gamma1).*(-6.*iomega.*lambda.*n0.*(1 + lambda.*n0.^2).^2 + tau0.*(4.*lambda.^2.*n0.^3 + 3.*k.^2.*(1 - gamma.*k.^2).*(3.*lambda.*n0.^2 - 1))) + 4.*lambda.*n0.^2.*tau0.*(lambda.*n0 + 12.*k.^2.*(1 + lambda.*n0.^2));
Pn = D1 + 6.*D2.*k.^2 - n0.*Gamma1.*(alpha2 + 6.*alpha2.*k.^2) - 2.*i.*k.*n0.*(alpha1 + 2.*alpha2.*k.^2).*Gamma21;
Qn = k.^2.*B.*Gamma1.*(alpha1 + alpha2.*k.^2) + 2.*k.^2.*Gamma2.*(alpha1 + 4.*alpha2.*k.^2) - k.^2.*B2.*conj(Gamma1).*(alpha1 + alpha2.*k.^2);

P1 = -((-Pro.*n0.*k.^2.*(alpha1 + alpha2.*k.^2) - Pn.*(tau0.*(1 - gamma.*k.^2).*(K30 - n0.*k.^2.*K0) - iomega.*K10)));
Q1 = -(-Qro.*n0.*k.^2.*(alpha1 + alpha2.*k.^2) - Qn.*(iomega.*K10 + tau0.*(1 - gamma.*k.^2).*(K30 - n0.*K0.*k.^2)));

P = P1./(k.^2.*(D1 + D2.*k.^2).*(1 + lambda.*n0.^2).^2 - tau0.*(1 - gamma.*k.^2).*(1 - lambda.*n0.^2 - n0.*k.^2.*(1 + lambda.*n0.^2)));
Q = Q1./(k.^2.*(D1 + D2.*k.^2).*(1 + lambda.*n0.^2).^2 - tau0.*(1 - gamma.*k.^2).*(1 - lambda.*n0.^2 - n0.*k.^2.*(1 + lambda.*n0.^2)));

% parameters of thediscretized model
c1n = -D1./(12.*h.^2) - D2./h.^4; c2n = 16.*D1./(12.*h.^2) + 4.*D2./h.^4;
c3n = -30.*D1./(12.*h.^2) - 6.*D2./h.^4; c4n = -alpha1./(12.*h.^2) - alpha2./h.^4;
c5n = 16.*alpha1./(12.*h.^2) + 4.*alpha2./h.^4; 
c6n = -30.*alpha1./(12.*h.^2) - 6.*alpha2./h.^4; c7n = -1./(12.*h); 
c8n = alpha1./(12.*h) + alpha2./(2.*h.^3);
c9n = -8.*alpha1./(12.*h) - 2.*alpha2./(2.*h.^3);

d1r = 2.*tau0.*lambda./(12.*h); d2r = 1 - 30.*gamma./(12.*h.^2); d3r = 1./(12.*h.^2);
d4r = -tau0./(12.*h.^2); d5r = -2.*tau0./(12.*h); d6r = gamma./(2.*h.^3) - 1./(12.*h);
d7r = 8./(12.*h) - 2.*gamma./(2.*h.^3); d8r = gamma./h.^4 - 1./(12.*h.^2); 
d9r = 16./(12.*h.^2) - 4.*gamma./h.^4; d10r = 6.*gamma./h.^4 - 30./(12.*h.^2);

nn = zeros(1,N); Rho = zeros(1,N); varphi = @(t) -log(1 + sin(iomega.*epsilon.^2.*t));
for kk = 1:M
    for j = 1:N
%         varphi(kk) = -log(1 + sin(iomega.*epsilon.^2.*t(kk)));
        nn(j,kk) = n0 + epsilon.*real((mu.*(exp(-(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*(exp((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + exp((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002)) + mu.^2.*(Q./P).*exp(-3.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*((h1 - h2)./(h1 + h2)).^2.*(exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./h1.^2 + exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002) + ((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001))./h2.^2))./(1 + mu.^2.*(Q./P).*exp(-2.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*(exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001))./(8.*h1.^2) + exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./(8.*h2.^2) + exp(((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./(h1 + h2).^2) + mu.^4.*(Q./(8.*h1.*h2.*P)).^2.*exp(-4.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*((h1 - h2)./(h1 + h2)).^4.*exp(2.*(((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))))).*exp(i.*k.*x(j) - iomega.*t(kk)));
        n0_in(j,kk) = nn(j,kk);

        Rho(j,kk) = rho0 + epsilon.*real(Gamma1.*(mu.*(exp(-(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*(exp((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + exp((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002)) + mu.^2.*(Q./P).*exp(-3.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*((h1 - h2)./(h1 + h2)).^2.*(exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./h1.^2 + exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002) + ((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001))./h2.^2))./(1 + mu.^2.*(Q./P).*exp(-2.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*(exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001))./(8.*h1.^2) + exp(2.*((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./(8.*h2.^2) + exp(((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))./(h1 + h2).^2) + mu.^4.*(Q./(8.*h1.*h2.*P)).^2.*exp(-4.*(varphi(t(kk)) - 4.*h1.*h2.*(epsilon.^2.*t(kk)).*P.*(h1 + h2)./(h1 - h2))).*((h1 - h2)./(h1 + h2)).^4.*exp(2.*(((epsilon.*(x(j) - Vg.*t(kk))).*h1 + (epsilon.^2.*t(kk)).*P.*h1.*(h1.^2 + 4.*h2.^2 + 3.*h1.*h2)./(h2 - h1) + h001) + ((epsilon.*(x(j) - Vg.*t(kk))).*h2 + (epsilon.^2.*t(kk)).*P.*h2.*(4.*h1.^2 - h2.^2 + 5.*h1.*h2)./(h2 - h1) + h002))))).*exp(i.*k.*x(j) - iomega.*t(kk)));
        rho0_in(j,kk) = Rho(j,kk);
    end
end

n(1:N) = nn(1:N,1) + 0.01.*rand.*(max(nn(1:N,1)) - min(nn(1:N,1)))./max(nn(1:N,1)); 
% rho_in(1:N) = Rho(1:N,1); 
ro(1:N) = Rho(1:N,1) + 0.01.*rand.*(max(Rho(1:N,1)) - min(Rho(1:N,1)))./max(Rho(1:N,1));

% figure
% subplot(1,2,1)
% plot(x, n, 'linewidth',3)
% xlabel 'x'; ylabel 'n(x,0)'
% 
% subplot(1,2,2)
% plot(x, ro,'linewidth',3)
% xlabel 'x'; ylabel '\rho(x,0)'

kk1n = zeros(1,N); kk2n = zeros(1,N); kk3n = zeros(1,N); kk4n = zeros(1,N);
kk1ro = zeros(1,N); kk2ro = zeros(1,N); kk3ro = zeros(1,N); kk4ro = zeros(1,N);
n1= zeros(1,N); n2= zeros(1,N); n3= zeros(1,N);
ro1= zeros(1,N); ro2= zeros(1,N); ro3= zeros(1,N);

Nf = []; Rof = [];

% Periodic boundary conditions on the domain.
for kk = 1:M
    for j = 1:N
        if(j == 1)
            jp = j + 1; jp2 = jp;
%             jp2 = j + 2;
            
            jm = N; jm2 = jm;
%             jm2 = 1;
        elseif(j == N)
            jp = 1; jp2 = jp;
%             jp2 = 2;
            
            jm = j - 1; jm2 = jm;
%             jm2 = j - 2; 
        else
            jp = j + 1; jp2 = jp;
%             jp2 = j + 2;
            
            jm = j - 1; jm2 = jm;
%             jm2 = j - 2;
%             if(jm2 < 0);
%                 jm2 = 1;
%             end
        end     
%         kk1n(j) = n(jp2) + n(jm2);
        kk1n(j) = c1n.*(n(jp2) + n(jm2)) + c2n.*(n(jp) + n(jm)) + ...
            c3n.*n(j) - n(j).*(c4n.*(ro(jp2) + ro(jm2)) + ...
            c5n.*(ro(jp) + ro(jm)) + c6n.*ro(j)) + c7n.*(n(jm2) -...
            n(jp2) + 8.*(n(jp) - n(jm))).*(c8n.*(ro(jm2) - ro(jp2)) +...
            c9n.*(ro(jp) - ro(jm)));
        
        kk1ro(j) = d1r.*(1 - lambda.*n(j).^2).*(d2r.*ro(j) - d3r.*(ro(jp2) +...
            ro(jm2) - 16.*(ro(jp) + ro(jm)))).*(n(jm2) - n(jp2) +...
            8.*(n(jp) - n(jm))).^2./(1 + lambda.*n(j).^2).^3 + ...
            d4r.*(d2r.*ro(j) - d3r.*(ro(jp2) + ro(jm2) - 16.*(ro(jp) +...
            ro(jm)))).*(1 - lambda.*n(j).^2).*(16.*(ro(jp) + ro(jm)) - ...
            30.*n(j) - (ro(jp2) + ro(jm2)))./(1 + lambda.*n(j).^2).^2 +...
            d5r.*(d6r.*(ro(jp2) - ro(jm2)) + d7r.*(ro(jp) - ...
            ro(jm))).*(n(jm2) - n(jp2) + 8.*(n(jp) - n(jm))).*(1 - ...
            lambda.*n(j).^2)./(1 + lambda.*n(j).^2).^2 - ...
            tau0.*n(j).*(d8r.*(ro(jp2) + ro(jm2)) + d9r.*(ro(jp) + ...
            ro(jm)) + d10r.*ro(j))./(1 + lambda.*n(j).^2);
                
        n1(j) = n(j) + dt.*kk1n(j)./2;
        ro1(j) = ro(j) + dt.*kk1ro(j)./2;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
        kk2n(j) = c1n.*(n1(jp2) + n1(jm2)) + c2n.*(n1(jp) + n1(jm)) + ...
            c3n.*n1(j) - n1(j).*(c4n.*(ro1(jp2) + ro1(jm2)) + ...
            c5n.*(ro1(jp) + ro1(jm)) + c6n.*ro1(j)) + c7n.*(n1(jm2) -...
            n1(jp2) + 8.*(n1(jp) - n1(jm))).*(c8n.*(ro1(jm2) - ro1(jp2)) +...
            c9n.*(ro1(jp) - ro1(jm)));
        
        kk2ro(j) = d1r.*(1 - lambda.*n1(j).^2).*(d2r.*ro1(j) - d3r.*(ro1(jp2) +...
            ro1(jm2) - 16.*(ro1(jp) + ro1(jm)))).*(n1(jm2) - n1(jp2) +...
            8.*(n1(jp) - n1(jm))).^2./(1 + lambda.*n1(j).^2).^3 + ...
            d4r.*(d2r.*ro1(j) - d3r.*(ro1(jp2) + ro1(jm2) - 16.*(ro1(jp) +...
            ro1(jm)))).*(1 - lambda.*n1(j).^2).*(16.*(ro1(jp) + ro1(jm)) - ...
            30.*n1(j) - (ro1(jp2) + ro1(jm2)))./(1 + lambda.*n1(j).^2).^2 +...
            d5r.*(d6r.*(ro1(jp2) - ro1(jm2)) + d7r.*(ro1(jp) - ...
            ro1(jm))).*(n1(jm2) - n1(jp2) + 8.*(n1(jp) - n1(jm))).*(1 - ...
            lambda.*n1(j).^2)./(1 + lambda.*n1(j).^2).^2 - ...
            tau0.*n1(j).*(d8r.*(ro1(jp2) + ro1(jm2)) + d9r.*(ro1(jp) + ...
            ro1(jm)) + d10r.*ro1(j))./(1 + lambda.*n1(j).^2);
        
        n2(j) = n(j) + dt.*kk2n(j)./2;
        ro2(j) = ro(j) + dt.*kk2ro(j)./2;
%--------------------------------------------------------------------------     
%--------------------------------------------------------------------------
        kk3n(j) = c1n.*(n2(jp2) + n2(jm2)) + c2n.*(n2(jp) + n2(jm)) + ...
            c3n.*n2(j) - n2(j).*(c4n.*(ro2(jp2) + ro2(jm2)) + ...
            c5n.*(ro2(jp) + ro2(jm)) + c6n.*ro2(j)) + c7n.*(n2(jm2) -...
            n2(jp2) + 8.*(n2(jp) - n2(jm))).*(c8n.*(ro2(jm2) - ro2(jp2)) +...
            c9n.*(ro2(jp) - ro2(jm)));
        
        kk3ro(j) = d1r.*(1 - lambda.*n2(j).^2).*(d2r.*ro2(j) - d3r.*(ro2(jp2) +...
            ro2(jm2) - 16.*(ro2(jp) + ro2(jm)))).*(n2(jm2) - n2(jp2) +...
            8.*(n2(jp) - n2(jm))).^2./(1 + lambda.*n2(j).^2).^3 + ...
            d4r.*(d2r.*ro2(j) - d3r.*(ro2(jp2) + ro2(jm2) - 16.*(ro2(jp) +...
            ro2(jm)))).*(1 - lambda.*n2(j).^2).*(16.*(ro2(jp) + ro2(jm)) - ...
            30.*n2(j) - (ro2(jp2) + ro2(jm2)))./(1 + lambda.*n2(j).^2).^2 +...
            d5r.*(d6r.*(ro2(jp2) - ro2(jm2)) + d7r.*(ro2(jp) - ...
            ro2(jm))).*(n2(jm2) - n2(jp2) + 8.*(n2(jp) - n2(jm))).*(1 - ...
            lambda.*n2(j).^2)./(1 + lambda.*n2(j).^2).^2 - ...
            tau0.*n2(j).*(d8r.*(ro2(jp2) + ro2(jm2)) + d9r.*(ro2(jp) + ...
            ro2(jm)) + d10r.*ro2(j))./(1 + lambda.*n2(j).^2);
        
        n3(j) = n(j) + dt.*kk3n(j);
        ro3(j) = ro(j) + dt.*kk3ro(j);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
        kk4n(j) = c1n.*(n3(jp2) + n3(jm2)) + c2n.*(n3(jp) + n3(jm)) + ...
            c3n.*n3(j) - n3(j).*(c4n.*(ro3(jp2) + ro3(jm2)) + ...
            c5n.*(ro3(jp) + ro3(jm)) + c6n.*ro3(j)) + c7n.*(n3(jm2) -...
            n3(jp2) + 8.*(n3(jp) - n3(jm))).*(c8n.*(ro3(jm2) - ro3(jp2)) +...
            c9n.*(ro3(jp) - ro3(jm)));
        
        kk4ro(j) = d1r.*(1 - lambda.*n3(j).^2).*(d2r.*ro3(j) - d3r.*(ro3(jp2) +...
            ro3(jm2) - 16.*(ro3(jp) + ro3(jm)))).*(n3(jm2) - n3(jp2) +...
            8.*(n3(jp) - n3(jm))).^2./(1 + lambda.*n3(j).^2).^3 + ...
            d4r.*(d2r.*ro3(j) - d3r.*(ro3(jp2) + ro3(jm2) - 16.*(ro3(jp) +...
            ro3(jm)))).*(1 - lambda.*n3(j).^2).*(16.*(ro3(jp) + ro3(jm)) - ...
            30.*n3(j) - (ro3(jp2) + ro3(jm2)))./(1 + lambda.*n3(j).^2).^2 +...
            d5r.*(d6r.*(ro3(jp2) - ro3(jm2)) + d7r.*(ro3(jp) - ...
            ro3(jm))).*(n3(jm2) - n3(jp2) + 8.*(n3(jp) - n3(jm))).*(1 - ...
            lambda.*n3(j).^2)./(1 + lambda.*n3(j).^2).^2 - ...
            tau0.*n3(j).*(d8r.*(ro3(jp2) + ro3(jm2)) + d9r.*(ro3(jp) + ...
            ro3(jm)) + d10r.*ro3(j))./(1 + lambda.*n3(j).^2);
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
        
        n(j) = n(j) + dt.*(kk1n(j) + 2.*kk2n(j) + 2.*kk3n(j) + kk4n(j))./6;
        ro(j) = ro(j) + dt.*(kk1ro(j) + 2.*kk2ro(j) + 2.*kk3ro(j) + kk4ro(j))./6;
        
        Nf(j,kk) = n(j); Rof(j,kk) = ro(j);
    end
end

figure
plot(x, n0_in(:,25000), 'g', x, Nf(:,25000), 'r', 'linewidth', 3); hold on
xlabel 'x'; ylabel 'n(x,t)'
legend('Analytical Solutions', 'Numerical solutions')

figure
plot(x, n0_in(:,50000), 'g', x, Nf(:,50000), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('Analytical Solutions', 'Numerical solutions')

figure
plot(x, n0_in(:,75000), 'g', x, Nf(:,75000), 'r', 'linewidth', 3); hold on
xlabel 'x'; ylabel 'n(x,t)'
legend('Analytical Solutions', 'Numerical solutions')

figure
plot(x, n0_in(:,99000), 'g', x, Nf(:,9900), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('Analytical Solutions', 'Numerical solutions')

% subplot(1,2,2)
% plot(x, Rho(:,1), 'k', 'linewidth', 3); hold on
% plot(x, Rof(:,1), 'g', x, Rof(:,50), 'r', x, Rof(:,500), 'b', x, Rof(:,990), 'c', 'linewidth', 3)
% xlabel 'x'; ylabel '\rho(x,t)'
% % legend( 't = dt', 't = 50dt', 't = 100dt', 't = 250*dt')
% 
figure
% subplot(1,2,1)
mesh(x,t,Nf'); colormap(jet)
xlabel('x'); ylabel('t'); zlabel('n(x,t)')
% 
% subplot(1,2,2)
% mesh(x,t,Rof'); colormap(jet)
% xlabel('x'); ylabel('t'); zlabel('\rho(x,t)')