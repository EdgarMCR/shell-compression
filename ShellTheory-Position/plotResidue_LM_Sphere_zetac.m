function [] = plotResidue_LM_Sphere_zetac(Np, Pp, zetacz)
close all;
global N;
global I;
N=round(Np/2);
I=Np-N;

global P;
P=Pp;

global zetac;
zetac=zetacz;

global lambda;
lambda=0;

global H;
H=0.9;
% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1;
nu=0.5;
h=0.05;

[X] = initialGuess_LM_Sphere(Np, zetac);
X(2*N+2*I+1,1)=zetac;
X(2*N+2*I+2,1)=Pp;

phi = X(1:N);
rs = X(N+1:2*N);

Phi = X(2*N+1:2*N+I);
R = X(2*N+I+1:2*N+2*I);
[v] = volumn_of_region(rs, phi, R, Phi, I, N);
fprintf('v=%f and difference is %f \n', v, abs(v-pi/4));
length(X);

for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
    x1(ii)= X(N+ii)*sin(X(ii));
    z1(ii)= X(N+ii)*cos(X(ii));
    x2(ii)= X(2*N+I+ii)*sin(X(2*N+ii));
    z2(ii)= X(2*N+I+ii)*cos(X(2*N+ii));
%     fprintf('\n %d \t theta = %.2e \t r = %.2e \t z = %.2e \t  \n', ii, X(ii), X(N+ii), z(ii));
end
figure()
plot(x1,z1, 'bo', 'MarkerFaceColor','b');
hold on;
plot(x2,z2, 'rs', 'MarkerFaceColor','r');
% scatter(x,z)
xlabel('x');
ylabel('z');

figure()
n=1:N;
i=N+1:N+I;
plot(n, X(1:N), 'ro-');
hold on;
plot(n, X(1*N+1:2*N), 'bs-.');

plot(i, X(2*N+1:2*N+I), 'rs-.');
plot(i, X(2*N+I+1:2*N+2*I), 'bs-.');

plot(2*N+I+1, X(2*N+2*I+1), 'ko');
plot(2*N+I+2, X(2*N+2*I+2), 'ks');

title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Value')

% [r] = residuals_LM_sphere(X);
[r]=residuals_LM_sphere_adpatedToCompressed_zetac(X);
% fprintf('r(2*N+2*I+1,1) = %f and r(2*N+2*I+2,1) = %f \n', r(2*N+2*I+1,1), r(2*N+2*I+2,1));

initialNorm=residueNorm(r, pi/2);
fprintf('\n Initial Residue Err = %.6e \n', initialNorm);

% area= H^2 * ( tan(X(N,1)) + pi/4) ;
% fprintf(' residue = %.4e \t area = %.4e \t residue p = %.4e \t residue zetac = %.4e  \n',  residueNorm(r, pi/2), area, r(N+2*I+1,1), r(N+2*I+2,1));

figure()
n=1:N;
i=N+1:N+I;
plot(n, r(1:N), 'ro-', 'MarkerFaceColor','r');
hold on;
plot(n, r(N+1:2*N), 'bs-', 'MarkerFaceColor','b');

plot(i, r(2*N+1:2*N+I), 'r>-', 'MarkerFaceColor','r');
plot(i, r(2*N+I+1:2*N+2*I), 'bh-', 'MarkerFaceColor','b');

plot(N+I+1, r(2*N+2*I+1), 'ko', 'MarkerFaceColor','k');
plot(N+I+2, r(N+2*I+2), 'ks', 'MarkerFaceColor','k');
% ylim([-1000 1000])
title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Residue')
legend('\zeta Contact', 'r Contact', '\zeta Free', 'r Free','\zeta_c', 'P','Location','southwest'); %, 'Pressure', '\zeta_c')
% minR=1*10^60;
% hmin=0;
% for Hh=0.999:-0.001:0.5
%     H=Hh;
%     [X] = initialGuess_compressedSphere3(Np, Hh, zetac);
%     [r] = residuals_compressedSphere(X);
%     if minR > residueNorm(r, pi/2)
%         minR=residueNorm(r, pi/2);
%         hmin=Hh;
%         area= H^2 * ( tan(X(N,1)) + pi/4) ;
%     end
% end
% 
% minR=1*10^60;
% zetacmin=0;
% for zetacT=1e-5:0.001:0.1
%     [X] = initialGuess_compressedSphere3(Np, Hh, zetacT);
%     [r] = residuals_compressedSphere(X);
%     if minR > residueNorm(r, pi/2)
%         minR=residueNorm(r, pi/2);
%         zetacmin=zetacT;
%         area= H^2 * ( tan(X(N,1)) + pi/4) ;
%     end
% end
% fprintf('\n minR = %.4e \t zetacmin = %.4e \t area = %.4e \n', minR, zetacmin, area);

% r(1:N)

