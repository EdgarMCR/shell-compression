function [] = plotResidue_split(Np, zetacs)
close all;
global N;
global I;
N=round(Np/2)
I=Np-N;

global zetac;
zetac=zetacs;

% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1000;
nu=0.5;
h=0.03;

[X] = initialGuess_splitSphere(Np, zetac);

length(X)

for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
    x1(ii)= X(N+ii)*sin(X(ii));
    z1(ii)= X(N+ii)*cos(X(ii));
    x2(ii)= X(3*N+I+ii)*sin(X(3*N+ii));
    z2(ii)= X(3*N+I+ii)*cos(X(3*N+ii));
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
plot(i, X(1*N+1:2*N), 'bs-.');

plot(i, X(3*N+1:3*N+I), 'rs-.');
plot(i, X(3*N+I+1:3*N+2*I), 'bs-.');

% plot(3*N+I+1, X(3*N+2*I+1), 'ko');
% plot(3*N+I+2, X(3*N+2*I+2), 'ks');

title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Value')

[r] = residuals_splitSphere(X);

% area= H^2 * ( tan(X(N,1)) + pi/4) ;
fprintf(' residue = %.4e  \n',  residueNorm(r, pi/2));

figure()
n=1:N;
i=N+1:N+I;
plot(n, r(1:N), 'ro-', 'MarkerFaceColor','r');
hold on;
plot(n, r(N+1:2*N), 'bs-', 'MarkerFaceColor','b');
plot(n, r(2*N+1:3*N), 'g<-', 'MarkerFaceColor','g');

plot(i, r(3*N+1:3*N+I), 'r>-', 'MarkerFaceColor','r');
plot(i, r(3*N+I+1:3*N+2*I), 'bh-', 'MarkerFaceColor','b');

% plot(N+I+1, r(N+2*I+1), 'ko', 'MarkerFaceColor','k');
% plot(N+I+2, r(N+2*I+2), 'ks', 'MarkerFaceColor','k');
% ylim([-1000 1000])
title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Residue')

