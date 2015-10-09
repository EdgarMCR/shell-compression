function [] = plotResidue_InflatingSphere(Np)
close all;
global N;
N=Np;
% global I;
% N=round(Np/2)
% I=Np-N;

% global zetac;
% zetac=zetacs;

% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness
global P;
P=1;

E=1000;
nu=0.5;
h=0.03;

[X] = initialGuess_inflatedSphere(N);

length(X)

for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
    x1(ii)= X(N+ii)*sin(X(ii));
    z1(ii)= X(N+ii)*cos(X(ii));
end
figure()
plot(x1,z1, 'bo', 'MarkerFaceColor','b');
% scatter(x,z)
xlabel('x');
ylabel('z');

figure()
n=1:N;

plot(n, X(1:N), 'ro-');
hold on;
plot(n, X(1*N+1:2*N), 'bs-.');


title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Value')

[r] = residuals_inflatingSphere_nonDim(X);

% area= H^2 * ( tan(X(N,1)) + pi/4) ;
fprintf(' residue = %.4e  \n',  residueNorm(r, pi/2));

figure()
n=1:N;
plot(n, r(1:N), 'ro-', 'MarkerFaceColor','r');
hold on;
plot(n, r(N+1:2*N), 'bs-', 'MarkerFaceColor','b');
% plot(n, r(2*N+1:3*N), 'g<-', 'MarkerFaceColor','g');
% 
% plot(i, r(3*N+1:3*N+I), 'r>-', 'MarkerFaceColor','r');
% plot(i, r(3*N+I+1:3*N+2*I), 'bh-', 'MarkerFaceColor','b');

% plot(N+I+1, r(N+2*I+1), 'ko', 'MarkerFaceColor','k');
% plot(N+I+2, r(N+2*I+2), 'ks', 'MarkerFaceColor','k');
% ylim([-1000 1000])
title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Residue')

