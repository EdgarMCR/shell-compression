function [] = plotResidue_contact(Np, Hh, zetac)
close all;
global N;
% global I;
N=round(Np)
% I=Np-N;

global H;
H=Hh;

% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1000;
nu=0.5;
h=0.1;

[X] = initialGuess_contactRegion(Np, Hh, zetac);
length(X)

for ii=1:N
    r(ii) = H / (cos(X(ii,1)));
    x1(ii)= r(ii)*sin(X(ii));
    z1(ii)= r(ii)*cos(X(ii));
%     fprintf('\n %d \t theta = %.2e \t r = %.2e \t z = %.2e \t  \n', ii, X(ii), X(N+ii), z(ii));
end
figure()
plot(x1,z1, 'bo', 'MarkerFaceColor','b');
hold on;
% plot(x2,z2, 'rs', 'MarkerFaceColor','r');
% scatter(x,z)
xlabel('x');
ylabel('z');

figure()
n=1:N;
% i=N+1:N+I;
plot(n, X(1:N), 'ro-', 'MarkerFaceColor','r');
hold on;
% plot(n, X(N+1:2*N), 'b-');
% plot(n, X(2*N+1:3*N), 'g');
% 
% plot(i, X(3*N+1:3*N+I), 'r-.');
% plot(i, X(3*N+I+1:3*N+2*I), 'b-.');

plot(N+1, X(N+1), 'ko', 'MarkerFaceColor','b');
plot(N+2, X(N+2), 'ks', 'MarkerFaceColor','g');

title('\phi (red) and r (blue)')
xlabel('Point')
ylabel('Value')

[r] = residuals_contactRegion(X);
r

figure()
n=1:N;
% i=N+1:N+I;
plot(n, r(1:N), 'ro-', 'MarkerFaceColor','r');
hold on;
% plot(n, r(N+1:2*N), 'bs-', 'MarkerFaceColor','b');
% plot(n, r(2*N+1:3*N), 'g<-', 'MarkerFaceColor','g');
% 
% plot(i, r(3*N+1:3*N+I), 'y>-', 'MarkerFaceColor','y');
% plot(i, r(3*N+I+1:3*N+2*I), 'ch-', 'MarkerFaceColor','c');

plot(N+1, r(N+1), 'ko', 'MarkerFaceColor','k');
plot(N+2, r(N+2), 'ks', 'MarkerFaceColor','k');
% ylim([-1000 1000])
title('\phi (red) and r (blue)')
xlabel('Point')
ylabel('Residue')

