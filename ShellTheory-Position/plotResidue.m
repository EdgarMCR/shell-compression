function [] = plotResidue(Np)
close all;
global P;
global lambda; %first lame parameter
global mu; %second lame parameter
global h; %membran thickness

lambda=1;
mu=1;
h=1;

P=0;

global N;
N=Np;

[X] = initialGuess_inflatedSphere(Np);



figure()
plot(X(1:N), 'r');
hold on;
plot(X(N+1:2*N), 'b');
title('\theta (red) and r (blue)')
xlabel('Point')
ylabel('Value')

[r] = residuals_inflatingSphere(X);
length(r);

figure()
plot(r(1:N), 'r')
hold on;
plot(r(N+1:2*N), 'b')
title('Residue \theta (red) and r (blue)')
xlabel('Point')
ylabel('Residue')