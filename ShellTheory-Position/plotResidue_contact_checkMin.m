function [] = plotResidue_contact_checkMin(Np, zetac)
close all;
global N;
% global I;
N=round(Np)
% I=Np-N;

global H;


% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1000;
nu=0.5;
h=0.1;
minR=1*10^60;
hmin=0;
for Hh=0.99:-0.001:0.1
    H=Hh;
    [X] = initialGuess_contactRegion(Np, Hh, zetac);
    [r] = residuals_contactRegion(X);
    if minR > residueNorm(r, zetac)
        minR=residueNorm(r, zetac);
        hmin=Hh;
    end
end

fprintf('\n minR = %.4e \t hmin = %.4e \n', minR, hmin);
