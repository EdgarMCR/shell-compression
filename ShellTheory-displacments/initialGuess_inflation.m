function[X] = initialGuess_inflation(Np)
% fprintf('Started %s ... \t', mfilename);

N=round(Np);
X=zeros([2*N,1]);
smallPerturbation=1e-4;
for ii=1:N
    zeta(ii)=(ii-1)/(N-1) * pi/2;
    X(ii,1)= (-1)^ii * smallPerturbation;
    X(N+ii,1)=(-1)^ii * smallPerturbation;
    
    xpos(ii)=(1+X(N+ii,1))* sin ( (zeta(ii) + X(ii,1)) );
    zpos(ii)=(1+X(N+ii,1))* cos ( (zeta(ii) + X(ii,1)) );
end

figure()
plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
hold on;
% plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
xlim([0,1.1]);
ylim([0,1.1]);
end

