function[X] = initialGuess_LM_Sphere(Np, zetac)
% fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([2*N+2*I,1]);
smallPerturbation=1e-7;
for ii=1:N
    zeta(ii)=(ii-1)/(N-1) * zetac;
    X(ii,1)= (-1)^ii * smallPerturbation;
    X(N+ii,1)=(-1)^ii * smallPerturbation;
    
    xpos(ii)=(1+X(N+ii,1))* sin ( (zeta(ii) + X(ii,1)) );
    zpos(ii)=(1+X(N+ii,1))* cos ( (zeta(ii) + X(ii,1)) );
end

for ii=1:I
    zeta(N+ii)=zetac + (ii-1)/(N-1) * (pi/2 - zetac);

    X(2*N+ii,1) = 0;
    X(2*N+I+ii,1) = 0;
    
    Xpos(ii)=(1+X(2*N+I+ii,1))* sin ( (zeta(N+ii) +X(2*N+ii,1)) );
    Zpos(ii)=(1+X(2*N+I+ii,1))* cos ( (zeta(N+ii) +X(2*N+ii,1)) );
end



figure()
plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
hold on;
plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
xlim([0,1.1]);
ylim([0,1.1]);
end

