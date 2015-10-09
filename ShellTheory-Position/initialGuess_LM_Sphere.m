function[X] = initialGuess_LM_Sphere(Np, zetac)
% fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([2*N+2*I,1]);
smallPerturbation=0 ;1e-7;
for ii=1:N
    zeta=(ii-1)/(N-1) * zetac;
    X(ii,1)=zeta + (-1)^ii * smallPerturbation;
    X(N+ii,1)=1 - (-1)^ii * smallPerturbation;
    
    xpos(ii)=X(N+ii,1)* sin ( X(ii,1) );
    zpos(ii)=X(N+ii,1)* cos ( X(ii,1) );
end

for ii=1:I
    zeta=zetac + (ii-1)/(N-1) * (pi/2 - zetac);

    X(2*N+ii,1) = zeta;
    X(2*N+I+ii,1) = 1;
    
    Xpos(ii)=X(2*N+I+ii,1)* sin ( X(2*N+ii,1) );
    Zpos(ii)=X(2*N+I+ii,1)* cos ( X(2*N+ii,1) );
end



% figure()
% plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
% hold on;
% plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
% xlim([0,1.2]);
% ylim([0,1]);
end

