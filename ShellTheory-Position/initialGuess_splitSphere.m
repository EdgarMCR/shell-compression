function[X] = initialGuess_splitSphere(Np, zetac)
% fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([3*N+2*I,1]);


for ii=1:N
    X(ii,1)=zetac * (ii-1)/(N-1);
    X(N+ii,1)=1;
    X(2*N+ii,1)=2;
    
    xpos(ii)=X(N+ii,1)*sin(X(ii,1));
    zpos(ii)=X(N+ii,1)*cos(X(ii,1));
end

for ii=1:I
    X(3*N+ii,1)=zetac + (pi/2 - zetac) * (ii-1)/(N-1);
    X(3*N+I+ii,1)=1;
    
    Xpos(ii)=X(3*N+I+ii,1)*sin(X(3*N+ii,1));
    Zpos(ii)=X(3*N+I+ii,1)*cos(X(3*N+ii,1));
end


figure()
plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
hold on;
plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
xlim([0,1.2]);
ylim([0,1]);