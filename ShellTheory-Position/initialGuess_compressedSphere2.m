function[X] = initialGuess_compressedSphere2(Np, H, zetac)
fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([N+2*I+2,1]);
fprintf('\n length(X) = %d \n', length(X));
% close all;

xc=H*tan(zetac)

for ii=1:N
    xpos(ii)=xc * (ii-1)/(N-1);
    X(ii,1)= atan(xpos(ii)/H);    
    x1(ii) = sqrt(H^2+xpos(ii)^2) * sin(X(ii,1));
    z1(ii) = sqrt(H^2+xpos(ii)^2) * cos(X(ii,1));
end

for jj=1:I
    
    psi=pi/2* (jj-1)/(I-1);
    
    X(N+I+jj, 1)=sqrt( (xc+H*sin(psi))^2+(H*cos(psi))^2 );
    
    x2(jj)=X(N+I+jj, 1)*sin(X(N+jj, 1));
    z2(jj)=X(N+I+jj, 1)*cos(X(N+jj, 1));
    
    X(N+jj,1)= atan(x2(jj)/z2(jj));
end

X(N+2*I+1,1)=zetac;
X(N+2*I+2,1)=1;
% plot(x1,z1, 'bs', 'MarkerFaceColor', 'b');
% hold on;
% plot(x2,z2, 'ro', 'MarkerFaceColor', 'r');
% xlim([0,1.2]);
% ylim([0,1]);