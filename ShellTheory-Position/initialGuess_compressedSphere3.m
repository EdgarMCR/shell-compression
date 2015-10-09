function[X] = initialGuess_compressedSphere3(Np, H, zetac)
% fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([2*N+2*I+2,1]);
% fprintf('\n length(X) = %d \n', length(X));
% close all;


% x_c= ( zetac);
x_c = H* zetac * (pi/2)/(pi/2 - zetac);
for ii=1:N
    xpos(ii)=(ii-1)*(x_c)/(N-1);
    zpos(ii)=H;
    X(ii,1)=atan(xpos(ii)/H);
    X(N+ii,1)=sqrt( xpos(ii)^2 + H^2);
%     X(2*N+ii,1)=1;
end

for ii=1:I
    Xpos(ii)=x_c+H*cos(pi/2*(1-(ii-1)/(I-1)));
    Zpos(ii)=H*sin(pi/2*(1-(ii-1)/(I-1)));
    X(2*N+ii,1) = atan(Xpos(ii)/Zpos(ii));
    X(2*N+I+ii,1) = sqrt( Xpos(ii)^2+Zpos(ii)^2 );
end

X(2*N+2*I+1,1)=0;
X(2*N+2*I+2,1)=zetac;


end

% figure()
% plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
% hold on;
% plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
% xlim([0,1.2]);
% ylim([0,1]);