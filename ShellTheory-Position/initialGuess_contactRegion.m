function[X] = initialGuess_contactRegion(Np, H, zetac)
fprintf('Started %s ... \t', mfilename);

N=round(Np);
X=zeros([N+2,1]);
fprintf('\n length(X) = %d \n', length(X));
% close all;
xc=H*tan(zetac);
phic=atan(xc/H);


for ii=1:N
    xpos(ii)=xc * (ii-1)/(N-1);
    X(ii,1)= atan(xpos(ii)/H);
    r(ii) = H / (cos(X(ii,1)));
%     X(N+ii,1)=H/cos(X(ii,1));
%     X(2*N+ii,1)=H;
    
    x1(ii) = r(ii) * sin(X(ii,1));
    z1(ii) = r(ii) * cos(X(ii,1));
end

if (phic - X(N,1)) < 1*10^(-5)
    fprintf('\n Correct final phi \n');
else
    fprintf('\n phix = %.2e \t final phi = %.2e \n', phic, X(N,1));
end
% for jj=1:I
%     X(3*N+jj,1)= zetac + (pi/2-zetac) * (jj-1)/(I-1);
%     psi=pi/2* (jj-1)/(I-1);
%     
%     X(3*N+I+jj, 1)=sqrt( (xc+H*sin(psi))^2+(H*cos(psi))^2 );
%     
%     x2(jj)=X(3*N+I+jj, 1)*sin(X(3*N+jj, 1));
%     z2(jj)=X(3*N+I+jj, 1)*cos(X(3*N+jj, 1));
% end

X(N+1,1)=zetac;
X(N+2,1)=1;
% plot(x1,z1, 'bs', 'MarkerFaceColor', 'b');
% hold on;
% plot(x2,z2, 'ro', 'MarkerFaceColor', 'r');
% xlim([0,1.2]);
% ylim([0,1]);