function[X] = initialGuess_compressedSphere(Np, H, zetac)
fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([3*N+2*I+2,1]);
fprintf('\n length(X) = %d \n', length(X));

for ii=1:N
    thetha=zetac * (ii-1)/(N-1);
    X(ii,1)=thetha;
    X(N+ii,1)=H/cos(thetha);
    X(2*N+ii,1)=H;
    z(ii)=X(N+ii)*cos(X(ii));
end
z
xc=H*tan(zetac)
for jj=1:I
    thetha2 = zetac + (pi - zetac) * (jj-1)/(I-1);
    psi=pi/2 * (jj-1)/(I-1);
    
    X(3*N+ii,1)=thetha2;
    X(3*N+I+ii,1)=sqrt( (xc+H*sin(psi))^2 + (H*cos(psi))^2);
end
X(3*N+2*I+1,1)=zetac;
X(3*N+2*I+2,1)=1;
fprintf('... Finished %s \n', mfilename);

for ii=1:N
    x(ii)= X(N+ii)*sin(X(ii));
    z(ii)= X(N+ii)*cos(X(ii));
%     fprintf('\n %d \t theta = %.2e \t r = %.2e \t z = %.2e \t  \n', ii, X(ii), X(N+ii), z(ii));
    Xx(ii)= X(3*N+I+ii)*sin(X(3*N+ii));
    Zz(ii)= X(3*N+I+ii)*cos(X(3*N+ii)); 
end


figure()
plot(x,z, 'bo', 'MarkerFaceColor','b');
hold on;
plot(Xx,Zz, 'rs', 'MarkerFaceColor','r');
% scatter(x,z)
xlabel('x');
ylabel('z');
    

% close all;
% h=figure(); %set(gcf,'Visible', 'off'); 
% plot(X(1:N,1), X(N+1:2*N,1), 'ro')

% hold on
% plot(X(N+1:N+9, 1), X(N+I+1:N+I+9, 1), 'gs')
% plot(X(N+10:N+I, 1), X(N+I+10:N+2*I, 1), 'bs')
% filename=sprintf('picss\\H=%.2f_p=%.1e.png',H);
% print(h,filename, '-dpng')
% fprintf('ended %s \n', mfilename);