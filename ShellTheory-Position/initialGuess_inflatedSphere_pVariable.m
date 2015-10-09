function[X] = initialGuess_inflatedSphere_pVariable(Np)
% fprintf('Started %s \t', mfilename);

N=Np;
X=zeros([2*Np+1,1]);

for ii=1:N
    thetha=pi/2 * (ii-1)/(N-1);
    X(ii,1)=thetha;
    X(N+ii,1)=1;
end
X(2*N+1,1)=1;

% close all;
% h=figure(); %set(gcf,'Visible', 'off'); 
% plot(X(1:N,1), X(N+1:2*N,1), 'ro')

% hold on
% plot(X(N+1:N+9, 1), X(N+I+1:N+I+9, 1), 'gs')
% plot(X(N+10:N+I, 1), X(N+I+10:N+2*I, 1), 'bs')
% filename=sprintf('picss\\H=%.2f_p=%.1e.png',H);
% print(h,filename, '-dpng')
% fprintf('ended %s \n', mfilename);