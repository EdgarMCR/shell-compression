function[X, p] = initialGuess_closeToSpherical(Np, H)
% fprintf('Started %s \t', mfilename);
N=round(Np/2);
I=Np-N;

X=zeros([N+2*I+1,1]);

zetac=pi/2*H;
%assume AE=1
% p=(H^2*pi^3)/(16*zetac^3); %What was that?
AE=1;
p=(AE*pi)/(4*zetac)  * ((H^2 * pi^2)/(4*zetac^2) - 1);
fprintf('zetac= %.4f  p=%.2e\n',zetac, p);
x_c=1-H;
for ii=1:N
    X(ii,1)=(ii-1)*(x_c)/(N-1);
end

for ii=1:I
    X(N+ii,1)=x_c+ H*cos(pi/2*(1-(ii-1)/(I-1)));
    X(N+I+ii,1)=H*sin(pi/2*(1-(ii-1)/(I-1)));
end
X(N+2*I+1,1)=zetac;

close all;
h=figure();set(gcf,'Visible', 'off'); 
plot(X(1:N,1), H, 'ro')
hold on
plot(X(N+1:N+2, 1), X(N+I+1:N+I+2, 1), 'gs')
plot(X(N+2:N+I, 1), X(N+I+2:N+2*I, 1), 'bs')
filename=sprintf('pics\\H=%.2f_p=%.1e.png',H, p);
print(h,filename, '-dpng')
% fprintf('ended %s \n', mfilename);