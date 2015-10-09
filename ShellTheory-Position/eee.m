close all;
clear all;
N=1000;
for ii=1:N
    zeta(ii)=pi/2 *(ii-1)/(N-1);
%     p(ii)=5*sin(2*zeta(ii));
    p(ii)=-cos(2*zeta(ii));
end
plot( zeta, p);