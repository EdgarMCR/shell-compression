function [normR] = testResidue(Np)
%Test Residue

[X,p]=initialGuess(Np, 0.5);

 % The number of discrete points in the first region
global N;
% The number of discrete points in the second region
global I; 
N=round(Np/2);
I=Np-N;
% Parameters:
%height H to which thing is compressed
% global H; 
%pressure inside region
global P; 
P=p;

[r] = residuals2(X);
% normR=norm(r);
[normR]=residueNorm(r, pi/2);
normString=sprintf('The norm is %.4e', normR);
figure();set(gcf,'Visible', 'off'); 
bar(1:N,r(1:N, 1), 'r', 'EdgeColor', 'r');

h=figure();set(gcf,'Visible', 'off'); 
% bar(r);

bar(1:N,r(1:N, 1), 'r', 'EdgeColor', 'r');
hold on;
bar(N+1:N+I, r(N+1:N+I, 1), 'b', 'EdgeColor', 'b');
bar(N+I+1:N+2*I, r(N+I+1:N+2*I, 1), 'g', 'EdgeColor', 'g');
bar(N+2*I+1, r(N+2*I+1, 1), 'k', 'EdgeColor', 'k');
legend('x', 'X', 'Y', '\zeta_c');
yminmax=ylim;

text(N/10,yminmax(1)/2,normString);
set(gca,'DefaultTextFontSize',16)
titleString=sprintf('Residue2 with Np = %d', Np);
title(titleString,'FontSize', 18)
% xlabel('Volumn Flux [ml/min]','FontSize', 16)
% ylabel('Distance from centreline [mm]','FontSize', 16)

filename=sprintf('residue2_Np=%d.png',Np);
if ispc
    print(h,strcat('residue\',filename), '-dpng')
elseif isunix
    print(h,strcat('residue/',filename), '-dpng')
else
    fprintf('\n Not saved, neither pc nor unix \n');
end
    