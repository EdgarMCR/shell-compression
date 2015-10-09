function[] = runScript_compressedSphere(Np, H, zetaca)
% fprintf('%.4f \n', H);
% Np = str2num(Np)
% p = str2num(p)
global N; global I; N=Np/2; I=Np-N;
close all;
% [X]=initialGuess_compressedSphere4(Np, H, zetac);
[X]=initialGuess_LM_Sphere(Np, zetaca);

%Plot initial guess
for ll=1:N;
    z(ll)=X(N+ll) * cos(X(ll));
    x(ll)=X(N+ll) * sin(X(ll));
end
for ll=1:I;
    z2(ll)=X(2*N + I + ll) * cos(X(2*N + ll));
    x2(ll)=X(2*N + I + ll) * sin(X(2*N + ll));
end
h1=figure();  plot(x, z, 'ro-'); hold on; plot(x2, z2, 'bs-');
titleString=sprintf('Inital Guess Compressed Sphere'); 
title(titleString,'FontSize', 12); xlabel('x','FontSize', 16); ylabel('z','FontSize', 16);
xlim([0 1.2]); ylim([0 1.2]);

global zetac;
zetac=zetaca;
global P;
P = 0.0;
fprintf('Fixing both zeta_c as well as p at %.2e, %.2e',  zetac, P);
% X(2*N+2*I+1,1)=0.00001;
% X(2*N+2*I+2,1)=zetaca;

path='sol.txt';
length(X);

% [r] = residuals_compressedSphere(X);
% area= H^2 * ( tan(X(round(Np/2),1)) + pi/4) ;
% fprintf(' residue = %.4e \t area = %.4e \n',  residueNorm(r, pi/2), area);

global lambda; %degree of flatness
lambda=0;

[r] = residuals_compressedSphere(X);
h1=figure(); x = linspace(1,2*N+2*I, 2*N+2*I);
plot(x(1:N), r(1:N), 'ob'); hold on; plot(x(N+1:2*N), r(N+1:2*N), 'or'); 
plot(x(2*N+1:2*N+I), r(2*N+1:2*N+I), 'sb'); plot(x(2*N+I+1:2*N+2*I), r(2*N+I+1:2*N+2*I), 'sr'); 
titleString=sprintf('Residue Inital Guess'); title(titleString,'FontSize', 12); xlabel('Entry in Solution Vector','FontSize', 16); ylabel('Residue','FontSize', 16);


% Y =problem_compressedSphere(Np,  H, X, path, 1, path);
fprintf('\n');
% 
% for ii =1:20;
%     P = P + 0.5;
%     fprintf('\n\nP = %.3e\n', P);
%     Y =problem_compressedSphere(Np,  H, Y, path, 1, path);
% end

% for ii =1:13;
%     lambda = lambda + 0.01;
%     fprintf('\n\nlambda = %.3f\n', lambda);
%     Y =problem_compressedSphere(Np,  H, Y, path, 1, path);
% end
% for ii =1:10;
%     lambda = lambda + 0.001;
%     fprintf('\n\nlambda = %.3f\n', lambda);
%     Y =problem_compressedSphere(Np,  H, Y, path, 1, path);
% end

% p=p+0.2;
% fprintf('\n\n run with p = %.2e \n',  p);
% problem_inflatingSphere(Np,  p, 0, path, 1, path);

% for ii=1:11
%     p=p+0.1;
%     fprintf('\n\n %d th run with p = %.2f \n', ii, p);
%     problem_inflatingSphere(Np,  p, 0, path, 1, path);
% end

% global L;
% 
% increment=0.01;
% leng=((10.0-p)/increment)+1;
% for ii=1:leng
%     close all;
%     L=ii;
%     p=p+increment;
%     fprintf('\n\n %d th run with p = %.2e \n', ii, p);
%     problem_inflatingSphere(Np,  p, 0, path, 1, path);
% end
fprintf('\n');
