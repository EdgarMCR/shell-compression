function[] = runScript_contactSphere(Np, H, zetac)
% fprintf('%.4f \n', H);
% Np = str2num(Np)
% p = str2num(p)

close all;
[X]=initialGuess_contactRegion(Np, H, zetac);
path='sol.txt';
problem_contactRegion(Np,  H, X, path, 1, path);
fprintf('\n');

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
