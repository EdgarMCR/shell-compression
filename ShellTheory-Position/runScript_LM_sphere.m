function[] = runScript_LM_sphere(Np, p, zetacz)
% fprintf('%.4f \n', H);
% Np = str2num(Np)
% p = str2num(p)
global H;
H=0.95;

global lambda;
lambda=0;


close all;
% [X]=initialGuess_compressedSphere4(Np, H, zetac);
[X]=initialGuess_LM_Sphere(Np, zetacz);
path='sol_LM_sphere.txt';
% length(X)

zetac=zetacz;
Y=problem_LM_sphere(Np,  p, X, path, 1, path);
fprintf('\n');



% for ii=1:2
%     if ii <3;
%         p=p+1;
%     else
%         p=p+1;
%     end
%     fprintf('%d th run with p=%f \n', ii, p);
%     Y=problem_LM_sphere(Np,  p, Y, path, 1, path);
%     fprintf('\n');
% end

% lambdaStep=1;
% for ii=1:20
%     lambda=lambda+lambdaStep;
%     fprintf('%d th run with lambda=%f \n', ii, lambda);
%     Y=problem_LM_sphere(Np,  p, Y, path, 1, path);
%     fprintf('\n');
% end

% diff=Y1-Y2;
% figure();
% plot(diff)
% p=p+0.2;
% fprintf('\n\n run with p = %.2e \n',  0.001);
% problem_LM_sphere(Np,  0.001, 0, path, 1, path);

% pStep=0.05;
% counter=0;
% while (p<1)
%     try
%         counter=counter+1;
%         p=p+pStep;
%         fprintf('\n\n %d th run with p = %.2f \n', counter, p);
%         problem_LM_sphere(Np,  p, 0, path, 1, path);
%     catch
%         p=p-pStep;
%         pStep=pStep/2;
%         fprintf('\n pStep halfed to %.2e', pStep);
%         if pStep < 1e-8
%             error('pStep to small');
%         end
%     end
% end


% target=1;
% increment=0.1;
% while (abs(target-flatnessFactor)>=0.1)
%     close all;
%     try
%         flatnessFactor=flatnessFactor+increment;
%         fprintf('\n\n %d th run with flatnessFactor = %.2e \n', ii, flatnessFactor);
%         problem_LM_sphere(Np,  p, 0, path, 1, path);
%     catch
%         flatnessFactor=flatnessFactor-increment;
%         increment=increment/2;
%         fprintf('increment halfed to %.2e \n', increment);
%         if increment < 1e-8
%             error('Increment to small!');
%         end
%     end
% end
fprintf('\n');
