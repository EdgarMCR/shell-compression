function[] = runScript_InflatedSphere_nonDim(Np, p)
% fprintf('%.4f \n', H);
% Np = str2num(Np)
% p = str2num(p)
global counter;

counter=1;
close all;
[X]=initialGuess_inflatedSphere(Np);
path='sol.txt';
Y=problem_inflatingSphere_nonDim(Np,  p, X, path, 1, path);
fprintf('\n');


[r1] = residuals_inflatingSphere_nonDim(X);
[r] = residuals_inflatingSphere_nonDim(Y);
fprintf(' residue X = %.4e \t residue Y = %.4e  \n',  residueNorm(r1, pi/2),  residueNorm(r, pi/2));
fprintf('\n');
for ii=1:20;
    p = p + 50;
    fprintf('\n \tii=%d', ii);
    Y=problem_inflatingSphere_nonDim(Np,  p, Y, path, 1, path);
end

% counter=counter+1;
% pStep=1.0;
% while p<100
%     try
%         p=p+pStep;
%         fprintf('p = %.2e \n', p);
%         Y=problem_inflatingSphere_nonDim(Np,  p, 0, path, 1, path);
%         counter=counter+1;
%     catch
%         p=p-pStep;
%         pStep=pStep/2.0;
%         fprintf('\n reducing pStep to %.2e \n', pStep);
%         if pStep < 1e-8
%             error('pStep to small! Less than 1e-8. ');
%         end
%     end
%     close all;
% end
    




% for kk=1:2*Np
%     diff(kk)=X(kk)- Y(kk);
% end
% 
% for ii = 1:Np
%     xpos(ii)= X(Np+ii) * cos(X(ii));
%     zpos(ii)= X(Np+ii) * sin(X(ii));
%     xpos2(ii)= Y(Np+ii) * cos(Y(ii));
%     zpos2(ii)= Y(Np+ii) * sin(Y(ii));
%     xdiff(ii)= diff(Np+ii) * cos(diff(ii));
%     zdiff(ii)= diff(Np+ii) * sin(diff(ii));
% end

% plot(xpos, zpos, 'ro', 'MarkerFaceColor', 'r')
% hold on;
% plot(xdiff, zdiff, 'bs', 'MarkerFaceColor', 'b')
% plot(xpos, zpos, 'kx', 'MarkerFaceColor', 'k')

% plot(diff)



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
