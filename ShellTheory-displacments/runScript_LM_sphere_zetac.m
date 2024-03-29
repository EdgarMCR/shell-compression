function[] = runScript_LM_sphere_zetac(Np, p, zetacz)
% fprintf('%.4f \n', H);
% Np = str2num(Np)
% p = str2num(p)
global H;
H=0.0;

global lambda;
lambda=0;

t0=clock;
global startTime;
startTime=t0(6)*1000;


close all;
% [X]=initialGuess_compressedSphere4(Np, H, zetac);
[X]=initialGuess_LM_Sphere(Np, zetacz);
path='LM_Sphere_zetac\\sol_LM_sphere_zetac_0.txt';
% length(X)

% [r] = residuals_compressedSphere(X);
% area= H^2 * ( tan(X(round(Np/2),1)) + pi/4) ;
% fprintf(' residue = %.4e \t area = %.4e \n',  residueNorm(r, pi/2), area);
global zetac;
zetac=zetacz;
Y=problem_LM_sphere_zetac(Np,  p, X, path, 1, path);
fprintf('\n');
count=0;

for ii=1:10;
    if ii < 5
        p = p +0.01;
    else
        p = p +0.001;
    end
    fprintf('\n \tii = %d', ii)
    Y=problem_LM_sphere_zetac(Np,  p, Y, path, 1, path);
end
% 
% N=round(Np/2); I=Np-N;
% for ll=1:N;
%     z1(1,ll)=Y(N+ll) * cos(Y(ll));
%     x1(1,ll)=Y(N+ll) * sin(Y(ll));
%     r1(1,ll)=Y(N+ll);
%     psi1(1,ll)=Y(ll);
% end
% for kk=1:I;
%     z2(1,kk)=Y(2*N+I+kk) * cos(Y(2*N+kk));
%     x2(1,kk)=Y(2*N+I+kk) * sin(Y(2*N+kk));
%     r2(1,kk)=Y(2*N+I+kk);
%     psi2(1,kk)=Y(2*N+kk);
% end
% count=count+1;
% 
% 
% lambdaStep=1e-7;
% for ii=1:10
% %     for ii=1:length(Y)
% %         Y(ii)= Y(ii) + (-1)^(ii) * 1e-8;
% %     end
%     lambda=lambda+lambdaStep;
%     fprintf('zetac = %e \n %d th run with lambda=%f ',Y(length(Y)), ii, lambda);
%     path=sprintf('LM_Sphere_zetac\\sol_LM_sphere_zetac_%d.txt', ii);
%     try
%         Y=problem_LM_sphere_zetac(Np,  p, Y, path, 1, path);
%     catch
%         break;
%     end
%     
%     for ll=1:N;
%         z1(ii+1,ll)=Y(N+ll) * cos(Y(ll));
%         x1(ii+1,ll)=Y(N+ll) * sin(Y(ll));
%         r1(ii+1,ll)=Y(N+ll);
%         psi1(ii+1,ll)=Y(ll);
%     end
%     for kk=1:I;
%         z2(ii+1,kk)=Y(2*N+I+kk) * cos(Y(2*N+kk));
%         x2(ii+1,kk)=Y(2*N+I+kk) * sin(Y(2*N+kk));
%         r2(ii+1,kk)=Y(2*N+I+kk);
%         psi2(ii+1,kk)=Y(2*kk);
%     end
%     count=count+1;
%     lastLambda=lambda;
%     fprintf('\n');
% end
% 
% h1=figure();
% plot(x2(1, :), z2(1,:), 'ro', 'MarkerFaceColor', 'r')
% hold on;
% plot(x1(1, :), z1(1,:), 'cs', 'MarkerFaceColor', 'c')
% title('first solution') 
% filename=sprintf('LM_Sphere_zetac\\First_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
% print(h1,filename, '-dpng')
% 
% h1=figure();
% plot(x2(count, :), z2(count,:), 'ro', 'MarkerFaceColor', 'r')
% hold on;
% plot(x1(count, :), z1(count,:), 'cs', 'MarkerFaceColor', 'c')
% title('last solution') 
% filename=sprintf('LM_Sphere_zetac\\Last_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
% print(h1,filename, '-dpng')
% 
% h1=figure();
% diffx1(1:N)=x1(1, 1:N)-x1(count, 1:N)
% diffz1(1:N)=z1(1, 1:N)-z1(count, 1:N)
% diffx2(1:I)=x2(1, 1:I)-x2(count, 1:I)
% diffz2(1:I)=z2(1, 1:I)-z2(count, 1:I)
% 
% diffr1(1:N)=r1(1, 1:N)-r1(count, 1:N)
% diffr2(1:I)=r2(1, 1:I)-r2(count, 1:I)
% diffpsi1(1:N)=psi1(1, 1:N)-psi1(count, 1:N)
% diffpsi2(1:I)=psi2(1, 1:I)-psi2(count, 1:I)
% 
% fprintf('Last solution psi = %.2e Psi= %.2e zetac = %.2e \t', psi1(count, N), psi2(count, 1), Y(2*N+2*I+1));
% fprintf('Last solution r = %.2e R= %.2e \n', r1(count, N), r2(count, 1));
% 
% plot(diffx1, diffz1, 'cs', 'MarkerFaceColor', 'c')
% hold on;
% plot(diffx2, diffz2, 'ro', 'MarkerFaceColor', 'r')
% title('Difference between initial and final displacement')
% legend('Contact Region', 'Free Region')
% filename=sprintf('LM_Sphere_zetac\\Difference_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
% print(h1,filename, '-dpng')
% 
% h1=figure();
% n=1:N;
% i=N+1:N+I;
% plot(n, diffr1, 'cs', 'MarkerFaceColor', 'c')
% hold on;
% plot(i, diffr2, 'ro', 'MarkerFaceColor', 'r')
% title('Difference between initial and final radius')
% legend('Contact Region', 'Free Region')
% filename=sprintf('LM_Sphere_zetac\\Difference_r_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
% print(h1,filename, '-dpng')
% 
% 
% 
% h1=figure();
% n=1:N;
% i=N+1:N+I;
% plot(n, diffpsi1, 'cs', 'MarkerFaceColor', 'c')
% hold on;
% plot(i, diffpsi2, 'ro', 'MarkerFaceColor', 'r')
% title('Difference between initial and final angle \psi')
% legend('Contact Region', 'Free Region')
% filename=sprintf('LM_Sphere_zetac\\Difference_psi_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
% print(h1,filename, '-dpng')

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
