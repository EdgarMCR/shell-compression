function [Xnew] = finite_diff_Jacobian(Xnew, func)
% fprintf('Started Newton Solve \t');
% the number of degrees of freedom
n = length(Xnew);
% fprintf(' n = %d \n', n);



%for plotting LM


global step;
global plotLM;
global startTime;
plotLM=false;


% the count for the number of iterations
step = 0;

% increments by which we change variables 
eps = 1e-8; 
% eps = 1e-3; 

% error on the size of the residuals for finite-differencing 
old_norm = 1e-7;

% the number of maximum iterations
r = 10;

Jac=zeros([n, n]);

plotStuff=false;

%track zetac
global lambda;
counter=1;
zetacValue(counter)=Xnew(n);

while 1
    % Solving the equation:
    % need to finite-difference the equation in order to obtain the jacobian
    [Residual] = feval(func, Xnew);
    if step==0
        initialNorm=residueNorm(Residual, pi/2);
        residueValue(counter)=initialNorm;
        fprintf('\n Initial Residue Err = %.6e \t', initialNorm);
%         if ((initialNorm < old_norm))
%             fprintf('Newton Steps=%d \t',step);
%             break;      
%         end
    end
%     fprintf('zeta_c residue = %.4e \t ', Residual(n,1));
    for ii=1:n;
        Xn = Xnew;
        Xn(ii) = Xnew(ii)+eps;
        [Residual_eps{ii}] = feval(func, Xn);
    end

    % Building the Jacobian
%     fprintf('\n length(Jac(:, 1)) =%d, length(Residual) =%d, length(Residual_eps{jj}) =%d,length(eps) =%d\n', length(Jac(:, 1)), length(Residual),length(Residual_eps{1}),length(eps));
    for jj=1:n;
        Jac(:, jj) = (Residual-Residual_eps{jj})/eps;
    end
    
%     dlmwrite('Jac.txt',Jac, ' ')
%     pause
    
    Xnew = Xnew + Jac\Residual;
    
    % the iteration number
	step = step+1;
    
    plotLM=true;
    [Residual] = feval(func, Xnew);
    plotLM=false;
    
%     if plotStuff
%         global N; global I;
%         for ll=1:N;
%             z1(ll)=Xnew(N+ll) * cos(Xnew(ll));
%             x1(ll)=Xnew(N+ll) * sin(Xnew(ll));
%         end
%         for kk=1:I;
%             z2(kk)=Xnew(2*N+I+kk) * cos(Xnew(2*N+kk));
%             x2(kk)=Xnew(2*N+I+kk) * sin(Xnew(2*N+kk));
%         end
%         h1=figure(); %set(gcf,'Visible', 'off'); 
%         % plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
%         plot(x1, z1, 'ro-')
%         hold on;
%         plot(x2, z2, 'bs-')
%         titleString=sprintf('Solution Sphere w/ LM %d lambda=%f NewtonStep=%d',startTime, lambda, step);
%         title(titleString,'FontSize', 16)
%         xlabel('x','FontSize', 16)
%         ylabel('z','FontSize', 16)
%         xlimR=xlim;
%         xmax=xlimR(2)*1.2;
%         xlim([0 xmax])
%         ylim([0 xmax])
%         filename=sprintf('zetac\\Plot_%d_lambda=%.2e_NewtonStep=%d.png',startTime,lambda, step);
%         print(h1,filename, '-dpng')
%     end 

    counter=counter+1;
    zetacValue(counter)=Xnew(n);
    
    % Displaying pressure and lambda1 at 0 
    fprintf('\n Newton Method Step =  %1d Err = %.6e \t',step, residueNorm(Residual, pi/2));
    residueValue(counter)=residueNorm(Residual, pi/2);

    
    % If the count for the number of iterations exceeds 100, the intial guess is really - really bad
    % -> need to stop and provide a better guess
    if (step == r)
        fprintf('\n Second to last %e and last entry in residue = %e \n',Residual(n-1),Residual(n));
        if plotStuff
            h1=figure(); %set(gcf,'Visible', 'off'); 
            % plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
            xplot=0:counter-1;
            [AX,H1,H2] = plotyy(xplot,zetacValue,xplot,residueValue,'plot');
            set(get(AX(1),'Ylabel'),'String','\zeta_c')
            set(get(AX(2),'Ylabel'),'String','Residue')
            xlabel('Time (\musec)')
            titleString=sprintf('zeta_c  and residue through Newton Iterations');
            set(H1,'LineStyle','--', 'Marker','o','MarkerFaceColor','blue')
            set(H2,'LineStyle',':', 'Marker','s','MarkerFaceColor','green')
            title(titleString,'FontSize', 14)
            xlabel('Newton Step','FontSize', 16)
            ylabel('\zeta_c','FontSize', 16)
            filename=sprintf('zetac\\Plot_zeta_c_%d_lambda_=%e_%d.png', startTime,lambda, randi(50));
            print(h1,filename, '-dpng')
        end
        error('Jac:NotConverging','Reached Newton Steps limit, problem is not converging');
    end
    
    clear('Jac');
    clear('Residual_eps');
    
    % The iterations stop when the accuracy is satisfactory or the method
    % has had too many iterations
    if ((residueNorm(Residual, pi/2) < old_norm))
        fprintf('Newton Steps=%d \t',step);
        break;      
    end
    
end

if plotStuff
    h1=figure(); %set(gcf,'Visible', 'off'); 
    % plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
    xplot=0:counter-1;
    [AX,H1,H2] = plotyy(xplot,zetacValue,xplot,residueValue,'plot');
    set(get(AX(1),'Ylabel'),'String','\zeta_c')
    set(get(AX(2),'Ylabel'),'String','Residue')
    xlabel('Time (\musec)')
    titleString=sprintf('zeta_c  and residue through Newton Iterations');
    set(H1,'LineStyle','--', 'Marker','o','MarkerFaceColor','blue')
    set(H2,'LineStyle',':', 'Marker','s','MarkerFaceColor','green')
    title(titleString,'FontSize', 14)
    xlabel('Newton Step','FontSize', 16)
    ylabel('\zeta_c','FontSize', 16)
    filename=sprintf('zetac\\Plot_zeta_c_%d_lambda_=%e_%d.png', startTime, lambda, randi(50));
    print(h1,filename, '-dpng')
end

end


