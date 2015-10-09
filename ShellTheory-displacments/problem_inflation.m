function [Y] = problem_inflation(Np, Pp, X, pathToSol,save_data,  savePath)
% fprintf('Started %s \t', mfilename);
if length(X)<Np 
    xGiven=0;
else
    xGiven=1;
end

if(nargin < 5)
    save_data=1;
end
if(nargin < 6)
    savePath =   'calculated_values_test.txt';
end


t = cputime;
% The number of discrete points in the first region
global N;
N=round(Np);

global P; P=Pp;

% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio
global h; %membran thickness

E=1;
nu=0.5;
h=0.05;



if(xGiven==0)
    % Initial guess, should be reading from a file
    Data = importdata(pathToSol);
    X=zeros([2*N,1]);
    X(1:N, 1) = Data(1:length(Data), 2); 
    X(N+1:2*N, 1) = Data(1:length(Data), 3); 
end


% func =@residuals_LM_sphere;
fprintf('\n Using %s as residue. \n', 'residuals_inflation');
func =@residuals_inflation;
fprintf('\n Before Newton solve\n');
[Y]= finite_diff_Jacobian(X, func);

if save_data == 1
    if(nargin < 6)
        savePath =   'calculated_values_InflatedSphere.txt';
    end

    filename = strcat(savePath); %num2str(N_p)

    % Open the file
    fvk = fopen(filename, 'w');
    
    for ii=1:N;
        zeta = pi/2 * (ii-1)/(N-1);
        fprintf(fvk, '%4.12f\t %4.12f\t  %4.12f\t %4.12f  \n', zeta, Y(ii, 1),  Y(N+ii, 1), 1);
    end
    % Close the file
    fclose(fvk); 
end
%Plot the solution, here the solution vector is in displacements
delta1 = 1/(N-1);
for ll=1:N;
    xi = (ll-1)*delta1; %is actually xi
    zeta= xi * (pi/2);
    z1(ll)=(1.0 + Y(N+ll)) * cos((zeta + Y(ll)));
    x1(ll)=(1.0 + Y(N+ll)) * sin((zeta + Y(ll)));
end

% p=Y(2*N+2*I+2);
p=P;
h1=figure(); %set(gcf,'Visible', 'off'); 
% plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
plot(x1, z1, 'ro-')
hold on;
titleString=sprintf('Solution Inflated Sphere w/ LM Np=%d p=%.2e', Np, p);
title(titleString,'FontSize', 16)
xlabel('x','FontSize', 16)
ylabel('z','FontSize', 16)
xlimR=xlim;
xmax=xlimR(2)*1.2;
xlim([0 1.49])
ylim([0 1.49])
% filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%04d.png',L);
global lambda;
global startTime;
filename=sprintf('InflatedSphere\\Plot_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
print(h1,filename, '-dpng')

% h1=figure(); %set(gcf,'Visible', 'off'); 
% % plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
% plot(Y(1:N), 'ro')
% hold on;
% plot(Y(N+1:2*N), 'bs')
% titleString=sprintf('Solution displacement Vector Np=%d p=%.2e', Np, p);
% title(titleString,'FontSize', 12)
% xlabel('u1 / u3','FontSize', 16)
% ylabel('value','FontSize', 16)
% ylim([0 max(Y)])



e = cputime-t;
fprintf('CPU Runtime = %.2f seconds ', e)
end


   
