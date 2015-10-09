function [Y] = problem_LM_sphere_zetac(Np, Pp, X, pathToSol,save_data,  savePath)
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
global I;
N=round(Np/2); I=Np-N;

global P;
P=Pp;
global zetac;


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
%     X=zeros([2*N+2*I+1,1]);
    X=zeros([2*N+2*I,1]);
    for ii = 1:length(Data);
        if(Data(ii, 4) == 0)
            Nh = ii;
            Ih = length(Data)-ii+1;
            if N ~= Nh || I ~= Ih
                fprintf('\n %s: Readin data will fail because of point mismatch N, Nh =  %d, %d or I, Ih =  %d, %d \n \n', mfilename, N, Nh, I, Ih);
            end
            break;
        end
    end
    
    X(1:N, 1) = Data(1:Nh, 2); 
    X(N+1:2*N, 1) = Data(1:Nh, 3); 
%     X(2*N+1:2*N, 1) = Data(1:length(Data), 4); 
    X(2*N + 1 :2 *N + I, 1) = Data(Nh:length(Data), 2); 
    X(2*N + I + 1 : 2*N + 2*I, 1) = Data(Nh:length(Data), 3); 
end
%  X(2*N+2*I+1,1)=zetac;
% X(2*N+2*I+2,1)=Pp;

% func =@residuals_LM_sphere;
fprintf('\n Using %s as residue. \n', 'residuals_LM_sphere_adpatedToCompressed_zetac');
func =@residuals_LM_sphere_adpatedToCompressed_zetac;
fprintf('\n Before Newton solve zetac is %f \n', zetac);
[Y]= finite_diff_Jacobian(X, func);
% zetac=Y(2*N+2*I+1,1);
fprintf('\n After Newton solve zetac is %f \n', zetac);

if save_data == 1
    if(nargin < 6)
        savePath =   'calculated_values_InflatedSphere.txt';
    end

    filename = strcat(savePath); %num2str(N_p)

    % Open the file
    fvk = fopen(filename, 'w');
    
    for ii=1:N-1;
        zeta = zetac * (ii-1)/(N-1);
        fprintf(fvk, '%4.12f\t %4.12f\t  %4.12f\t %4.12f  \n', zeta, Y(ii, 1),  Y(N+ii, 1), 1);
    end
    fprintf('\n I = %d \n', I);
    for jj=1:I;
        zeta = zetac + (pi/2-zetac) * (jj-1)/(I-1);
        fprintf(fvk, '%4.12f\t %4.12f\t  %4.12f\t %4.12f  \n', zeta, Y(2*N+jj, 1),  Y(2*N+I+jj, 1), 0);
    end

    % Close the file
    fclose(fvk); 
end
%Plot the solution, here the solution vector is in displacements

delta1 = 1/(N-1);
delta2 =  1/(I-1);
for ll=1:N;
    xi = (ll-1)*delta1; %is actually xi
    zeta= xi*zetac;
    z1(ll)=(1.0 + Y(N+ll)) * cos((zeta + Y(ll)));
    x1(ll)=(1.0 + Y(N+ll)) * sin((zeta + Y(ll)));
end
for kk=1:I;
    xi = (kk-1)*delta2; %is actually xi
    zeta=zetac + xi*(pi/2-zetac);
    z2(kk)=(1.0+Y(2*N+I+kk)) * cos(zeta +Y(2*N+kk));
    x2(kk)=(1.0+Y(2*N+I+kk)) * sin(zeta +Y(2*N+kk));
end

% p=Y(2*N+2*I+2);
p=P;
h1=figure(); %set(gcf,'Visible', 'off'); 
% plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
plot(x1, z1, 'ro-')
hold on;
plot(x2, z2, 'bs-')
titleString=sprintf('Solution Inflated Sphere w/ LM Np=%d p=%.2e', Np, p);
title(titleString,'FontSize', 16)
xlabel('x','FontSize', 16)
ylabel('z','FontSize', 16)
xlimR=xlim;
xmax=xlimR(2)*1.2;
xlim([0 1.5])
ylim([0 1.5])
% filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%04d.png',L);
global lambda;
global startTime;
filename=sprintf('NoPorSetac_Solution\\Plot_Solution_%d_lambda=%.2e_Np=%d_p=%.2e.png',startTime,lambda, Np, p);
print(h1,filename, '-dpng')

fprintf('\n u1(N),  U1(1) = %.2e, %.2e  and u3(N),  U3(1) = %.2e, %.2e \n \n', Y(N), Y(2*N+1), Y(2*N), Y(2*N+I+1));
e = cputime-t;
fprintf('CPU Runtime = %.2f seconds ', e)
end


   
