function [Y] = problem_inflatingSphere_nonDim_pVariable(Np, p, X, pathToSol,save_data,  savePath)
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
N=Np;
% global I;
% N=round(Np/2); I=Np-I;


% Parameters:
global P;
% global lambda; %first lame parameter
% global mu; %second lame parameter
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1000;
nu=0.5;
h=1;

P=p;




if(xGiven==0)
    % Initial guess, should be reading from a file
    Data = importdata(pathToSol);
    X=zeros([2*N+1,1]);
    X(1:N, 1) = Data(1:length(Data)-1, 2); 
    X(N+1:2*N, 1) = Data(1:length(Data)-1, 3); 
    X(2*N+1,1) = p; 
end

func =@residuals_inflatingSphere_nonDim_pVariable;

[Y]= finite_diff_Jacobian(X, func); 


if save_data == 1
    if(nargin < 6)
        savePath =   'calculated_values_InflatedSphere.txt';
    end

    filename = strcat(savePath); %num2str(N_p)

    % Open the file
    fvk = fopen(filename, 'w');

    for ii=1:N;
        zeta = pi/2 * (1-(ii-1)/(N-1));;
        fprintf(fvk, '%4.12f\t %4.12f\t %4.12f \n', zeta, Y(ii, 1),  Y(N+ii, 1));
    end
    fprintf(fvk, '%4.12f\n', Y(2*N+1, 1));
    % Close the file
    fclose(fvk); 
end
for ll=1:N;
    z(ll)=Y(N+ll) * sin(Y(ll));
    x(ll)=Y(N+ll) * cos(Y(ll));
end

global L;
h=figure(); %set(gcf,'Visible', 'off'); 
% plot(Y(1:N,1), Y(N+1:2*N,1), 'ro-')
plot(x, z, 'ro-')
titleString=sprintf('Solution Inflated Sphere Np=%d p=%.2e', Np, p);
title(titleString,'FontSize', 16)
xlabel('x','FontSize', 16)
ylabel('z','FontSize', 16)
xlimR=xlim;

xmax=xlimR(2);
xlim([0 xmax])
ylim([0 xmax])

% str1 = 'Using pressure profile P =-p*cos(2*\zeta)';
% text(xmax/10,xmax/4,str1)

global counter;
% filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_cosineP_Np=%d_%d.png',Np, counter);

% filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%04d.png',L);
filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%d_p=%.2e.png',Np, p)
% filename=sprintf('InflatedSphere/Plot_Solution_InflatedSphere_Np=%d_p=%.2e.png',Np, p);
print(h,filename, '-dpng')

e = cputime-t;
fprintf('CPU Runtime = %.2f seconds ', e)
end


   
