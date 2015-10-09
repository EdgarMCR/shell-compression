function [] = problem_compressedSphere(Np, Hh, X, pathToSol,save_data,  savePath)
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

global H;
H=Hh;

% Parameters:
global E; %Youngs Modulus
global nu; %Poission Ratio

global h; %membran thickness

E=1000;
nu=0.5;
h=0.05;



if(xGiven==0)
    % Initial guess, should be reading from a file
    Data = importdata(pathToSol);
    X=zeros([2*N+2*I+2,1]);
    
    for ii = 1:length(Data);
        if(Data(ii, 4) == 0)
            Nh = ii;
            Ih = length(Data)-ii+1;
            if N ~= Nh || I ~= Ih
                fprintf('\n %s: Readin data will fail because of point mismatch N, Nh =  %d, %d \n \n', mfilename, N, Nh);
            end
            break;
        end
    end
    
    X(1:N, 1) = Data(1:length(Data), 2); 
    X(N+1:2*N, 1) = Data(1:length(Data), 3); 
%     X(2*N+1:2*N, 1) = Data(1:length(Data), 4); 
    X(N +1:N +I, 1) = Data(1:length(Data), 2); 
    X(N +I+1:N +2*I, 1) = Data(1:length(Data), 3); 
    X(N +2*I+1, 1) = Data(N, 0);
    X(N +2*I+2, 1) = 1;
end

func =@residuals_LM_sphere_adpatedToCompressed(X);

[Y]= finite_diff_Jacobian(X, func); 


if save_data == 1
    if(nargin < 6)
        savePath =   'calculated_values_InflatedSphere.txt';
    end

    filename = strcat(savePath); %num2str(N_p)

    % Open the file
    fvk = fopen(filename, 'w');
    
    zetac=Y(2*N+2*I+1,1);
    for ii=1:N;
        zeta = zetac * (ii-1)/(N-1);
        fprintf(fvk, '%4.12f\t %4.12f\t  %4.12f\t %4.12f \n', zeta, Y(ii, 1),  Y(N+ii, 1),Y(2*N+ii, 1));
    end
    for jj=1:I;
        zeta = zetac + (pi/2-zetac) * (jj-1)/(I-1);
        fprintf(fvk, '%4.12f\t %4.12f\t  %4.12f\t %4.12f \n', zeta, Y(N+jj, 1),  Y(N+I+jj, 1));
    end

    % Close the file
    fclose(fvk); 
end
for ll=1:N;
    z(ll)=Y(N+ll) * cos(Y(ll));
    x(ll)=Y(N+ll) * sin(Y(ll));
end

h1=figure(); %set(gcf,'Visible', 'off'); 
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
% filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%04d.png',L);
filename=sprintf('InflatedSphere\\Plot_Solution_InflatedSphere_Np=%d_p=%.2e.png',Np, p);
print(h1,filename, '-dpng')

e = cputime-t;
fprintf('CPU Runtime = %.2f seconds ', e)
end


   
