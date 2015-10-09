function [r] = residuals_LM_sphere_adpatedToCompressed_zetac(X)
 % The number of discrete points in the first region
global N;
global I;

% fprintf('N=%d \n', N);

% Parameters:
%pressure inside region
global P;
% global zetac;

phi = X(1:N);
rs = X(N+1:2*N);

Phi = X(2*N+1:2*N+I);
R = X(2*N+I+1:2*N+2*I);

zetac = X(2*N+2*I+1);
% P = X(2*N+2*I+2);

delta1 = 1/(N-1);
delta2 =  1/(I-1);

dxidzeta1 = 1/zetac;
dxidzeta2 = 1/(pi/2-zetac);

r=zeros([2*N+2*I+1, 1]);

% rs=zeros([N, 1]);
% for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
% end


%##########################################################################
% Contact Region

% Boundary condition (first point in the first region): Symmetry dR/d\zeta = 0
% r(1, 1) = dxidzeta1*(-3*r(1) + 4*r(2) - r(3))/(2*delta1);
%Symmetry on y
r(1,1)=phi(1);
% r(N+1, 1) = dxidzeta1*(-3*rs(1) + 4*rs(2) - rs(3))/(2*delta1);
% r(N+1,1)=rs(1)-1;
r(N+1,1)=dxidzeta1 * (-3*rs(1) + 4*rs(2) - rs(3))/(2*delta1);

for ii = 2:N-1; %Using central difference
    xi = (ii-1)*delta1; %is actually xi
    zeta= xi*zetac;
    xiValue1(ii-1)=xi;
    zetaValue1(ii-1)=zeta;
    rc1=dxidzeta1*(rs(ii+1) -rs(ii-1))/(2*delta1);
    rc11=dxidzeta1^2 * (rs(ii+1) - 2*rs(ii) +rs(ii-1))/delta1^2;
    
    phic1=dxidzeta1*(phi(ii+1) -phi(ii-1))/(2*delta1);
    phic11= dxidzeta1^2 *(phi(ii+1) - 2*phi(ii) +phi(ii-1))/delta1^2;

    
    
    global lambda;
    global H;
    pContact =lambda * (rs(ii)*cos(phi(ii)) - H) + (1-lambda)* rs(ii)-1;
    pContactValue(ii-1)=pContact;
%     pContact=0.1;

    [reqn, phieqn , ~ , ~, ~, ~] = linearElasticity( rs(ii), rc1, rc11, phi(ii), phic1, phic11, zeta,  pContact);
	
    %equation 1 (for Phi)
    r(ii, 1) = phieqn;
    
    %equation 2 (for R)
    r(N+ii, 1) = reqn; 
end

%Boundary Conidtion - continuity
r(N,1) = phi(N)-Phi(1);

% r(2*N,1) = lambda * (rs(N)*cos(phi(N)) - H) + (1-lambda)*rs(N) - 1; %TODO What condition to impose here?
r(2*N,1) = dxidzeta1 * (3*rs(N) - 4*rs(N-1) + rs(N-2))/(2*delta1) - dxidzeta2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta2);




%##########################################################################
% Free Region

%boundary condition, dphi/dzeta = dPhi/dzeta
% r(2*N+1,1) = -dxidzeta1^2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta1) * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta1) * sin(Phi(1));
r(2*N+1,1) = dxidzeta1 * (3*phi(N) - 4*phi(N-1) + phi(N-2))/(2*delta1) - dxidzeta2 * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta2);
r(2*N+I+1,1)=rs(1)-R(1);

for jj = 2:I-1; %Using central difference
%     fprintf('jj = %d \t ', jj);
    xi = (jj-1)*delta2; %is actually xi
    zeta=zetac + xi*(pi/2-zetac);
    
    xiValue2(jj-1)=xi;
    zetaValue2(jj-1)=zeta;
    
    Rc1=dxidzeta2*(R(jj+1) -R(jj-1))/(2*delta2);
    Rc11=dxidzeta2^2 * (R(jj+1) - 2*R(jj) +R(jj-1))/delta2^2;
    
    Phic1=dxidzeta2*(Phi(jj+1) -Phi(jj-1))/(2*delta2);
    Phic11= dxidzeta2^2 *(Phi(jj+1) - 2*Phi(jj) +Phi(jj-1))/delta2^2;
    
    [reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity( R(jj), Rc1, Rc11, Phi(jj), Phic1, Phic11, zeta, P);
    
    %equation 1 (for Phi)
    r(2*N+jj, 1) = phieqn;
    
    %equation 2 (for R)
    r(2*N+I+jj, 1) =  reqn;
end

%bc , Phi at zeta= pi/2 = pi/2
r(2*N+I,1)=Phi(N)-pi/2;
%boundary condition,  symmetry on dR/d\zeta = 0
r(2*N+2*I,1) = dxidzeta2*(3*R(N) - 4*R(N-1) + R(N-2))/(2*delta2);

%##########################################################################
%zeta_c

% ii=NaN; jj=NaN; zeta=NaN; %to check for uuncorrected calls
clear ii; clear jj; clear zeta; clear phic11; clear Phic11;

rc1=dxidzeta1*(3*rs(N) - 4*rs(N-1) + rs(N-2))/(2*delta1);
rc11=dxidzeta1^2*(-2*rs(N) + 5*rs(N-1) - 4*rs(N-2) + rs(N-3))/(delta1^2);

phic1=dxidzeta1*(3*phi(N) - 4*phi(N-1) + phi(N-2))/(2*delta1);
phic11= dxidzeta1^2*(-2*phi(N) + 5*phi(N-1) - 4*phi(N-2) + phi(N-3))/(delta1^2);

% global lambda; global H;
pContact =lambda * (rs(N)*cos(phi(N)) - H) + (1-lambda)* rs(N)-1;

[~, ~ , sigma11_L, sigma11c1_L, sigma22s_L, e11_L ] = linearElasticity( rs(N), rc1, rc11, phi(N), phic1, phic11, zetac, pContact );

% fprintf('sigma11_L = %f , sigma11c1_L = %f , sigma22s_L = %f , e11_L = %f \n', sigma11_L, sigma11c1_L, sigma22s_L, e11_L);



Rc1=dxidzeta2*(-3*R(1) + 4*R(2) - R(3))/(2*delta2);
Rc11=dxidzeta2^2*(2*R(1) - 5*R(2) + 4*R(3) - R(4))/(delta2^2);

Phic1=dxidzeta2*(-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta2); %(Phi(jj+1) -Phi(jj-1))/(2*delta2);
Phic11=dxidzeta2^2*(2*Phi(1) - 5*Phi(2) + 4*Phi(3) - Phi(4))/(delta2^2);

[~, ~ , sigma11_R, sigma11c1_R, sigma22s_R, e11_R ] = linearElasticity( R(1), Rc1, Rc11, Phi(1), Phic1, Phic11, zetac, P );
% fprintf('sigma11_R = %f , sigma11c1_R = %f , sigma22s_R = %f , e11_R = %f \n', sigma11_R, sigma11c1_R, sigma22s_R, e11_R);

r(2*N+2*I+1,1)= (sigma11_R-sigma11_L)  + (sigma22s_R-sigma22s_L) ; %^ + (sigma11c1_R-sigma11c1_L); %TODO
% r(N+2*I+1,1)= sigma111-sigma11N; %TODO


% %##########################################################################
% %p
% %Ultimatly the volume constraint 
[v] = volumn_of_region(rs, phi, R, Phi, I, N);
% r(2*N+2*I+2,1)= v - pi/4; %TODO
% % r(N+2*I+2,1)= H^2 * ( tan(X(N,1)) + pi/4 ) - pi/4; %TODO


global plotLM;
global step;
global startTime;


if plotLM
    
%     figure()
%     n=2:N-1;
%     i=N+1:N+I-2;
% 
%     plot(n, xiValue1, 'ro-');
%     hold on;
%     plot(n, zetaValue1, 'bs-');
% 
%     plot(i, xiValue2, 'ro-');
%     plot(i, zetaValue2, 'bs-');
%     legend('\xi 1', '\zeta 1', '\xi 2' , '\zeta 2','Location','northwest');

    j1=figure();set(gcf,'Visible', 'off');
    n=2:N-1;
    i=N+1:N+I-2;

    plot(n, pContactValue, 'ro-');
    hold on;
    legend('pContact','Location','northwest');
    filename=sprintf('zetac\\%s_%d_step=%d_lambda=%.2e_Np=%d_p=%.2e.png',mfilename, startTime, step, lambda,   N, P);
    print(j1,filename, '-dpng')
end


end

