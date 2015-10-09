function [r] = residuals_compressedSphere(X)
 % The number of discrete points in the first region
global N;
global I;

% fprintf('N=%d \n', N);

% Parameters:
%pressure inside region
global H; 


global E; %Youngs Modulus
global nu; %Poission ratio

global h; %membran thickness
% the distance between the discrete points first region

phi = X(1:N);
rs = X(N+1:2*N);
% lm = X(2*N+1:3*N);

Phi = X(2*N+1:2*N+I);

R = X(2*N+I+1:2*N+2*I);

%fixing P and zeta_c
% P = X(2*N+2*I+1);
% zetac = X(2*N+2*I+2);
global zetac;
global P;

delta1 = 1/(N-1);
delta2 =  1/(I-1);

dxidzeta1 = 1/zetac;
dxidzeta2 = 1/(pi/2-zetac);

%fixing P and zeta_c
% r=zeros([2*N+2*I+2, 1]);
r=zeros([2*N+2*I, 1]);

% rs=zeros([N, 1]);
% for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
% end

sigma11N=0;
sigma111=0;


global lambda; %the degree of flatness

%##########################################################################
% Contact Region

% Boundary condition (first point in the first region): Symmetry dR/d\zeta = 0
% r(1, 1) = dxidzeta1*(-3*r(1) + 4*r(2) - r(3))/(2*delta1);
%Symmetry on y
r(1,1)=phi(1);
% r(N+1, 1) = rs(1) - H;
r(N+1, 1) = dxidzeta1 * (-3*rs(1) + 4*rs(2) - rs(3))/(2*delta1);

function [pContact] = enforceContactThroughPressure(rad, phia, lambda, H)
%     pContact= (1-lambda)*(rad-1) + lambda* (rad*cos(phia) - H);
    pContact=P;
end

for ii = 2:N-1; %Using central difference
    xi = (ii-1)*delta1; %is actually xi
    zeta= xi*zetac;
%     fprintf('\nzeta and phi: %.3f, %.3f', zeta, phi(ii));
%     xiValue1(ii-1)=xi;
%     zetaValue1(ii-1)=zeta;
    rc1=dxidzeta1*(rs(ii+1) -rs(ii-1))/(2*delta1);
    rc11=dxidzeta1^2 * (rs(ii+1) - 2*rs(ii) +rs(ii-1))/delta1^2;
    
    phic1=dxidzeta1*(phi(ii+1) -phi(ii-1))/(2*delta1);
    phic11= dxidzeta1^2 *(phi(ii+1) - 2*phi(ii) +phi(ii-1))/delta1^2;
%     fprintf('rs(ii) = %.2e \t rc1 = %.2e \t phi(ii) = %.2e \t phic1 = %.2e \n ', rs(ii), rc1, phi(ii), phic1);
%     fprintf('phic1 = %.2e \t rc1 = %.2e \t ', phic1, rc1);
    [pContact] = enforceContactThroughPressure(rs(ii), phi(ii), lambda, H);
    
	[reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( rs(ii), rc1, rc11, phi(ii), phic1, phic11, zeta, pContact );
    %equation 1 (for Phi)
    %  h  \left[\left[ \sigma^{11} \left( \Phi - \zeta - R_{,1}  \right)  + \sigma^{22} \left( \left( (\Phi - \zeta) \frac{ \cos^2(\zeta) }{ \sin^2(\zeta)} + R \frac{ \cos(\zeta) }{ \sin(\zeta)} \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(\Phi_{,1} + R - 1 \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p \ (\Phi- \zeta)  \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(ii, 1) = phieqn;
    
    %equation 2 (for R)
    % h  \left[\left[ \sigma^{11} \left( \Phi_{,1} + R - 1  \right)  +   \sigma^{22} \left( \left( R + (\Phi - \zeta) \frac{ \cos (\zeta) }{ \sin (\zeta)}   \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(R_{,1} + \zeta - \Phi  \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p  R \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(N+ii, 1) = reqn;
end

%Last point, do backward differencing
ii=N; zeta=zetac;
rc1=dxidzeta1*(3/2 * rs(ii) -2* rs(ii-1) + 1/2 * rs(ii-2))/(delta1);
rc11=dxidzeta1^2 * (2 * rs(ii) -5* rs(ii-1) + 4 * rs(ii-2) - rs(ii-3))/delta1^2;
    
phic1=dxidzeta1*(3/2 * phi(ii) -2* phi(ii-1) + 1/2 * phi(ii-2))/(delta1);
phic11= dxidzeta1^2 *(2 * phi(ii) -5* phi(ii-1) + 4 * phi(ii-2) - phi(ii-3))/delta1^2;

[pContact] = enforceContactThroughPressure(rs(ii), phi(ii), lambda, H);

[reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( rs(ii), rc1, rc11, phi(ii), phic1, phic11, zeta, pContact );
%equation 1 (for Phi)
% r(ii, 1) = phieqn; %Overwritten by boundary condition

%equation 2 (for R)
% r(N+ii, 1) = reqn;
    
%Boundary Conidtion - continuity
% r(N,1) = phi(N)-Phi(1);

% r(2*N,1) = rs(N)-R(1);
% r(2*N,1) = rs(N)*cos(phi(N)) - H; %TODO What condition to impose here?
% r(2*N,1) = (1-lambda)*(rs(N)-1) + lambda* (rs(N)*cos(phi(N)) - H);
r(N,1) = dxidzeta1 * (3*phi(N) + 4*phi(N-1) - phi(N-2))/(2*delta1) - dxidzeta2 * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta2);
r(2*N,1) = dxidzeta1 * (3*rs(N) + 4*rs(N-1) - rs(N-2))/(2*delta1) - dxidzeta2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta2);

%##########################################################################
% Free Region

%boundary condition, dZ/dzeta = 0
% r(2*N+1,1) = -dxidzeta2^2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta2) * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta2) * sin(Phi(1));
r(2*N+1,1) = phi(N) - Phi(1);
r(2*N+I+1,1)=rs(N)-R(1);


%first point with forward differencing
jj=1; zeta=zetac;

Rc1=dxidzeta2*(-3/2 * R(jj) +2* R(jj+1) - 1/2 * R(jj+2))/(delta2);
Rc11=dxidzeta2^2 * (2 * R(jj) -5* R(jj+1) + 4 * R(jj+2) - R(jj+3))/delta2^2;

Phic1=dxidzeta2*(-3/2 * Phi(jj) +2* Phi(jj+1) - 1/2 * Phi(jj+2))/(delta2);
Phic11= dxidzeta2^2 *(2 * Phi(jj) -5* Phi(jj+1) + 4 * Phi(jj+2) - Phi(jj+3))/delta2^2;

[reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( R(ii), Rc1, Rc11, Phi(ii), Phic1, Phic11, zeta, P );
%equation 1 (for Phi)
% r(2*N+jj, 1) = phieqn;%Boundary condition

%equation 2 (for R)
% r(2*N+I+jj, 1) = reqn; %Boundary condition

for jj = 2:I-1; %Using central difference
%     fprintf('jj = %d \t ', jj);
    xi = (jj-1)*delta2; %is actually xi
    zeta=zetac + xi*(pi/2-zetac);
%     fprintf('\nzeta and Phi: %.3f, %.3f', zeta, Phi(jj));
    Rc1=dxidzeta2*(R(jj+1) -R(jj-1))/(2*delta2);
    Rc11=dxidzeta2^2 * (R(jj+1) - 2*R(jj) +R(jj-1))/delta2^2;
    
    Phic1=dxidzeta2*(Phi(jj+1) -Phi(jj-1))/(2*delta2);
    Phic11= dxidzeta2^2 *(Phi(jj+1) - 2*Phi(jj) +Phi(jj-1))/delta2^2;
%     fprintf('R(jj) = %.2e \t Rc1 = %.2e \t Phi(jj) = %.2e \t Phic1 = %.2e \n ', R(jj), Rc1, Phi(jj), Phic1);
%     fprintf('Phic1 = %.2e \t Rc1 = %.2e \t ', Phic1, Rc1);

    [reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( R(ii), Rc1, Rc11, Phi(ii), Phic1, Phic11, zeta, P );
    %equation 1 (for Phi)
    %  h  \left[\left[ \sigma^{11} \left( \Phi - \zeta - R_{,1}  \right)  + \sigma^{22} \left( \left( (\Phi - \zeta) \frac{ \cos^2(\zeta) }{ \sin^2(\zeta)} + R \frac{ \cos(\zeta) }{ \sin(\zeta)} \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(\Phi_{,1} + R - 1 \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p \ (\Phi- \zeta)  \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(2*N+jj, 1) = phieqn;
    
    %equation 2 (for R)
    % h  \left[\left[ \sigma^{11} \left( \Phi_{,1} + R - 1  \right)  +   \sigma^{22} \left( \left( R + (\Phi - \zeta) \frac{ \cos (\zeta) }{ \sin (\zeta)}   \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(R_{,1} + \zeta - \Phi  \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p  R \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(2*N+I+jj, 1) = reqn;
end

%bc , Phi ( at zeta= pi/2) = pi/2
r(2*N+I,1)=Phi(N)-pi/2;
%boundary condition,  symmetry on dR/d\zeta = 0
r(2*N+2*I,1) = dxidzeta2*(3*R(N) - 4*R(N-1) + R(N-2))/(2*delta2);







%##########################################################################
%p
%Ultimatly the volume constraint 
% [v] = volumn_of_region(Phi, R);
% r(2*N+2*I+1,1)= v - pi/4; %TODO
% r(N+2*I+2,1)= H^2 * ( tan(X(N,1)) + pi/4 ) - pi/4; %TODO




%##########################################################################
%zeta_c

% rc1=dxidzeta1*(3*rs(N) - 4*rs(N-1) + rs(N-2))/(2*delta1);
% % rc11=dxidzeta1^2 * (rs(ii+1) - 2*rs(ii) +rs(ii-1))/delta1^2;
% 
% phic1=dxidzeta1*(3*phi(N) - 4*phi(N-1) + phi(N-2))/(2*delta1);
% % phic11= dxidzeta1^2 *(phi(ii+1) - 2*phi(ii) +phi(ii-1))/delta1^2;
% 
% 
% %e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
% e11= phic1*rs(N) - phi(N)*rc1 - zeta*rc1 - zeta * phi(N) - phic1 - rs(N) +...
%     1/2*( (phi(N))^2 + zeta^2 + (rs(N))^2 + phic1^2 + rc1^2 );
% 
% % e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
% e22s = ( phi(N) - zeta ) * rs(N) * ( cos(zeta)/sin(zeta) ) + ...
%     1/2 * ( phi(N) - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (rs(N))^2 - 1 );
% 
% %\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
% sigma11_LHS = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
% 
% 
% Rc1=dxidzeta2*(-3*R(1) + 4*R(2) - R(3))/(2*delta2);
% 
% Phic1=dxidzeta2*(-3*Phi(N) + 4*Phi(N-1) - Phi(N-2))/(2*delta2); %(Phi(jj+1) -Phi(jj-1))/(2*delta2);
% 
% %e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
% e11= Phic1*R(1) - Phi(1)*Rc1 - zetac*Rc1 - zeta * Phi(1) - Phic1 - R(1) +...
%     1/2*( (Phi(1))^2 + zetac^2 + (R(1))^2 + Phic1^2 + Rc1^2 );
% 
% % e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
% e22s = ( Phi(1) - zetac ) * R(jj) * ( cos(zetac)/sin(zetac) ) + ...
%     1/2 * ( Phi(1) - zetac)^2 * ( cos(zetac)/sin(zetac) )^2 + 1/2 * ( (R(1))^2 - 1 );
% 
% %\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
% sigma11_RHS = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
%sigma^{11} = sigma^{11} at zeta_c
% r(2*N+2*I+2,1)= (sigma11_RHS-sigma11_LHS) ; %TODO
% r(N+2*I+1,1)= sigma111-sigma11N; %TODO

function [v] = volumn_of_region(Phi, R)
    x_critical=R(1)*sin(Phi(1));
    z_critical=R(1)*cos(Phi(1));
%         if( le((abs(z_critical-H)),  10^(-8)))
%             fprintf('\n %.1e \t diff = %.4e  large \n', 1e-8, abs(z_critical-H));
%         end
    v=x_critical*z_critical + pi/4 * z_critical^2;
end

% n=1:N-2;
% size(pContactValue)
% 
% % figure(5)
% % plot(n, pContactValue, 'ro');
% % title('pContactValue')

end

