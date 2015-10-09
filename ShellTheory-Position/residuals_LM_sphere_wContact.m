function [r] = residuals_LM_sphere_wContact(X)
 % The number of discrete points in the first region
global N;
global I;

% fprintf('N=%d \n', N);

% Parameters:
%pressure inside region
global P;
global zetac;

global lambda;
global H;

global E; %Youngs Modulus
global nu; %Poission ratio

global h; %membran thickness
% the distance between the discrete points first region

phi = X(1:N);
rs = X(N+1:2*N);

Phi = X(2*N+1:2*N+I);
R = X(2*N+I+1:2*N+2*I);


% P = X(2*N+2*I+1);
% zetac = X(2*N+1*I+2);

delta1 = 1/(N-1);
delta2 =  1/(I-1);

dxidzeta1 = 1/zetac;
dxidzeta2 = 1/(pi/2-zetac);

r=zeros([2*N+2*I, 1]);

% rs=zeros([N, 1]);
% for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
% end

    function [pContact] = contactPressure(r, phi, lambda, H, P)
        if lambda <=1;
            pContact =  + lambda * (r*cos(phi) - H) + ( 1-lambda)*(r-1);
        else
            pContact = + lambda * (r*cos(phi) - H);
        end
%         pContact = P;
    end
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
    
%     [pContact] = contactPressure(rs(ii), phi(ii), lambda, H, P);
    pContact = (rs(ii)*cos(phi(ii)) - H);
    [reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( rs(ii), rc1, rc11, phi(ii), phic1, phic11, zeta, pContact );
    
    r(ii, 1) = phieqn;
    r(N+ii, 1) = reqn;
end

%Boundary Conidtion - continuity
r(N,1) = phi(N)-Phi(1);

r(2*N,1) = rs(N) - R(1); %TODO What condition to impose here?
% r(2*N,1) = dxidzeta1 * (3*rs(N) - 4*rs(N-1) + rs(N-2))/(2*delta1) - dxidzeta2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta2);




%##########################################################################
% Free Region

%boundary condition, dphi/dzeta = dPhi/dzeta
% r(2*N+1,1) = -dxidzeta1^2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta1) * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta1) * sin(Phi(1));
r(2*N+1,1) = dxidzeta1 * (3*phi(N) - 4*phi(N-1) + phi(N-2))/(2*delta1) - dxidzeta2 * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta2);
r(2*N+I+1,1)=dxidzeta1 * (3*rs(N) - 4*rs(N-1) + rs(N-2))/(2*delta1) - dxidzeta2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta2);

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

    [reqn, phieqn , ~, ~, ~, ~ ] = linearElasticity_Shell( R(jj), Rc1, Rc11, Phi(jj), Phic1, Phic11, zeta, P );
    
    r(2*N+jj, 1) = phieqn;
    r(2*N+I+jj, 1) = reqn;

end

%bc , Phi at zeta= pi/2 = pi/2
r(2*N+I,1)=Phi(N)-pi/2;
%boundary condition,  symmetry on dR/d\zeta = 0
r(2*N+2*I,1) = dxidzeta2*(3*R(N) - 4*R(N-1) + R(N-2))/(2*delta2);







% %##########################################################################
% %p
% %Ultimatly the volume constraint 
% [v] = volumn_of_region(Phi, R);
% r(2*N+2*I+1,1)= v - pi/4; %TODO
% % r(N+2*I+2,1)= H^2 * ( tan(X(N,1)) + pi/4 ) - pi/4; %TODO




% %##########################################################################
% %zeta_c
% 
% rc1=dxidzeta1*(3*rs(1) - 4*rs(2) + rs(3))/(2*delta1);
% % rc11=dxidzeta1^2 * (rs(jj+1) - 2*rs(jj) +rs(jj-1))/delta1^2;
% 
% phic1=dxidzeta1*(3*phi(1) - 4*phi(2) + phi(3))/(2*delta1);
% % phic11= dxidzeta1^2 *(phi(jj+1) - 2*phi(jj) +phi(jj-1))/delta1^2;
% 
% %e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
% e11= phic1*r(jj) - phi(jj)*rc1 - zeta*rc1 - zeta * phi(jj) - phic1 - rs(jj) +...
%     1/2*( (phi(jj))^2 + zeta^2 + (rs(jj))^2 + phic1^2 + rc1^2 );
% 
% % e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
% e22s = ( phi(jj) - zeta ) * rs(jj) * ( cos(zeta)/sin(zeta) ) + ...
%     1/2 * ( phi(jj) - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (rs(jj))^2 - 1 );
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
% %sigma^{11} = sigma^{11} at zeta_c
% r(2*N+2*I+1,1)= (sigma11_RHS-sigma11_LHS) ; %TODO
% % r(N+2*I+1,1)= sigma111-sigma11N; %TODO
% 
% %     function [v] = volumn_of_region(Phi, R)
% %         x_critical=R(1)*sin(Phi(1));
% %         z_critical=R(1)*cos(Phi(1));
% %         
% % %         if( le((abs(z_critical-H)),  10^(-8)))
% % %             fprintf('\n %.1e \t diff = %.4e  large \n', 1e-8, abs(z_critical-H));
% % %         end
% %         
% %         %volumn
% %         v=x_critical*z_critical + pi/4 * z_critical^2;
% % 
% %     end

% figure()
% n=2:N-1;
% i=N+1:N+I-2;
% 
% plot(n, xiValue1, 'ro-');
% hold on;
% plot(n, zetaValue1, 'bs-');
% 
% plot(i, xiValue2, 'ro-');
% plot(i, zetaValue2, 'bs-');
% legend('\xi 1', '\zeta 1', '\xi 2' , '\zeta 2','Location','northwest');


end

