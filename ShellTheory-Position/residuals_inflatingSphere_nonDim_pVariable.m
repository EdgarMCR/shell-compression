function [r] = residuals_inflatingSphere_nonDim(X)
 % The number of discrete points in the first region
global N;
% fprintf('N=%d \n', N);

% Parameters:
%pressure inside region
global P; 

global E; %Youngs Modulus
global nu; %Poission ratio

global h; %membran thickness
% the distance between the discrete points first region
delta =  1/(N-1); 

dxidzeta=1/(pi/2);

Phi = X(1:N);
R = X(N+1:2*N);
p=X(2*N+1,1);


r=zeros([2*N+1, 1]);



% Boundary condition (first point in the first region): Symmetry dR/d\zeta = 0
r(1, 1) = dxidzeta*(-3*R(1) + 4*R(2) - R(3))/(2*delta);
%Symmetry on y
r(N+1,1)=Phi(1);

for jj = 2:N-1; %Using central difference
%     fprintf('jj = %d \t ', jj);
    xi = (jj-1)*delta; %is actually xi
    zeta=xi*pi/2;
    
%     p=-P;
%     p=-P*cos(2*zeta);
%     p=P+P*cos(2*zeta);
    
    zetaValue(jj)=zeta;
    pValue(jj)=p;
    
    Rc1=dxidzeta*(R(jj+1) -R(jj-1))/(2*delta);
    Rc11=dxidzeta^2 * (R(jj+1) - 2*R(jj) +R(jj-1))/delta^2;
    
    Phic1=dxidzeta*(Phi(jj+1) -Phi(jj-1))/(2*delta);
    Phic11= dxidzeta^2 *(Phi(jj+1) - 2*Phi(jj) +Phi(jj-1))/delta^2;
%     fprintf('R(jj) = %.2e \t Rc1 = %.2e \t Phi(jj) = %.2e \t Phic1 = %.2e \n ', R(jj), Rc1, Phi(jj), Phic1);
%     fprintf('Phic1 = %.2e \t Rc1 = %.2e \t ', Phic1, Rc1);
    
    sqrtA = (( Rc1 + zeta - Phi(jj))^2 * (cos(zeta))^2 + ...
        (R(jj) + Phic1 - 1  )^2 - 2 *(Rc1 + Phi(jj) - ...
        zeta)*(R(jj) + Phic1 -1)*sin(zeta)*cos(zeta)  )^(1/2) *...
        sin(zeta)*(R(jj) + (Phi(jj) - zeta)* (cos(zeta)/sin(zeta)) );
    
    %e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
    e11= Phic1*R(jj) - Phi(jj)*Rc1 - zeta*Rc1 - zeta * Phi(jj) - Phic1 - R(jj) +...
        1/2*( (Phi(jj))^2 + zeta^2 + (R(jj))^2 + Phic1^2 + Rc1^2 );
    
    %e_{11,1} = R_{,1} \left( R + R_{,11} -1 \right) - \Phi \left( R_{,11} +1 \right) +\Phi_{,1} \left(  \Phi + \Phi_{,11} - \zeta \right) + \Phi_{,11} \left( R - 1 \right) + \zeta
    e11c1 = Rc1 * ( R(jj) + Rc11 - 1 ) - Phi(jj) * (  Rc11 + 1 ) + ...
        Phic1 * ( Phi(jj) + Phic11 - zeta ) + Phic11 * ( R(jj) - 1 ) + zeta;
    
   % e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
    e22s = ( Phi(jj) - zeta ) * R(jj) * ( cos(zeta)/sin(zeta) ) + ...
        1/2 * ( Phi(jj) - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (R(jj))^2 - 1 );
    
    %\left( \frac{e_{22}}{\sin^2(\zeta)}  \right)_{,1}  
    %=   \left( \Phi_{,1}  - 1  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} +\left( \Phi  - \zeta  \right) R_{,1}  \frac{\cos(\zeta)}{\sin(\zeta)} - \left( \Phi  - \zeta  \right) R  \frac{1}{\sin^2(\zeta)} +  \left( \Phi  - \zeta  \right) \left( \Phi_{,1}  - 1  \right)  \frac{\cos^2(\zeta)}{\sin^2(\zeta)} -  \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos(\zeta)}{\sin^3(\zeta)}+R R_{,1}
    e22sc1 = ( Phic1 - 1 ) * R(jj) * ( cos(zeta)/sin(zeta) ) + ( Phi(jj) - zeta ) * Rc1 * ( cos(zeta)/sin(zeta) ) - ( Phi(jj) - zeta ) * R(jj) * ( 1/(sin(zeta))^2) + ( Phi(jj) - zeta ) * ( Phic1 - 1 ) * ( cos(zeta)/sin(zeta) )^2 - ( Phi(jj) - zeta )^2 * ( cos(zeta)/(sin(zeta))^3 ) + R(jj) * Rc1;
    
    %\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
    sigma11 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
    sigma11c1 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11c1 + ( nu/(1-nu)  ) * e22sc1 );

    %\sigma^{22} = \frac{h \ E}{(1+\nu)} \left[  \left(  \frac{ \nu}{1- \nu}   \right) e_{11} +   \left( 1 + \frac{ \nu}{1- \nu}  \right) \frac{1}{\sin^2(\zeta)}  e_{22} \right] \frac{1}{\sin^2(\zeta)}
    sigma22s = ( E / (1+nu)) * (  nu/(1-nu) * e11 + ( 1+ nu/(1-nu) ) * e22s  );
    
%     fprintf('size(sqrtA) = %d, %d \t', size(sqrtA));
%     fprintf('size(e11) = %d, %d \t', size(e11));
%     fprintf('size(e22s) = %d, %d \t', size(e22s));
%     fprintf('size(e11c1) = %d, %d \t', size(e11c1));
%     fprintf('size(e22sc1) = %d, %d \t', size(e22sc1));
%     fprintf('size(sigma11) = %d, %d \t', size(sigma11));
%     fprintf('size(sigma11c1) = %d, %d \t', size(sigma11c1));
%     fprintf('size(sigma22s) = %d, %d \n', size(sigma22s));
%     fprintf('jj = %d \t zeta = %.2e \t sqrtA = %.2e \t e11 = %.2e \t e22s = %.2e \t sigma11 = %.2e \t sigma11c1 = %.2e \t sigma22s = %.2e  \n',  jj, zeta, sqrtA, e11, e22s, sigma11, sigma11c1, sigma22s);

    %equation 1 (for Phi)
    %  h  \left[\left[ \sigma^{11} \left( \Phi - \zeta - R_{,1}  \right)  + \sigma^{22} \left( \left( (\Phi - \zeta) \frac{ \cos^2(\zeta) }{ \sin^2(\zeta)} + R \frac{ \cos(\zeta) }{ \sin(\zeta)} \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(\Phi_{,1} + R - 1 \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p \ (\Phi- \zeta)  \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(jj, 1) = h * ( ( sigma11*( Phi(jj) - zeta - Rc1 ) + sigma22s * ( (Phi(jj) - zeta)* ( (cos(zeta))^2/(sin(zeta))^2 ) + R(jj) * ( (cos(zeta))/(sin(zeta) ) )) )* sin(zeta)  - sigma11c1 * ( Phic1 + R(jj) -1 ) * sin(zeta) - sigma11 * ( Phic11 + Rc1) * sin(zeta) - sigma11 * ( Phic1 + R(jj) -1 ) * cos(zeta) )  - p * ( Phi(jj) - zeta) * sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
    
    %equation 2 (for R)
    % h  \left[\left[ \sigma^{11} \left( \Phi_{,1} + R - 1  \right)  +   \sigma^{22} \left( \left( R + (\Phi - \zeta) \frac{ \cos (\zeta) }{ \sin (\zeta)}   \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(R_{,1} + \zeta - \Phi  \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p  R \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(N+jj, 1) =  h * (  ( sigma11 * ( Phic1 + R(jj) - 1 ) + sigma22s * ( R(jj) + ( Phi(jj) - zeta ) * ( cos(zeta) / sin(zeta) ) ) ) * sin(zeta) - sigma11c1 * (Rc1 + zeta - Phi(jj) ) * sin(zeta) - sigma11 * (Rc11 + 1 - Phic1)  * sin(zeta) - sigma11 * (Rc1 + zeta - Phi(jj) ) * cos(zeta) )   - p * R (jj) *sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
end

%boundary condition,  symmetry on dR/d\zeta = 0
r(N,1) = dxidzeta*(3*R(N) - 4*R(N-1) + R(N-2))/(2*delta);
%bc , Phi at zeta= pi/2 = pi/2
r(2*N,1)=Phi(N)-pi/2;

r(2*N+1,1)=P-p;
% h1=figure();
% plot(zetaValue,pValue); 
end

