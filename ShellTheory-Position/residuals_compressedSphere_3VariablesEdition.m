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

Phi = X(3*N+1:3*N+I);

R = X(3*N+I+1:3*N+2*I);


P = X(3*N+2*I+1);
zetac = X(3*N+2*I+2);

delta1 = 1/(N-1);
delta2 =  1/(I-1);

dxidzeta1 = 1/zetac;
dxidzeta2 = 1/(pi/2-zetac);

r=zeros([3*N+2*I+2, 1]);

% rs=zeros([N, 1]);
% for ii=1:N
%     rs(ii)= H/(cos(X(ii)));
% end

sigma11N=0;
sigma111=0;

%##########################################################################
% Contact Region

% Boundary condition (first point in the first region): Symmetry dR/d\zeta = 0
% r(1, 1) = dxidzeta1*(-3*r(1) + 4*r(2) - r(3))/(2*delta1);
%Symmetry on y
r(1,1)=phi(1);
r(N+1, 1) = rs(1) - H;

for jj = 2:N-1; %Using central difference
%     fprintf('jj = %d \t ', jj);
    zeta = (jj-1)*delta1; %is actually xi
    
    rc1=dxidzeta1*(rs(jj+1) -rs(jj-1))/(2*delta1);
    rc11=dxidzeta1^2 * (rs(jj+1) - 2*rs(jj) +rs(jj-1))/delta1^2;
    
    phic1=dxidzeta1*(phi(jj+1) -phi(jj-1))/(2*delta1);
    phic11= dxidzeta1^2 *(phi(jj+1) - 2*phi(jj) +phi(jj-1))/delta1^2;
%     fprintf('rs(jj) = %.2e \t rc1 = %.2e \t phi(jj) = %.2e \t phic1 = %.2e \t rc1 = %.2e \n ', rs(jj), rc1, phi(jj), phic1, rc1);
%     fprintf('Phic1 = %.2e \t Rc1 = %.2e \t ', Phic1, Rc1);
    
    sqrtA = (( rc1 + zeta - phi(jj))^2 * (cos(zeta))^2 + ...
        (rs(jj) + phic1 - 1  )^2 - 2 *(rc1 + phi(jj) - ...
        zeta)*(rs(jj) + phic1 -1)*sin(zeta)*cos(zeta)  )^(1/2) *...
        sin(zeta)*(rs(jj) + (phi(jj) - zeta)* (cos(zeta)/sin(zeta)) );
    
    %e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
    e11= phic1*r(jj) - phi(jj)*rc1 - zeta*rc1 - zeta * phi(jj) - phic1 - rs(jj) +...
        1/2*( (phi(jj))^2 + zeta^2 + (rs(jj))^2 + phic1^2 + rc1^2 );
    
    %e_{11,1} = R_{,1} \left( R + R_{,11} -1 \right) - \Phi \left( R_{,11} +1 \right) +\Phi_{,1} \left(  \Phi + \Phi_{,11} - \zeta \right) + \Phi_{,11} \left( R - 1 \right) + \zeta
    e11c1 = rc1 * ( rs(jj) + rc11 - 1 ) - phi(jj) * (  rc11 + 1 ) + ...
        phic1 * ( phi(jj) + phic11 - zeta ) + phic11 * ( rs(jj) - 1 ) + zeta;
    
   % e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
    e22s = ( phi(jj) - zeta ) * rs(jj) * ( cos(zeta)/sin(zeta) ) + ...
        1/2 * ( phi(jj) - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (rs(jj))^2 - 1 );
    
    %\left( \frac{e_{22}}{\sin^2(\zeta)}  \right)_{,1}  
    %=   \left( \Phi_{,1}  - 1  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} +\left( \Phi  - \zeta  \right) R_{,1}  \frac{\cos(\zeta)}{\sin(\zeta)} - \left( \Phi  - \zeta  \right) R  \frac{1}{\sin^2(\zeta)} +  \left( \Phi  - \zeta  \right) \left( \Phi_{,1}  - 1  \right)  \frac{\cos^2(\zeta)}{\sin^2(\zeta)} -  \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos(\zeta)}{\sin^3(\zeta)}+R R_{,1}
    e22sc1 = ( phic1 - 1 ) * rs(jj) * ( cos(zeta)/sin(zeta) ) + ( phi(jj) - zeta ) * rc1 * ( cos(zeta)/sin(zeta) ) - ( phi(jj) - zeta ) * rs(jj) * ( 1/(sin(zeta))^2) + ( phi(jj) - zeta ) * ( phic1 - 1 ) * ( cos(zeta)/sin(zeta) )^2 - ( phi(jj) - zeta )^2 * ( cos(zeta)/(sin(zeta))^3 ) + rs(jj) * rc1;
    
    %\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
    sigma11 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
    if jj == N-1
        sigma11N=sigma11;
    end
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
%     fprintf('jj = %d \t zeta = %.2e \t  e11 = %.2e \t e22s =t %.2e \t sigma11 = %.2e \t sigma11c1 = %.2e \t sigma22s = %.2e  \n',  jj, zeta,  e11, e22s, sigma11, sigma11c1, sigma22s);
    
    %equation relating phi and r
    pContact = rs(jj)*cos(phi(jj)) - H;
    
    %equation 1 (for Phi)
    %  h  \left[\left[ \sigma^{11} \left( \Phi - \zeta - R_{,1}  \right)  + \sigma^{22} \left( \left( (\Phi - \zeta) \frac{ \cos^2(\zeta) }{ \sin^2(\zeta)} + R \frac{ \cos(\zeta) }{ \sin(\zeta)} \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(\Phi_{,1} + R - 1 \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p \ (\Phi- \zeta)  \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(jj, 1) = h * ( ( sigma11*( phi(jj) - zeta - rc1 ) + sigma22s * ( (phi(jj) - zeta)* ( (cos(zeta))^2/(sin(zeta))^2 ) + rs(jj) * ( (cos(zeta))/(sin(zeta) ) )) )* sin(zeta)  - sigma11c1 * ( phic1 + rs(jj) -1 ) * sin(zeta) - sigma11 * ( phic11 + rc1) * sin(zeta) - sigma11 * ( phic1 + rs(jj) -1 ) * cos(zeta) )   -pContact * ( Phi(jj) - zeta) * sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
    
    %equation 2 (for R)
    % h  \left[\left[ \sigma^{11} \left( \Phi_{,1} + R - 1  \right)  +   \sigma^{22} \left( \left( R + (\Phi - \zeta) \frac{ \cos (\zeta) }{ \sin (\zeta)}   \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(R_{,1} + \zeta - \Phi  \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p  R \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(N+jj, 1) =  h * (  ( sigma11 * ( phic1 + rs(jj) - 1 ) + sigma22s * ( rs(jj) + ( phi(jj) - zeta ) * ( cos(zeta) / sin(zeta) ) ) ) * sin(zeta) - sigma11c1 * (rc1 + zeta - phi(jj) ) * sin(zeta) - sigma11 * (rc11 + 1 - phic1)  * sin(zeta) - sigma11 * (rc1 + zeta - phi(jj) ) * cos(zeta) ) - pContact * R (jj) *sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
    
    
end
%Boundary Conidtion - continuity
r(N,1) = phi(N)-Phi(1);

r(2*N,1) = rs(N)*cos(phi(N)) - H; %TODO What condition to impose here?




%##########################################################################
% Free Region

%boundary condition, dZ/dzeta = 0
r(3*N+1,1) = -dxidzeta1^2 * (-3*R(1) + 4*R(2) - R(3))/(2*delta1) * (-3*Phi(1) + 4*Phi(2) - Phi(3))/(2*delta1) * sin(Phi(1));
r(3*N+I+1,1)=rs(1)-R(1);

for jj = 2:I-1; %Using central difference
%     fprintf('jj = %d \t ', jj);
    zeta = (jj-1)*delta2; %is actually xi
    
    Rc1=dxidzeta2*(R(jj+1) -R(jj-1))/(2*delta2);
    Rc11=dxidzeta2^2 * (R(jj+1) - 2*R(jj) +R(jj-1))/delta2^2;
    
    Phic1=dxidzeta2*(Phi(jj+1) -Phi(jj-1))/(2*delta2);
    Phic11= dxidzeta2^2 *(Phi(jj+1) - 2*Phi(jj) +Phi(jj-1))/delta2^2;
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
    if jj == 2
        sigma111=sigma11;
    end
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
    r(3*N+jj, 1) = h * ( ( sigma11*( Phi(jj) - zeta - Rc1 ) + sigma22s * ( (Phi(jj) - zeta)* ( (cos(zeta))^2/(sin(zeta))^2 ) + R(jj) * ( (cos(zeta))/(sin(zeta) ) )) )* sin(zeta)  - sigma11c1 * ( Phic1 + R(jj) -1 ) * sin(zeta) - sigma11 * ( Phic11 + Rc1) * sin(zeta) - sigma11 * ( Phic1 + R(jj) -1 ) * cos(zeta) )  - P * ( Phi(jj) - zeta) * sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
    
    %equation 2 (for R)
    % h  \left[\left[ \sigma^{11} \left( \Phi_{,1} + R - 1  \right)  +   \sigma^{22} \left( \left( R + (\Phi - \zeta) \frac{ \cos (\zeta) }{ \sin (\zeta)}   \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)     -  \left\{ \left( \sigma^{11} \left(R_{,1} + \zeta - \Phi  \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]    -   \frac{p  R \sqrt{A}}{\sqrt{R^2 + \Phi^2 + \zeta^2 - 2 \Phi \zeta}}   = 0
    r(3*N+I+jj, 1) =  h * (  ( sigma11 * ( Phic1 + R(jj) - 1 ) + sigma22s * ( R(jj) + ( Phi(jj) - zeta ) * ( cos(zeta) / sin(zeta) ) ) ) * sin(zeta) - sigma11c1 * (Rc1 + zeta - Phi(jj) ) * sin(zeta) - sigma11 * (Rc11 + 1 - Phic1)  * sin(zeta) - sigma11 * (Rc1 + zeta - Phi(jj) ) * cos(zeta) )   - P * R (jj) *sqrtA/( (R(jj))^2 + (Phi(jj))^2 + zeta^2 - 2 * Phi(jj) * zeta ) ;
end

%bc , Phi at zeta= pi/2 = pi/2
r(3*N+I,1)=Phi(N)-pi/2;
%boundary condition,  symmetry on dR/d\zeta = 0
r(3*N+2*I,1) = dxidzeta2*(3*R(N) - 4*R(N-1) + R(N-2))/(2*delta2);


%##########################################################################
%zeta_c

rc1=dxidzeta1*(3*rs(1) - 4*rs(2) + rs(3))/(2*delta1);
% rc11=dxidzeta1^2 * (rs(jj+1) - 2*rs(jj) +rs(jj-1))/delta1^2;

phic1=dxidzeta1*(3*phi(1) - 4*phi(2) + phi(3))/(2*delta1);
% phic11= dxidzeta1^2 *(phi(jj+1) - 2*phi(jj) +phi(jj-1))/delta1^2;

%e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
e11= phic1*r(jj) - phi(jj)*rc1 - zeta*rc1 - zeta * phi(jj) - phic1 - rs(jj) +...
    1/2*( (phi(jj))^2 + zeta^2 + (rs(jj))^2 + phic1^2 + rc1^2 );

% e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
e22s = ( phi(jj) - zeta ) * rs(jj) * ( cos(zeta)/sin(zeta) ) + ...
    1/2 * ( phi(jj) - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (rs(jj))^2 - 1 );

%\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
sigma11_LHS = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );


Rc1=dxidzeta2*(-3*R(1) + 4*R(2) - R(3))/(2*delta2);

Phic1=dxidzeta2*(-3*Phi(N) + 4*Phi(N-1) - Phi(N-2))/(2*delta2); %(Phi(jj+1) -Phi(jj-1))/(2*delta2);

%e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
e11= Phic1*R(1) - Phi(1)*Rc1 - zetac*Rc1 - zeta * Phi(1) - Phic1 - R(1) +...
    1/2*( (Phi(1))^2 + zetac^2 + (R(1))^2 + Phic1^2 + Rc1^2 );

% e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
e22s = ( Phi(1) - zetac ) * R(jj) * ( cos(zetac)/sin(zetac) ) + ...
    1/2 * ( Phi(1) - zetac)^2 * ( cos(zetac)/sin(zetac) )^2 + 1/2 * ( (R(1))^2 - 1 );

%\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
sigma11_RHS = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );




%##########################################################################
%p
%Ultimatly the volume constraint but not for now. 
[v] = volumn_of_region(Phi, R);
r(3*N+2*I+1,1)= v - pi/4; %TODO
% r(N+2*I+2,1)= H^2 * ( tan(X(N,1)) + pi/4 ) - pi/4; %TODO

%##########################################################################
%zetac
%sigma^{11} = sigma^{11} at zeta_c
r(3*N+2*I+2,1)= (sigma11_RHS-sigma11_LHS) ; %TODO
% r(N+2*I+1,1)= sigma111-sigma11N; %TODO

    function [v] = volumn_of_region(Phi, R)
        x_critical=R(1)*sin(Phi(1));
        z_critical=R(1)*cos(Phi(1));
        
%         if( le((abs(z_critical-H)),  10^(-8)))
%             fprintf('\n %.1e \t diff = %.4e  large \n', 1e-8, abs(z_critical-H));
%         end
        
        %volumn
        v=x_critical*z_critical + pi/4 * z_critical^2;

    end

end

