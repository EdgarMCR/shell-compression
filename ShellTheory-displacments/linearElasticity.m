function [u1eqn, u3eqn , sigma11, sigma11c1, sigma22s, e11 ] = linearElasticity( u3, u3c1, u3c11, u1, u1c1, u1c11, zeta, P )
%Evaluates goiverning equations given the displacments and their
%derivatives
%   Detailed explanation goes here

global E; %Youngs Modulus
global nu; %Poission ratio
global h; %membran thickness

if E==0 || nu==0 || h==0
    warning('\n E or nu or h is zero! E = %f, nu=%f, h=%f.\n', E, nu, h);
end

sqrtA = (( u3c1 - u1)^2 * (cos(zeta))^2 + ...
    (1 + u1c1 + u3)^2 - 2 *(u3c1 - u1)*(1 + u1c1 + u3)*sin(zeta)*cos(zeta)  )^(1/2) *... 
    sin(zeta)*(1+(u1)* (cos(zeta)/sin(zeta)) + u3 );

%e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
e11= u1c1 + u3 + u1c1*u3 - u1*u3c1 +...
    1/2*( (u1)^2 + u1c1^2 + (u3)^2 + u3c1^2 );

%e_{11,1} = R_{,1} \left( R + R_{,11} -1 \right) - \Phi \left( R_{,11} +1 \right) +\Phi_{,1} \left(  \Phi + \Phi_{,11} - \zeta \right) + \Phi_{,11} \left( R - 1 \right) + \zeta
e11c1 = u1c11 + u3c1 + u1c11*u3 - u1 * u3c11 + u1 * u1c1 + u1c1*u1c11 + u3 * u3c1 + u3c1*u3c11;

% e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
e22s = (u1c1 + u1c1*u3 + u1*u3c1)* (cos(zeta)/sin(zeta)) - (u1 +u1*u3)* (1/(sin(zeta))^2) + u3c1 + u1 * u1c1 * ((cos(zeta))^2/(sin(zeta))^2) - (u1)^2 * ((cos(zeta))/(sin(zeta))^3) + u3 * u3c1 ;

%\left( \frac{e_{22}}{\sin^2(\zeta)}  \right)_{,1}  
%=   \left( \Phi_{,1}  - 1  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} +\left( \Phi  - \zeta  \right) R_{,1}  \frac{\cos(\zeta)}{\sin(\zeta)} - \left( \Phi  - \zeta  \right) R  \frac{1}{\sin^2(\zeta)} +  \left( \Phi  - \zeta  \right) \left( \Phi_{,1}  - 1  \right)  \frac{\cos^2(\zeta)}{\sin^2(\zeta)} -  \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos(\zeta)}{\sin^3(\zeta)}+R R_{,1}
e22sc1 = ( u1c1 + u1c1 * u3 + u1*u3c1 ) * ( cos(zeta)/sin(zeta) ) - ( u1 + u1* u3 ) * ( 1/(sin(zeta))^2 ) + ( u1 * u1c1 ) * ( cos(zeta)/sin(zeta) )^2 - ( u1 )^2 * ( cos(zeta)/(sin(zeta))^3 ) + u3 * u3c1;

%\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
sigma11 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
sigma11c1 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11c1 + ( nu/(1-nu)  ) * e22sc1 );

%\sigma^{22} = \frac{h \ E}{(1+\nu)} \left[  \left(  \frac{ \nu}{1- \nu}   \right) e_{11} +   \left( 1 + \frac{ \nu}{1- \nu}  \right) \frac{1}{\sin^2(\zeta)}  e_{22} \right] \frac{1}{\sin^2(\zeta)}
sigma22s = ( E / (1+nu)) * (  nu/(1-nu) * e11 + ( 1+ nu/(1-nu) ) * e22s  );



%governing equations
%h  \left[\left[ \sigma^{11} \left( u^1 - u^3_{,1}   \right)  + \sigma^{22} \left( \left( \frac{ \cos(\zeta) }{ \sin(\zeta)} + u^1 \frac{ \cos^2(\zeta) }{ \sin^2(\zeta)} + u^3 \frac{ \cos(\zeta) }{ \sin(\zeta)} \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)   -  \left\{ \left( \sigma^{11} \left( 1 + u^1_{,1} + u^3 \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]     -   \frac{p \ u^1   \sqrt{A}}{\sqrt{\left(1 + u^3\right)^2 + \left(u^1\right)^2}}   = 0\label{eq:GovEqR}
u1eqn   = h * ( ( sigma11*( u1-u3c1 ) + sigma22s * (  ( u1 * (cos(zeta))^2/(sin(zeta))^2 ) + ((u3 +1 ) * (cos(zeta))/(sin(zeta) ) )) )*sin(zeta) - sigma11c1 * ( 1 + u1c1 + u3 ) * sin(zeta) - sigma11 * ( u1c11 + u3c1 ) * sin(zeta) - sigma11 * ( 1 + u1c1 + u3  ) * cos(zeta) )  - P * ( u1 ) * sqrtA/( (1+ u3 )^2 + ( u1 )^2) ;

%h  \left[\left[ \sigma^{11} \left( 1 + u^1_{,1} + u^3  \right)  + \sigma^{22} \left( \left( 1 + u^1 \frac{ \cos (\zeta) }{ \sin (\zeta)} + u^3  \right) \sin^2(\zeta)  \right) \right] \sin(\zeta)-  \left\{ \left( \sigma^{11} \left( - \left( u^1 - u^3_{,1} \right) \right)   \right)  sin(\zeta) \right\}_{,1} \ \right]     -   \frac{p  \left( 1 + u^3 \right)  \sqrt{A}}{\sqrt{\left(1 + u^3\right)^2 + \left(u^1\right)^2}}   = 0
u3eqn = h * (  ( sigma11 * ( 1 + u1c1 + u3 ) + sigma22s * ( 1 + u3 + u1 * ( cos(zeta) / sin(zeta) ) ) ) * sin(zeta) + sigma11c1 * (u1 - u3c1) * sin(zeta) + sigma11 * (u1c1 - u3c11)  * sin(zeta) + sigma11 * (u1 - u3c1 ) * cos(zeta) )   - P * (1+ u3 ) *sqrtA/( (1+ u3 )^2 + ( u1 )^2) ;

end

