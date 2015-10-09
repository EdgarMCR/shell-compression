function [reqn, phieqn , sigma11, sigma11c1, sigma22s, e11 ] = linearElasticity( R, Rc1, Rc11, Phi, Phic1, Phic11, zeta, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global E; %Youngs Modulus
global nu; %Poission ratio
global h; %membran thickness

if E==0 || nu==0 || h==0
    warning('\n E or nu or h is zero! E = %f, nu=%f, h=%f.\n', E, nu, h);
end

sqrtA = (( Rc1 + zeta - Phi)^2 * (cos(zeta))^2 + ...
    (R + Phic1 - 1  )^2 - 2 *(Rc1 + Phi - ...
    zeta)*(R + Phic1 -1)*sin(zeta)*cos(zeta)  )^(1/2) *...
    sin(zeta)*(R + (Phi - zeta)* (cos(zeta)/sin(zeta)) );

%e_{11} = \Phi_{,1} R - \Phi R_{,1} - \zeta R_{,1} - \zeta \Phi -  \Phi_{,1} - R + \frac{1}{2} \left[ \Phi^2 + \zeta^2 + R^2 + \left( \Phi_{,1} \right)^2 + \left( R_{,1} \right)^2\right]
e11= Phic1*R - Phi*Rc1 - zeta*Rc1 - zeta * Phi - Phic1 - R +...
    1/2*( (Phi)^2 + zeta^2 + (R)^2 + Phic1^2 + Rc1^2 );

%e_{11,1} = R_{,1} \left( R + R_{,11} -1 \right) - \Phi \left( R_{,11} +1 \right) +\Phi_{,1} \left(  \Phi + \Phi_{,11} - \zeta \right) + \Phi_{,11} \left( R - 1 \right) + \zeta
e11c1 = Rc1 * ( R + Rc11 - 1 ) - Phi * (  Rc11 + 1 ) + ...
    Phic1 * ( Phi + Phic11 - zeta ) + Phic11 * ( R - 1 ) + zeta;

% e_{22} = \left( \Phi  - \zeta  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} + \frac{1}{2} \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos^2(\zeta)}{\sin^2(\zeta)} +\frac{1}{2} \left(   R^2 - 1  \right) 
e22s = ( Phi - zeta ) * R * ( cos(zeta)/sin(zeta) ) + ...
    1/2 * ( Phi - zeta)^2 * ( cos(zeta)/sin(zeta) )^2 + 1/2 * ( (R)^2 - 1 );

%\left( \frac{e_{22}}{\sin^2(\zeta)}  \right)_{,1}  
%=   \left( \Phi_{,1}  - 1  \right) R  \frac{\cos(\zeta)}{\sin(\zeta)} +\left( \Phi  - \zeta  \right) R_{,1}  \frac{\cos(\zeta)}{\sin(\zeta)} - \left( \Phi  - \zeta  \right) R  \frac{1}{\sin^2(\zeta)} +  \left( \Phi  - \zeta  \right) \left( \Phi_{,1}  - 1  \right)  \frac{\cos^2(\zeta)}{\sin^2(\zeta)} -  \left( \left( \Phi  - \zeta  \right) \right)^2 \frac{\cos(\zeta)}{\sin^3(\zeta)}+R R_{,1}
e22sc1 = ( Phic1 - 1 ) * R * ( cos(zeta)/sin(zeta) ) + ( Phi - zeta ) * Rc1 * ( cos(zeta)/sin(zeta) ) - ( Phi - zeta ) * R * ( 1/(sin(zeta))^2) + ( Phi - zeta ) * ( Phic1 - 1 ) * ( cos(zeta)/sin(zeta) )^2 - ( Phi - zeta )^2 * ( cos(zeta)/(sin(zeta))^3 ) + R * Rc1;

%\sigma^{11} = h \left( \frac{E}{(1+\nu)} \left( 1 + \frac{ \nu}{1- \nu}  \right) \right) e_{11} + h \left( \frac{E}{(1+\nu)} \left(  \frac{ \nu}{1- \nu}  \right)  \right) \frac{1}{\sin^2(\zeta)} e_{22}
sigma11 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11 + ( nu/(1-nu)  ) * e22s  );
sigma11c1 = ( E / (1+nu)) *  ( ( ( 1 + nu/(1-nu)) ) * e11c1 + ( nu/(1-nu)  ) * e22sc1 );

%\sigma^{22} = \frac{h \ E}{(1+\nu)} \left[  \left(  \frac{ \nu}{1- \nu}   \right) e_{11} +   \left( 1 + \frac{ \nu}{1- \nu}  \right) \frac{1}{\sin^2(\zeta)}  e_{22} \right] \frac{1}{\sin^2(\zeta)}
sigma22s = ( E / (1+nu)) * (  nu/(1-nu) * e11 + ( 1+ nu/(1-nu) ) * e22s  );



%governing equations

reqn   = h * ( ( sigma11*( Phi - zeta - Rc1 ) + sigma22s * ( (Phi - zeta)* ( (cos(zeta))^2/(sin(zeta))^2 ) + R * ( (cos(zeta))/(sin(zeta) ) )) )* sin(zeta)  - sigma11c1 * ( Phic1 + R -1 ) * sin(zeta) - sigma11 * ( Phic11 + Rc1) * sin(zeta) - sigma11 * ( Phic1 + R -1 ) * cos(zeta) )  - P * ( Phi - zeta) * sqrtA/( (R)^2 + (Phi)^2 + zeta^2 - 2 * Phi * zeta ) ;
phieqn = h * (  ( sigma11 * ( Phic1 + R - 1 ) + sigma22s * ( R + ( Phi - zeta ) * ( cos(zeta) / sin(zeta) ) ) ) * sin(zeta) - sigma11c1 * (Rc1 + zeta - Phi ) * sin(zeta) - sigma11 * (Rc11 + 1 - Phic1)  * sin(zeta) - sigma11 * (Rc1 + zeta - Phi ) * cos(zeta) )   - P * R *sqrtA/( (R)^2 + (Phi)^2 + zeta^2 - 2 * Phi * zeta ) ;

end

