function [r] = residuals_inflation(X)
 % The number of discrete points in the first region
global N;

% Parameters:
%pressure inside region
global P;

u1 = X(1:N);
u3 = X(N+1:2*N);

delta1 = 1/(N-1);
dxidzeta1 = 1/(pi/2);
r=zeros([2*N, 1]);

%##########################################################################
% Contact Region

% Boundary condition (first point in the first region): Symmetry dR/d\zeta = 0
r(1,1)=u1(1); 
r(N+1,1)=dxidzeta1 * (-3*u3(1) + 4*u3(2) - u3(3))/(2*delta1); %Symmetry on y

for ii = 2:N-1; %Using central difference
    xi = (ii-1)*delta1; %is actually xi
    zeta= xi*pi/2;

    u1c1=dxidzeta1*(u1(ii+1) -u1(ii-1))/(2*delta1);
    u1c11=dxidzeta1^2 * (u1(ii+1) - 2*u1(ii) +u1(ii-1))/delta1^2;
    
    u3c1=dxidzeta1*(u3(ii+1) -u3(ii-1))/(2*delta1);
    u3c11= dxidzeta1^2 *(u3(ii+1) - 2*u3(ii) +u3(ii-1))/delta1^2;

    [reqn, phieqn , ~ , ~, ~, ~] = linearElasticity( u3(ii), u3c1, u3c11, u1(ii), u1c1, u1c11, zeta,  P);
	
    %equation 1 (for Phi)
    r(ii, 1) = phieqn;
    
    %equation 2 (for R)
    r(N+ii, 1) = reqn; 
end

%bc , Phi at zeta= pi/2 = pi/2
r(N,1)=u1(N);
%boundary condition,  symmetry on dR/d\zeta = 0
r(2*N,1) = dxidzeta1*(3*u3(N) - 4*u3(N-1) + u3(N-2))/(2*delta1);

end

