function[X, zetac] = initialGuess_compressedSphere4(Np, H, zetac)
% fprintf('Started %s ... \t', mfilename);

N=round(Np/2); I=Np-N;
X=zeros([2*N+2*I+2,1]);
% fprintf('\n length(X) = %d \n', length(X));
% close all;
v=0;

% x_c= ( zetac);
% fprintf('just before loop v-v0 = %.2e \n', abs(v-pi/4))
while abs(v-pi/4) >= 0.01
%     fprintf('entered loop \n')
    x_c = H* zetac * (pi/2)/(pi/2 - zetac);
    for ii=1:N
        xpos(ii)=(ii-1)*(x_c)/(N-1);
        zpos(ii)=H;
        X(ii,1)=atan(xpos(ii)/H);
        X(N+ii,1)=sqrt( xpos(ii)^2 + H^2);
    %     X(2*N+ii,1)=1;
    end

    for ii=1:I
        Xpos(ii)=x_c+H*cos(pi/2*(1-(ii-1)/(I-1)));
        Zpos(ii)=H*sin(pi/2*(1-(ii-1)/(I-1)));
        X(2*N+ii,1) = atan(Xpos(ii)/Zpos(ii));
        X(2*N+I+ii,1) = sqrt( Xpos(ii)^2+Zpos(ii)^2 );
    end

    X(2*N+2*I+1,1)=0;
    X(2*N+2*I+2,1)=zetac;
    
    [v] = volumn_of_region(X(2*N+1:2*N+I), X(2*N+I+1:2*N+2*I));
%     fprintf('\n v = %.2e \t zetac = %.2e ', v, zetac);
    if v<pi/4
        zetac=zetac+0.001;
    else
        zetac=zetac-0.001;
    end
    
end


function [v] = volumn_of_region(Phi, R)
    x_critical=R(1)*sin(Phi(1));
    z_critical=R(1)*cos(Phi(1));

%         if( le((abs(z_critical-H)),  10^(-8)))
%             fprintf('\n %.1e \t diff = %.4e  large \n', 1e-8, abs(z_critical-H));
%         end

    %volumn
    v=x_critical*z_critical + pi/4 * z_critical^2;

end

% figure()
% plot(xpos,zpos, 'bs', 'MarkerFaceColor', 'b');
% hold on;
% plot(Xpos,Zpos, 'ro', 'MarkerFaceColor', 'r');
% xlim([0,1.2]);
% ylim([0,1]);
end