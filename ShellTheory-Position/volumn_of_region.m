function [v] = volumn_of_region(rs, phi, R, Phi, I, N)
    for iii=1:N
        z1(iii) = rs(iii)*cos(phi(iii));
        x1(iii) = rs(iii)*sin(phi(iii));
    end
    Q1 = trapz(x1,z1);
    
    for jjj=1:I
        z2(jjj) = R(jjj)*cos(Phi(jjj));
        x2(jjj) = R(jjj)*sin(Phi(jjj));
    end
    Q2 = trapz(x2,z2);
    v=Q1+Q2;
end

