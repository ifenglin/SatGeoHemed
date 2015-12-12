function E2 = E_from_t( t,to,a,e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    GM = 398600.44 * 10^9;
    n=sqrt(GM/a^3);
    M = (t-to)*n;
    E2 = M;
    E1 = M+1;
    while(abs(E2-E1)<10^-5)
        E1 = E2;
        E2 = M + e*sin(E1);
    end

end

