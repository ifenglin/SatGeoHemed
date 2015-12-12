function ret = t2E(t,t0,a,e)
% This function returns Eccevtric anomaly of an object circulating around
% the Earth in a specific time t and time in pregee and a and e
GM = 398600.44 * 10^9; % Meter/second_square
n=sqrt(GM/a^3);
M = (t-t0)*n;
E0 = M;
E = M + e*sin(E0);
while(abs(E-E0)>10^-3)
    E0=E;
    E = M + e*sin(E0)
end% end of While loop
ret = E;    

