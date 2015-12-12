function ret = Sat2loc(a,e,inc,OMEGA,omega,t,t0,X_Obs,Y_Obs,Z_Obs)
% This function returns the position[x;y;z] of the Sat. in observertory 
% with the given position.
% a: semi-majour axis
% e: eccentricity
% inc: inclination of the satellite's orbit
% t: time
% t0: time of being in the pregee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


observer_pos = [X_Obs;Y_Obs;Z_Obs];
GM = 398600.44 * 10^9; % Meter/second_square
omega_earth = 2 * pi / 86164; % 1/second
T = 2*pi*sqrt(a^3/GM);
E = t2E(t,t0,a,e);
temp = xyz2blh(observer_pos(1),observer_pos(2),observer_pos(3),'grs80');
phi_observer = temp(1)*pi/180;
lambda_observer = temp(2)*pi/180;
p = a * (1 - e^2);
c = sqrt(p*GM);
K = [cos(OMEGA); sin(OMEGA); 0];
C0 = [sin(OMEGA)*sin(inc); -cos(OMEGA)*sin(inc); cos(inc)];
P = cos(omega)*K + sin(omega)*(cross(C0,K));
Q = -sin(omega)*K + cos(omega)*(cross(C0,K));
M2 = [1  0  0 ; 0  -1  0 ; 0  0  1];
f = 2 * atan( tan(E/2 ) * sqrt( (1+e)/(1+e) ) );
if(f<0)
    f= f + 2*pi;
end % end of if
r = a*(1-e^2)/( 1+e*cos(f) );
Rp = r * ( cos(f)*P +  sin(f)*Q );
Rp_rot =( Rp'*rotz(t*omega_earth))';
R_local =  roty(pi/2 - phi_observer) * M2 * ...
    rotz(lambda_observer-pi) * (Rp_rot-observer_pos );
% R_loc_pol = Azz(R_local);

ret = R_local;
