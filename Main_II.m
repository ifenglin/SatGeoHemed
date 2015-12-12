 %%%   Given Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear; clc; close all
format long
% Sats_name = ['ERSA','GPS','Beidou','Champ','Molniya'];
sat_parameters = [
    6838000   ,  0.004 ,  deg2rad(87)  ,  deg2rad(0)   ,  deg2rad(0); 
    6838000   ,  0.004 ,  deg2rad(87)  ,  deg2rad(0)   ,  deg2rad(0)];
GM = 398600.44 * 10^9; % Meter/second_square
sat_parameters(3,1) = (86164^2*GM/(4*pi^2))^(1/3);
omega_earth = 2 * pi / 86164; % 1/second
n = 720; % number of epochs
[n_sat , n_par] = size(sat_parameters);
%%%  calculting the Period of each Satellite T
for (i=1:n_sat)
    T(i,1) = 2*pi*sqrt(sat_parameters(i,1)^3/GM);
end % end of for-loop for i
colors = ['.r'; '.g'; '.b'; '.c'; '.m'];
pos_berlin = [3783.26649 ; 901.64960 ; 5035.24814]*1000;
whatE = linspace(0,4*pi,n)';
temp = xyz2blh(pos_berlin(1),pos_berlin(2),pos_berlin(3));
phi_berlin = temp(1)*pi/180;
lambda_berlin = temp(2)*pi/180;
h_berlin = temp(3);
a_earth = 6378140;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i=1)% i is the index of satellite
    a = sat_parameters(i,1);
    e = sat_parameters(i,2);
    t(:,i) = (whatE-e*sin(whatE))*T(i)/(2*pi);
    inc = sat_parameters(i,3);% angle i
    OMEGA = sat_parameters(i,4);
    omega = sat_parameters(i,5);
    J2 = 0.00108263;
    nn = sqrt(GM/a^3);
    dt = 2 * T(i) / n;
    diff_space_fixed_p = zeros(3,n); 
    diff_earth_fixed_L = zeros(n); 
    diff_earth_fixed_B = zeros(n);
    diff_topocentric_fixed_p = zeros(3,n);
    diff_topocentric_fixed_v = zeros(3,n);
    diff_topocentric_fixed_Az = zeros(n);
    diff_topocentric_fixed_El = zeros(n);
    for j=1:n% j is index of epoch
        delta_t = j * dt;
        d_omega = (-1 * 3 * nn * J2 * a_earth^2 ) * (1 - 5 * cos(inc)^2 ) / ( 4 * ( 1 - e^2 ) ^ 2 * a^2 );
        d_OMEGA = (-1 * 3 * nn * J2 * a_earth^2 ) * cos(inc) / ( 2 * ( 1 - e^2 ) ^ 2 * a^2 );
        d_M = nn + (-1 * 3 * nn * J2 * a_earth^2 ) * ( 3 * cos(i)^2 - 1 ) / ( 4 * ( 1 - e^2 ) ^ (3/2) * a^2);

        disp(j);
        disp(rad2deg(omega));
        disp(rad2deg(OMEGA));

        p_a = a;
        p_e = e;
        p_M = j + d_M * delta_t
        p_omega = rad2deg(omega + d_omega * delta_t)
        p_OMEGA = rad2deg(OMEGA + d_OMEGA * delta_t)
        p_i = rad2deg(inc);
        
        p_E = 0;
        E0 = j;
        E(1) = j + e * sin(E0);
        for k = 1:50
            E(k+1) = j + e * sin(E(k));
            if abs(E(k+1) - E(k)) < 10^-8
                break
            end 
            p_E = E(k+1)
        end 
        
        [unperturbed_space_fixed_p,unperturbed_space_fixed_v] = position_velocity(a,e,i,OMEGA,omega, j);
        [space_fixed_p,space_fixed_v] = position_velocity(p_a,p_e,p_i,p_OMEGA,p_omega, p_E);
        diff_space_fixed_p(:,j) = unperturbed_space_fixed_p - space_fixed_p;
        
        
        theta = deg2rad(-1 * delta_t * 360 / 86164);

        R3 = [cos(theta) -1*sin(theta) 0;
              sin(theta) cos(theta) 0;
              0 0 1];
        unperturbed_earth_fixed_p = R3 * unperturbed_space_fixed_p;
        earth_fixed_p = R3 * space_fixed_p;
        unperturbed_earth_fixed_v = R3 * unperturbed_space_fixed_v;
        earth_fixed_v = R3 * space_fixed_v;
        
        unperturbed_ret = xyz2blh(unperturbed_space_fixed_p(1),unperturbed_space_fixed_p(2),unperturbed_space_fixed_p(3));
        unperturbed_B = unperturbed_ret(1)
        unperturbed_L = unperturbed_ret(2)
        
        ret = xyz2blh(earth_fixed_p(1),earth_fixed_p(2),earth_fixed_p(3));
        B = ret(1)
        L = ret(2)
        
        diff_earth_fixed_L(j) = unperturbed_L - L;
        diff_earth_fixed_B(j) = unperturbed_B - B;
        

        R2 = [cos(pi/2 - phi_berlin) 0 sin(pi/2 - phi_berlin);
              0 1 0;
              -1*sin(pi/2 - phi_berlin) 0 cos(pi/2 - phi_berlin)];
        R3 = [cos(lambda_berlin) -1*sin(lambda_berlin) 0;
              sin(lambda_berlin) cos(lambda_berlin) 0;
              0 0 1];
        M1 = [ -1 0 0;
                0 1 0;
                0 0 1];
            
        topocentric_fixed_p = M1 * R2 * R3 * earth_fixed_p;
        unperturbed_topocentric_fixed_p = M1 * R2 * R3 * unperturbed_earth_fixed_p;
        topocentric_fixed_v = M1 * R2 * R3 * earth_fixed_v;
        unperturbed_topocentric_fixed_v = M1 * R2 * R3 * unperturbed_earth_fixed_v;
        topocentric_fixed_Az = rad2deg(atan2(topocentric_fixed_p(2), topocentric_fixed_p(1)));
        unperturbed_topocentric_fixed_Az = rad2deg(atan2(unperturbed_topocentric_fixed_p(2), unperturbed_topocentric_fixed_p(1)));
        topocentric_fixed_El = 90 - rad2deg(asin(topocentric_fixed_p(3)));
        unperturbed_topocentric_fixed_El = 90 - rad2deg(asin(unperturbed_topocentric_fixed_p(3)));

        diff_topocentric_fixed_p(:,j) = unperturbed_topocentric_fixed_p - topocentric_fixed_p;
        diff_topocentric_fixed_v(:,j) = unperturbed_topocentric_fixed_v - topocentric_fixed_v;
        diff_topocentric_fixed_Az(j) = unperturbed_topocentric_fixed_Az - topocentric_fixed_Az;
        diff_topocentric_fixed_El(j) = unperturbed_topocentric_fixed_El - topocentric_fixed_El;
        
        % perturbed 
        p(i) = p_a * (1 - p_e^2);%?
        c(i) = sqrt(p(i)*GM);%
        K(:,i) = [cos(p_OMEGA); sin(p_OMEGA); 0];
        C0(:,i) = [sin(p_OMEGA)*sin(p_i); -cos(p_OMEGA)*sin(p_i); cos(p_i)];
        P(:,i) = cos(p_omega)*K(:,i) + sin(p_omega)*(cross(C0(:,i),K(:,i)));
        Q(:,i) = -sin(p_omega)*K(:,i) + cos(p_omega)*(cross(C0(:,i),K(:,i)));
        
        
        %unperturbed
        p(i+1) = a * (1 - e^2);%?
        c(i+1) = sqrt(p(i+1)*GM);%
        K(:,i+1) = [cos(OMEGA); sin(OMEGA); 0];
        C0(:,i+1) = [sin(OMEGA)*sin(inc); -cos(OMEGA)*sin(inc); cos(inc)];
        P(:,i+1) = cos(omega)*K(:,i+1) + sin(omega)*(cross(C0(:,i+1),K(:,i+1)));
        Q(:,i+1) = -sin(omega)*K(:,i+1) + cos(omega)*(cross(C0(:,i+1),K(:,i+1)));
        
          
        for k=i:i+1
            f(k,j) = 2 * atan( tan(whatE(j)/2 ) * sqrt( (1+e)/(1+e) ) );
            % true anomaly of i.th satellite in j.th epoch
            if(f(k,j)<0)
                f(k,j) = f(k,j) + 2*pi;
            end % end of if
            
            r2(k,j) = p_a * ( 1 - p_e * cos( whatE(j) ) );
            r(k,j) = p_a*(1 - p_e^2)/( 1 + p_e*cos(f(k,j)) );
            Rp(:,k,j) = r(k,j) * ( cos(f(k,j))*P(:,k) +  sin(f(k,j))*Q(:,k) );
            Vp(:,k,j) = sqrt(c(k)/p(k)) * ( -sin(f(k,j))*P(:,k) +  ...
                ( e+cos(f(k,j)) ) *Q(:,k) );
            Rp_rot(:,k,j) = Rp(:,k,j)'*rotz(t(j,i)*omega_earth);
            Rp_rot_pol(:,k,j) = xyz2blh(Rp_rot(1,k,j),Rp_rot(2,k,j),...
                Rp_rot(3,k,j));
            XY_Hammer(:,k,j) = Hammer(Rp_rot_pol(1,k,j)*pi/180,...
                Rp_rot_pol(2,k,j)*pi/180);
        end
    end% end of for-loop for j
end% end of for-loop for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grs80 = referenceEllipsoid('grs80','km');
figure('Renderer','opengl')


ax = axesm('globe','Geoid',grs80,'Grid','on', ...
    'GLineWidth',1,'GLineStyle','-',...
    'Gcolor',[0.9 0.9 0.1],'Galtitude',100);
ax.Position = [0 0 1 1];
axis equal off
view(3)
load topo
geoshow(topo,topolegend,'DisplayType','texturemap')
demcmap(topo)
land = shaperead('landareas','UseGeoCoords',true);
plotm([land.Lat],[land.Lon],'Color','black')
rivers = shaperead('worldrivers','UseGeoCoords',true);
plotm([rivers.Lat],[rivers.Lon],'Color','blue')

% disp('------------------------------------------------------------------')
% disp('------------------------------------------------------------------')
%%-------------------------------------------------------------------------
for i=1:2
    for j=1:n
        x(j,1)=Rp(1,i,j);
        y(j,1)=Rp(2,i,j);
        z(j,1)=Rp(3,i,j);
    end
    hold on
    pk(i) = plot3( x/1000 , y/1000 , z/1000,colors(i,:));
    
end % end of for-loop for i

legend(pk,'CHAMP - perturbed','CHAMP - unperturbed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(1:n, diff_space_fixed_p);
axis([0 720 -inf inf]);
title('space-fixed: r - r_{perturbed}');
legend('x_{diff}','y_{diff}','z_{diff}');
ylabel('position in m');
xlabel('mean anomaly in degree');

figure;
plot(1:n, diff_earth_fixed_L, 1:n, diff_earth_fixed_B);
axis([0 720 -inf inf]);
title('Earth-fixed: L - L_{perturbed}, B - B_{perturbed}');
legend('L_{diff}','B_{diff}');
ylabel('longtitude and latitude in degree');
xlabel('mean anomaly in degree');

figure;
plot(1:n, diff_topocentric_fixed_p);
axis([0 720 -inf inf]);
title('Topocentric: r - r_{perturbed}');
legend('x_{diff}','y_{diff}', 'z_{diff}');
ylabel('position in m');
xlabel('mean anomaly in degree');

figure;
plot(1:n, diff_topocentric_fixed_v);
axis([0 720 -inf inf]);
title('Topocentric: V - V_{perturbed}');
legend('Vx_{diff}','Vy_{diff}', 'Vz_{diff}');
ylabel('velocity in m/s');
xlabel('mean anomaly in degree');

figure;
plot(1:n, diff_topocentric_fixed_Az, 1:n, diff_topocentric_fixed_El);
axis([0 720 -inf inf]);
title('Topocentric: Az - Az_{perturbed}, El - El_{perturbed}');
legend('Az_{diff}','El_{diff}');
ylabel('Azimuth and Elevation in degree');
xlabel('mean anomaly in degree');