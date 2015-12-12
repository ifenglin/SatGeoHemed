function [r1,v1] = position_velocity(a,e,i,OMEGA,om,E)
GM=398600.44;

   E =E * pi / 180;
   i =i *pi / 180;
   OMEGA =OMEGA*pi / 180;
   om =om*pi / 180;
   f=2*atan(tan(E/2)*sqrt(1+e/1-e));
   r=a*(1-e*cos(E));
   p=r*(1+e*cos(f));
   c=sqrt(p*GM);
   K=[cos(OMEGA);sin(OMEGA);0];
   C0=[sin(OMEGA)*sin(i);-cos(OMEGA)*sin(i);cos(i)];
   P=cos(om)*K+sin(om)*cross(C0,K);
   Q=-sin(om)*K+cos(om)*cross(C0,K);
  
    
    r1=r*(cos(f)*P+sin(f)*Q);
    v1 = c/p*(-sin(f)*P+(e+cos(f))*Q);
  
end
