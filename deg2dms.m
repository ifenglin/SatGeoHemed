function r=deg2dms(d)

% Calculations
g=fix(d);
m=fix((d-g)*60);
s=((d-g)*60-m)*60;
s=roundn(s,-4);
r=[g,m,s];
% r
