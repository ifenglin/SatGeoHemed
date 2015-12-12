function ret= Azz(X)
% This function returns the Zenith distanc and Azimuth of vector
x = X(1);
y = X(2);
z = X(3);
nr=norm([x,y,z]);
if(nr~=0)
    x=x/nr;
    y=y/nr;
    z=z/nr;
end
zen = acos(z);
if(x~=0)
    temp = abs(atan(y/x));
    if(x>0&&y>=0)
        az = temp;
    elseif(x>0&&y<0)
        az = 2*pi - temp;
    elseif(x<0&&y<0)
        az = pi + temp;
    elseif(x<0&&y>0)
        az = pi - temp;
    end
else
    if(y>0)
        az = pi/2;
    elseif(y<0)
        az = 3*pi/2;
    else
        az = 0;
    end% end if
end% end if

ret=[zen;az];
end% end function













% -----------------------------
% x = X(1);
% y = X(2);
% z = X(3);
% nr=norm([x,y,z]);
% if(nr~=0)
%     x=x/nr;
%     y=y/nr;
%     z=z/nr;
% end
% Zz= acos(z);
% if(sin(Zz)~=0)
%     temp = abs(acos(x/sin(Zz)));
%     if(x>=0&&y>=0)
%         Az = temp;
%     elseif(x>0&&y<0)
%         Az = 2*pi - temp;
%     elseif(x<0&&y<0)
%         Az = pi + temp;
%     elseif(x<0&&y>0)
%         Az = pi - temp;
%     end
% else
%     Az=0;
% end
% ret=[Zz;Az];
% end% end function