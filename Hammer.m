function ret = Hammer( phi ,lambda )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    temp = sqrt( 1 + cos(phi)*cos(lambda/2) );
    x = 2 * sqrt(2) * cos(phi) * sin(lambda/2) / temp;
    y = sqrt(2) * sin(phi)/temp;
    ret = [x;y];
end

