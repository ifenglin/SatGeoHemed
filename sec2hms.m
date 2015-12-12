function ret = sec2hms( t )
%convert second t hour minute second
% Hamed 11/28/2015 Berlin
t = t/3600;%  to hour
hour = floor(t);
temp = (t - hour)*60;
min = floor(temp);
sec = (temp - min)*60;
ret = [hour , min ,sec];
end

