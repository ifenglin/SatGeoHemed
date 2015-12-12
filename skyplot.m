function ret = skyplot( Az, zen ,cut_angle,time)
% % This function draws Skyplot with given arguments:
% % Az        : Azimuth         (in radian)
% % Zen       : Zenith angle    (in radian)
% % cut_angle : Cun angle       (in radian)
% % time      : time interval from 0h time(hour)
% %             if the time is 0 then two plots will not show up
% % Auther : Hamed 11/26/2015 berlin
% % Test   : ok
figure;
t =[ 0, 2*pi];
p = polar(t, [90,90],'w');
Az2 = [];
zen2 = [];
for i=1:length(Az)
    if(zen(i)<=(pi/2-cut_angle))
        Az2 = [Az2; Az(i)];
        zen2 = [zen2; zen(i)];
    end % end if
end% end of for
set(p, 'Visible', 'off')
hold on
p2 = polar(Az2,zen2*180/pi,'.b');
view([90 -90])

hHiddenText = findall(gca,'type','text');
Angles = 0 : 30 : 330;
hObjToDelete = zeros( length(Angles)-4, 1 );
k = 0;
for ang = Angles
   hObj = findall(hHiddenText,'string',num2str(ang));
   switch ang
   case 0
      set(hObj,'string','North','HorizontalAlignment','Left');
   case 90
      set(hObj,'string','East','VerticalAlignment','Bottom');
   case 180
      set(hObj,'string','South','HorizontalAlignment','Right');
   case 270
      set(hObj,'string','West','VerticalAlignment','Top');
   otherwise
      k = k + 1;
      hObjToDelete(k) = hObj;
   end
end
% delete( hObjToDelete(hObjToDelete~=0) );

hold on;
temp = [0:.1:2*pi+.1];
p3 = polar(temp,(pi/2-cut_angle)*ones(size(temp))*180/pi,'r');
legend([p2,p3],'skyplot','cut angle')
%%%%%%%%%%%%%%%%%%%%%%%% Visibility %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(time~=0)
    figure;
    plot(linspace(0,time,length(Az)) , Az*180/pi)
    grid on
    title('Azimuth in time')
    xlabel('Time (hour)')
    ylabel('Azimuth')
    
    figure;
    x = [0, 0, time,time];
    y = [cut_angle*180/pi, -95,-95, 0];
    col = [.9 .9 .9];
    p4 =fill(x, y, col);
 
    
    hold on
    p5 = plot(linspace(0,time,length(zen)) , (pi/2-zen)*180/pi);
    grid on
    title('Elevation angle in time')
    xlabel('Time (hour)')
    ylabel('Elevation Angle')
    legend([p4,p5],'Invisibility','Elevation');

end% end if

ret = 1;


end% end of function

