function ret =FiL2xy_merc(Fi,L)
dx=544.14;dy=296.28;mx=3.00243;my=3.11425;
x=(L*mx+dx);
y=(-Fi*my+dy);
ret = [x;y];