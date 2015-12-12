function ret =FiL2xy_hamm(Fi,L)
dx=640;dy=321;mx=224.540;my=228.582;
x=(L*mx+dx);
y=(-Fi*my+dy);
ret = [x;y];