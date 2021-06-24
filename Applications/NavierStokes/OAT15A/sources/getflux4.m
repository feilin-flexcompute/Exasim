function f = getflux4(udg,odg,param)
gam = param(1);
gam1 = gam - 1.0;
Re = param(2);
Pr = param(3);
Minf = param(4);
pmin = 0.2;
rmin = 0.2;
alpha = 1.0e3;
Re1 = 1/Re;
M2 = Minf^2;
c23 = 2.0/3.0;
avb = odg(1);
avr = odg(2);
avs = odg(3);
Re1 = Re1 + avs;
fc = Re1/(gam1*M2*Pr);
r = udg(1);
ru = udg(2);
rv = udg(3);
rE = udg(4);
rx = udg(5);
rux = udg(6);
rvx = udg(7);
rEx = udg(8);
ry = udg(9);
ruy = udg(10);
rvy = udg(11);
rEy = udg(12);
r = rmin + lmax(r-rmin,alpha);
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
rx = rx*dr;
ry = ry*dr;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv);
p = gam1*(rE-r*q);
p = pmin + lmax(p-pmin,alpha);
dp = atan(alpha*(p - pmin))/pi + (alpha*(p - pmin))/(pi*(alpha^2*(p - pmin)^2 + 1)) + 1/2;
h = E+p*r1;
fi = [ru, ru*uv+p, rv*uv, ru*h, ...
        rv, ru*vv, rv*vv+p, rv*h];
ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
qx = uv*ux + vv*vx;
px = gam1*(rEx - rx*q - r*qx);
px = px*dp;
Tx = gam*M2*(px*r - p*rx)*r1^2;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
qy = uv*uy + vv*vy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;
Ty = gam*M2*(py*r - p*ry)*r1^2;
% artificial viscosity for continuity equation
txx = (Re1)*c23*(2*ux - vy) + (avb+avr)*(ux+vy);
txy = (Re1)*(uy + vx);
tyy = (Re1)*c23*(2*vy - ux) + (avb+avr)*(ux+vy);
fv = [avr*rx, txx, txy, uv*txx + vv*txy + (fc+(avb+avr)*gam)*Tx, ...
      avr*ry, txy, tyy, uv*txy + vv*tyy + (fc+(avb+avr)*gam)*Ty];
f = fi+fv;
