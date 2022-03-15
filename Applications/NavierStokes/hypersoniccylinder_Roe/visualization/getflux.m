function f = getflux(udg,param)
gam = param(1);
gam1 = gam - 1.0;
Re = param(2);
Pr = param(3);
Minf = param(4);
Tref = param(10);
muRef = 1/Re;
Tinf = 1/(gam*gam1*Minf^2);
c23 = 2.0/3.0;

r = udg(:,1);
ru = udg(:,2);
rv = udg(:,3);
rE = udg(:,4);
rx = udg(:,5);
rux = udg(:,6);
rvx = udg(:,7);
rEx = udg(:,8);
ry = udg(:,9);
ruy = udg(:,10);
rvy = udg(:,11);
rEy = udg(:,12);

r1 = 1./r;
uv = ru.*r1;
vv = rv.*r1;
q = 0.5*(uv.*uv+vv.*vv);
p = gam1*(rE-r.*q);

ux = (rux - rx.*uv).*r1;
vx = (rvx - rx.*vv).*r1;
qx = uv.*ux + vv.*vx;
px = gam1*(rEx - rx.*q - r.*qx);
Tx = 1/gam1*(px.*r - p.*rx).*r1.^2;
uy = (ruy - ry.*uv).*r1;
vy = (rvy - ry.*vv).*r1;
qy = uv.*uy + vv.*vy;
py = gam1*(rEy - ry.*q - r.*qy);
Ty = 1/gam1*(py.*r - p.*ry).*r1.^2;


T = p./(gam1.*r);
Tphys = Tref/Tinf * T;
mu = getViscosity(muRef,Tref,Tphys,1);
fc = mu*gam/Pr;

txx = mu.*c23.*(2*ux - vy);
txy = mu.*(uy + vx);
tyy = mu.*c23.*(2*vy - ux);
sz = size(txx);


f = [zeros(sz), txx, txy, fc.*Tx, ...
      zeros(sz), txy, tyy, fc.*Ty];
  
  
  
