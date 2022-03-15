function ftau = stabRoe2D(u1,u2,n,mu)
 
gam  = mu(1);
entropyfix = mu(18);
gam1 = gam-1.0;
ftau = 0*u1;

% regularization parameters
alpha = 1.0e3;
rmin = 1.0e-2;
hmin = 1.0e-4;

nx   = n(1);              
ny   = n(2);

% /* Conserved variables */
rr   = u1(1);            
rur  = u1(2);
rvr  = u1(3);
rEr  = u1(4);


rl   = u2(1);
rul  = u2(2);
rvl  = u2(3);
rEl  = u2(4);

rr = rmin + lmax(rr-rmin,alpha);
rr1  = 1.0/rr;
ur   = rur*rr1;
vr   = rvr*rr1;
Er   = rEr*rr1;
u2r  = ur*ur+vr*vr;
hr = gam*Er - 0.5*gam1*u2r;
hr = hmin + lmax(hr-hmin,alpha);

rl = rmin + lmax(rl-rmin,alpha);
rl1  = 1.0/rl;
ul   = rul*rl1;
vl   = rvl*rl1;
El   = rEl*rl1;
u2l  = ul*ul+vl*vl;
hl = gam*El - 0.5*gam1*u2l;
hl = hmin + lmax(hl-hmin,alpha);

%/* Roe's Averaging */
di    = sqrt(rr*rl1);      
d1    = 1.0/(di+1.0);
ui    = (di*ur+ul)*d1;
vi    = (di*vr+vl)*d1;
hi    = (di*hr+hl)*d1;
af    = 0.5*(ui*ui+vi*vi);
ci2   = (gam-1.0)*(hi-af);
ci    = sqrt(ci2);
uni   = ui*nx+vi*ny;

dr    = rr-rl;
dru   = rur-rul;
drv   = rvr-rvl;
drE   = rEr-rEl;

%     rlam1 = sqrt((uni+ci)*(uni+ci) + entropyfix);
%     rlam2 = sqrt((uni-ci)*(uni-ci) + 1e-6);
%     rlam3 = sqrt(uni*uni + 1e-6);

lam1 = uni+ci;
lam2 = uni-ci;
lam3 = uni;

rlam1 = 0.5*(lam1+entropyfix) + 0.5*(lam1-entropyfix).*tanh(1e3*(lam1-entropyfix));
rlam2 = 0.5*(lam2+entropyfix) + 0.5*(lam2-entropyfix).*tanh(1e3*(lam2-entropyfix));
rlam3 = 0.5*(lam3+entropyfix) + 0.5*(lam3-entropyfix).*tanh(1e3*(lam3-entropyfix));

s1    = 0.5*(rlam1+rlam2);
s2    = 0.5*(rlam1-rlam2);
al1x  = (gam-1.0)*(af*dr-ui*dru-vi*drv+drE);
al2x  = -uni*dr+dru*nx+drv*ny;
cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci);
cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x;

ftau(1)  = 0.5*(rlam3*dr+cc1);
ftau(2)  = 0.5*(rlam3*dru+cc1*ui+cc2*nx);
ftau(3)  = 0.5*(rlam3*drv+cc1*vi+cc2*ny);
ftau(4)  = 0.5*(rlam3*drE+cc1*hi+cc2*uni);