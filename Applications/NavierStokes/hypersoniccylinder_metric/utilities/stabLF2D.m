function ftau = stabLF2D(u1,u2,n,mu)
 
gam  = mu(1);
gam1 = gam-1.0;
nx   = n(1);              
ny   = n(2);

% regularization parameters
alpha = 1.0e3;
rmin = 1.0e-2;
hmin = 1.0e-4;

% define intermediate state
uhat = 0.5*(u1+u2);

% conserved variables
rh = uhat(1);
ruh = uhat(2);
rvh = uhat(3);
rEh = uhat(4);

rh = rmin + lmax(rh-rmin,alpha);
rh1  = 1.0/rh;
uh   = ruh*rh1;
vh   = rvh*rh1;
Eh = rEh*rh1;
u2h  = uh*uh+vh*vh;
hh = gam*Eh - 0.5*gam1*u2h;
hh = hmin + lmax(hh-hmin,alpha);
ch = sqrt(gam1*(hh-0.5*u2h));

vnh = uh*nx + vh*ny;
absvnh = vnh*tanh(alpha*vnh); 
lam_max = absvnh + ch;

ftau = 0.5*lam_max*(u1-u2);

