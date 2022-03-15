% Add Exasim to Matlab search path
load('../HypersonicCylinder');

pde.dt = 1e-4*ones(1,2700)*2;   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = 200;          % solution is saved every 10 time steps
pde.soltime = 200:200:length(pde.dt); % steps at which solution are collected

gam = pde.physicsparam(1);
Minf = pde.physicsparam(4);
Tinf = pde.physicsparam(9);
Tref = pde.physicsparam(10);
Twall = pde.physicsparam(11);
Tw = Twall/Tref * Tinf;
Tt = (1+0.5*(gam-1)*Minf^2)*Tinf;
deltaT = (Tt - Tw);

% get solution from output files in dataout folder
sol = getsolution(['../dataout/out_t' num2str(pde.soltime(end))],dmd,master.npe);

[Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2(master,mesh,sol(:,:,:,end),pde.physicsparam,1,0,deltaT,1);

figure(1)
plot(atand(x(:,2)./x(:,1)),-Cp,'o');
figure(2)
plot(atand(x(:,2)./x(:,1)),Ch,'o');
figure(3)
plot(atand(x(:,2)./x(:,1)),-Cf,'o');