close all

% Add Exasim to Matlab search path
load('../HypersonicCylinder_p3');
% Get sets of Experimental data for comparisons
ExpSt = csvread('LAURA_St.txt');
ExpCf = csvread('LAURA_Cf.txt');
ExpCp = csvread('LAURA_Cp.txt');

timeend = 80000;

gam = pde.physicsparam(1);
Minf = pde.physicsparam(4);
Tinf = pde.physicsparam(9);
Tref = pde.physicsparam(10);
Twall = pde.physicsparam(11);
Tw = Twall/Tref * Tinf;
Tt = (1+0.5*(gam-1)*Minf^2)*Tinf;
deltaT = (Tt - Tw);

% get solution from output files in dataout folder
sol = getsolution(['../dataout/out_t' num2str(timeend)],dmd,master.npe);

[Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2(master,mesh,sol(:,:,:,end),pde.physicsparam,1,0,deltaT);

theta = atand(x(:,2)./x(:,1));

figure(1)
hold on
plot(theta,-Cp,'-r',ExpCp(:,1),ExpCp(:,2),'k--','LineWidth',2);
figure(2)
hold on
plot(theta,Ch,'-r',ExpSt(:,1),ExpSt(:,2),'k--','LineWidth',2);
figure(3)
hold on
plot(theta,-Cf,'-r',ExpCf(:,1),ExpCf(:,2),'k--','LineWidth',2);