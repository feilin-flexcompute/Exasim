% Add Exasim to Matlab search path
close all
clear all

cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");
addpath('utilities');
load('HypersonicCylinder_p3_n3');

timestepend = 50000;
pde.soltime = timestepend:1000:timestepend; % steps at which solution are collected

% % get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd,'dataout');

udg = sol(:,1:4,:,:);
qdg = sol(:,5:12,:,:);
vdg = mesh.vdg;

npe = size(sol,1);
nelem = size(sol,3);
nt = length(pde.soltime);
sb = zeros(npe,nelem,nt);
AV = zeros(npe,nelem,nt);


gam = pde.physicsparam(1);
gam1 = gam-1;
Minf = pde.physicsparam(4);
Tinf = pde.physicsparam(9);
Tref = pde.physicsparam(10);
Twall = pde.physicsparam(11);
Tw = Twall/Tref * Tinf;
Tt = (1+0.5*gam1*Minf^2)*Tinf;
deltaT = (Tt - Tw);
pinf = 1/(gam*Minf^2);
rEinf = pinf/gam1 + 0.5;
Hinf = (rEinf + pinf);

rho = udg(:,1,:,:);
uv = udg(:,2,:,:)./rho;
vv = udg(:,3,:,:)./rho;
E = udg(:,4,:,:)./rho;
u2 = sqrt(uv.^2+vv.^2);
p = abs((gam-1)*(udg(:,4,:,:) - 0.5*(udg(:,2,:,:).*uv + udg(:,3,:,:).*vv)));
Ma = u2./sqrt(gam*p./udg(:,1,:,:));
T = (gam*Minf^2)*p./udg(:,1,:,:);
Tstagnation = ((1+0.5*gam1*Ma.^2).*T)/(1+0.5*gam1*Minf^2);
H = (E + p./rho)/Hinf;

for iElem = 1:nelem
    for iT = 1:nt
        for ipe = 1:npe
            [sb(ipe,iElem,iT),AV(ipe,iElem,iT)] = getavparams(udg(ipe,:,iElem,iT),qdg(ipe,:,iElem,iT),vdg(ipe,:,iElem),pde.physicsparam);
            
        end
    end
end

solvis = zeros(size(sol,1),11,size(sol,3),size(sol,4));
solvis(:,1:4,:,:) = sol(:,1:4,:,:);
solvis(:,5,:,:) = Ma;
solvis(:,6,:,:) = T;          %T/Tinf
solvis(:,7,:,:) = p; %Cp
solvis(:,8,:,:) = H;
solvis(:,9,:,:) = Tstagnation;
solvis(:,10,:,:) = sb;
solvis(:,11,:,:) = AV;

% % visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density",1,"energy",4,"Mach", 5, "Temperature", 6, "Pressure", 7, "Total Enthalpy", 8, "Stagnation Temperature", 9, "sensor", 10, "AV", 11};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2 3]}; % list of vector fields for visualization
xdg = vis(solvis,pde,mesh); % visualize the numerical solution
disp("Done!");

addpath('visualization')
[Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2(master,mesh,sol(:,:,:,end),pde.physicsparam,1,0,deltaT);
theta = atand(x(:,2)./x(:,1));

% Get sets of Experimental data for comparisons
ExpSt = csvread('visualization/LAURA_St.txt');
ExpCf = csvread('visualization/LAURA_Cf.txt');
ExpCp = csvread('visualization/LAURA_Cp.txt');

figure(1)
hold on
plot(theta,-Cp,'-b',ExpCp(:,1),ExpCp(:,2),'k--','LineWidth',2);
figure(2)
hold on
plot(theta,Ch,'-b',ExpSt(:,1),ExpSt(:,2),'k--','LineWidth',2);
figure(3)
hold on
plot(theta,-Cf,'-b',ExpCf(:,1),ExpCf(:,2),'k--','LineWidth',2);

