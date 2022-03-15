% Add Exasim to Matlab search path
load('HypersonicCylinder');

timestepend = 2600;

pde.dt = 1e-4*ones(1,timestepend)*2;   % time step sizes
pde.soltime = 200:200:timestepend; % steps at which solution are collected

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);
udg = sol(:,1:4,:,:);
qdg = sol(:,5:12,:,:);
vdg = mesh.vdg;

npe = master.npe;
nelem = dmd{1}.elempartpts;
nt = length(pde.soltime);
avField = zeros(npe,nelem,nt);

gam = pde.physicsparam(1);
Minf = pde.physicsparam(4);
Tinf = pde.physicsparam(9);
Tref = pde.physicsparam(10);
Twall = pde.physicsparam(11);
Tw = Twall/Tref * Tinf;
Tt = (1+0.5*(gam-1)*Minf^2)*Tinf;
deltaT = (Tt - Tw);
pinf = 1/(gam*Minf^2);

uv = udg(:,2,:,:)./udg(:,1,:,:);
vv = udg(:,3,:,:)./udg(:,1,:,:);
u2 = sqrt(uv.^2+vv.^2);
p = abs((gam-1)*(udg(:,4,:,:) - 0.5*(udg(:,2,:,:).*uv + udg(:,3,:,:).*vv)));
Ma = u2./sqrt(gam*p./udg(:,1,:,:));
T = (gam*Minf^2)*p./udg(:,1,:,:);

for iElem = 1:nelem
    for iT = 1:nt
        for ipe = 1:npe
            avField(ipe,iElem,iT) = getavfield2d(udg(ipe,:,iElem,iT),qdg(ipe,:,iElem,iT),vdg(ipe,:,iElem),pde.physicsparam);
        end
    end
end

solvis = zeros(size(sol,1),4,size(sol,3),size(sol,4));
solvis(:,1,:,:) = Ma;
solvis(:,2,:,:) = T;          %T/Tinf
solvis(:,3,:,:) = 2*(p-pinf); %Cp
solvis(:,4,:,:) = avField;

% % visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"Mach", 1, "Temperature", 2, "Pressure", 3, "AV", 4};  % list of scalar fields for visualization
pde.visvectors = {}; % list of vector fields for visualization
xdg = vis(solvis,pde,mesh); % visualize the numerical solution
disp("Done!");

% [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2(master,mesh,sol(:,:,:,end),pde.physicsparam,1,0,deltaT,1);
% 
% figure(1)
% plot(atand(x(:,2)./x(:,1)),-Cp,'o');
% figure(2)
% plot(atand(x(:,2)./x(:,1)),Ch,'o');
% figure(3)
% plot(atand(x(:,2)./x(:,1)),-Cf,'o');
% 
