% Add Exasim to Matlab search path
load('../HypersonicCylinder');

gam = pde.physicsparam(1);
Minf = pde.physicsparam(4);
Tref = pde.physicsparam(10);
Twall = pde.physicsparam(11);


% Get some values and dimension first :
pinf = 1/(gam*Minf^2);
% Total temperature
Tinf = pinf/(gam-1);
Ttinf = Tinf * (1 + (gam-1)/2 * Minf^2);
% Delta between Freestream total and Wall temperature
TisoW = Twall/Tref * Tinf;
deltaT = Ttinf - TisoW;
% Id of the isothermal wall
wid = 1;
% Some dimensions
npv = (pde.porder+1)^2;

% % Get the Master structure
% master = mkmaster(mesh,2*porder);
master.nd = master.dim;

% Path to directory


% Initial and final times, and dsnap
tstart = 5000;
tend   = 5000;
dsnap  = 25000;
tshift = 0;

set(0,'defaultlinelinewidth',2.)

% Get sets of Experimental data for comparisons
ExpSt = csvread('LAURA_St.txt');
ExpCf = csvread('LAURA_Cf.txt');
ExpCp = csvread('LAURA_Cp.txt');

for itime=tstart:dsnap:tend

    fn1 =sprintf('%06d.png',itime);
    solfile = ['../dataout/out_t' num2str(itime)];
%     u = getsolution(solfile, length(elempart), npv, 12, elempart);
    u = getsolution(['../dataout/out_t' num2str(timeend)],dmd,master.npe);

    
    % Get Cp, Cf and Ch Coefficients
    [Cp,Cf,x,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata2(master,mesh,u,pde.physicsparam,wid,0,deltaT);
    
    % Get polar coordinates
    theta = atan2(x(:,2),-x(:,1))*180./pi;
    
    % Plot Cp (Prerssure Coeff) comparisons
    fig = figure('visible','off'); clf;
    set(fig, 'Position', [0 0 600 400])
    clf; plot(theta,-Cp,ExpCp(:,1),ExpCp(:,2),'k--'); xlim([-100 100]); ylim([0. 2.]); grid on;
    xlabel('\theta'); ylabel('Cp(\theta)');title('Pressure Coefficient');
    legend('Numerical p3','LAURA','Location','south');
%     fn = [path_dir 'pics/Cp' fn1];
%     print(fig,fn,'-dpng');
        
    % Plot Cf (Friction Coeff) comparisons
    fig = figure('visible','off'); clf;
    set(fig, 'Position', [0 0 600 400])
    clf; plot(theta,Cf,ExpCf(:,1),ExpCf(:,2),'k--'); xlim([-100 100]); ylim([-6e-3 6e-3]); grid on;
    xlabel('\theta'); ylabel('Cf(\theta)');title('Friction Coefficient');
    legend('Numerical p3','LAURA','Location','southeast');
%     fn = [path_dir 'pics/Cf' fn1];
%     print(fig,fn,'-dpng');
    
    % Plot Ch (Stanton Number) comparisons
    fig = figure('visible','off'); clf;
    set(fig, 'Position', [0 0 600 400])
    clf; plot(theta,Ch,ExpSt(:,1),ExpSt(:,2),'k--'); xlim([-100 100]); ylim([0e-3 9e-3]); grid on;
    xlabel('\theta'); ylabel('Ch(\theta)');title('Stanton number');
    legend('Numerical p3','LAURA','Location','south');
%     fn = [path_dir 'pics/Ch' fn1];
%     print(fig,fn,'-dpng');
    
    %close all;
   
end

% close all;

return;

% Read data in text file
%fileID = fopen('compar_data/Exp_St.txt','r');
%[ExpSt,~] = fscanf(fileID,['%f' ',']);
%fclose(fileID);

