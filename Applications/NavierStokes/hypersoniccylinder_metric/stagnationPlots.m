close all

% addpath('export_fig/');
addpath('utilities');
data = readtable('data.csv'); %export data along x line from paraview

%the indices may need to be adapted
LineMa = table2array(data(:,1));
Linep = table2array(data(:,2));
LineT = table2array(data(:,3));
LineH = table2array(data(:,4));
Linerho = table2array(data(:,5));
Linev = table2array(data(:,6))./Linerho;
x = -1 - table2array(data(:,10));

figure
hold on
plot(x,Linev,'k','LineWidth',2);
plot(x,Linep,'b','LineWidth',2);
plot(x,LineH,'r','LineWidth',2);
axis([0 0.6 0 2])

legend('$v_x/v_{\infty}$','$p/p_{\infty}$','$H/H_{\infty}$','location','northeast','Interpreter','latex','FontSize',18)
grid on
grid minor
xlabel('$d/L$','Fontsize',18,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',18)
set(gcf,'color','w');
box on
% export_fig 'StagnationQuantities1Roe100.pdf' -r500

figure
plot(x,LineT,'r','LineWidth',2);
axis([0 0.6 0 70])
legend('$T/T_{\infty}$','location','northeast','Interpreter','latex','FontSize',18)
grid on
% grid minor
xlabel('$d/L$','Fontsize',18,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',18)
set(gcf,'color','w');
box on

% create smaller axes in top right, and plot on it
axes('Position',[.2 .2 .25 .3])
box on
plot(x,LineT,'r','LineWidth',2);
axis([0 0.012 0 65])
grid on
grid minor
set(gca,'TickLabelInterpreter','latex','FontSize',12)
title('Zoom view','Fontsize',14,'interpreter','latex')
% export_fig 'StagnationQuantities2Roe100.pdf' -r500
