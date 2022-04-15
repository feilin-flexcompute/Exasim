addpath('../');
close all
[p,t,dgnodes] = mkmesh_circincirc_adapted(3,81,241,1,2,3);
% [p,t,dgnodes] = mkmesh_circincirc_Ma17b(3,201,201,1,2,3);
meshplot2D(p,t)
% mesh.p = mesh.p';
% meshplot(mesh);