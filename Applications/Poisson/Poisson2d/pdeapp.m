% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
L = 1e4; %domain size;
pde.physicsparam = [1,L];    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.scaling = L;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(17,17,1,1);
mesh.p = L*mesh.p;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-L)<1e-8, @(p) abs(p(2,:)-L)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master] = exasim(pde,mesh);


%Compute error
dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
mesh.dgnodes = dgnodes;

uh = sol(:,1,:);
qh = sol(:,2:3,:);

ue = @(x) sin(pi*x(:,1)/L).*sin(pi*x(:,2)/L);
qe = @(x) -pi*[cos(pi*x(:,1)/L).*sin(pi*x(:,2)/L), cos(pi*x(:,2)/L).*sin(pi*x(:,1)/L)]/L;

[~,erru] = computeerror(uh,mesh,master,ue);
[~,errq] = computeerror(qh,mesh,master,qe);

fprintf('L2 error u: %g\n',erru);
fprintf('L2 error q: %g\n',norm(errq));

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

        
