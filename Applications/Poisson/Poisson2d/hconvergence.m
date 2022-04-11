% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel2"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.physicsparam = 1;    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.HDGflag = 0;
pde.extUhat = 0;

% Choose computing platform and set number of processors
%pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 2;           % number of MPI processors

ue = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
qe = @(x) -pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)), cos(pi*x(:,2)).*sin(pi*x(:,1))];

nOrders = 1;
nMeshes = 4;

ErrorU = zeros(nOrders,nMeshes);
ErrorQ = zeros(nOrders,nMeshes);

ErrorU2 = zeros(nOrders,nMeshes);
ErrorQ2 = zeros(nOrders,nMeshes);

for iOrder = 1:nOrders
    for jMesh = 1:nMeshes
        mesh = initializemesh("src");
        
        pde.porder = iOrder;             % polynomial degree
        N = 2^(jMesh) + 1;
        [mesh.p,mesh.t] = squaremesh(N,N,N,0);
        % expressions for domain boundaries
        mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
        mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

        % call exasim to generate and run C++ code to solve the PDE model
        [sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

        dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
        x = dgnodes(:,1,:); y = dgnodes(:,2,:);
        uexact = sin(pi*x).*sin(pi*y);           % exact solution
        uh = sol(:,1,:);
        qh = sol(:,2:3,:);
        fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
        
        mesh.dgnodes = dgnodes;
        [erru,rerru] = computeerror(uh,mesh,master,ue);
        [errq,rerrq] = computeerror(qh,mesh,master,qe);
        ErrorU(iOrder,jMesh) = erru;
        ErrorU2(iOrder,jMesh) = rerru;
        errq = norm(errq);
        rerrq = norm(rerrq);
        ErrorQ(iOrder,jMesh) = errq;
        ErrorQ2(iOrder,jMesh) = rerrq;
        
        fprintf('L2 error u: %g\n',erru);
        fprintf('L2 relative error u: %g\n',rerru);
        fprintf('L2 error q: %g\n',errq);
        fprintf('L2 relative error q: %g\n',rerrq);
        
        rmdir dataout s
        rmdir datain s
        rmdir app s
    end
end

slopeU = log(ErrorU2(:,end)./ErrorU2(:,end-1))/log(0.5);
slopeQ = log(ErrorQ2(:,end)./ErrorQ2(:,end-1))/log(0.5);

disp('Orders of convergence u:'); disp(slopeU');
disp('Orders of convergence q:'); disp(slopeQ');
