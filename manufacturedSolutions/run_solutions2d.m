% Test Exact solutions
addpath('../common/');
global gamma mu Pr %#ok<GVMIS>

% Aceptable numerical error to stablish solution validity:
TOL = 1E-12;

% Domain length
L = 1.0;

% Fluid parameters
gamma = 1.4;
mu = 0.025;
Pr = 0.72;

% 2D - mesh
[x,y] = meshgrid(linspace(0,L,25),linspace(-L/2,L/2,25));

%% Compute Manufactured solution
figure(1); t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd
    
    % Compute exact solution
    [q,qt,qx,qy,qxx,qyy,s] = manufacturedSolution2d(x(:),y(:),L,t);
    
    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns2d(q,qt,qx,qy,qxx,qyy,s,TOL)), return; end

    % Update figure
    mesh(x,y,reshape(q(:,4),size(x))); axis square; zlim([0,2]);
    % plotNavierStokesEqns2d(x,y,q,t);
    % plotVorticityNavierStokesEqns2d(x,y,q,t);
    drawnow;
end

%% Compute Taylor-Green vortex solution
figure(2); t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd
    
    % Compute exact solution
    [q,qt,qx,qy,qxx,qyy,s] = taylorGreenVortex2d(x(:),y(:),L,t);
    
    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns2d(q,qt,qx,qy,qxx,qyy,s,TOL)), return; end

    % Update figure
    % plotNavierStokesEqns2d(x,y,q,t);
    plotVorticityNavierStokesEqns2d(x,y,q,t);
    drawnow;
end