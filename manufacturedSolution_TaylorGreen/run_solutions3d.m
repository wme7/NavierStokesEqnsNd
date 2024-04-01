% Test Exact solutions
addpath('../common/');
global gamma mu Pr %#ok<GVMIS>

% Aceptable numerical error to stablish solution validity:
TOL = 1E-11;

% Domain length
L = 1.0;

% Fluid parameters
gamma = 1.4;
mu = 0.025;
Pr = 0.72;

% 3D - mesh
[x,y,z] = meshgrid(linspace(0,L,25));

%% Compute Manufactured solution
figure(3); t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qt,qx,qy,qz,qxx,qyy,qzz,s] = manufacturedSolution3d(x(:),y(:),z(:),L,t);

    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns3d(q,qt,qx,qy,qz,qxx,qyy,qzz,s,TOL)), return; end

    % Update figure
    plotNavierStokesEqns3d(x,y,z,q,t);
    % plotVorticityNavierStokesEqns3d(x,y,z,q,t)
    drawnow;
end

%% Compute Taylor-Green vortex (Pseudo-)solution
figure(4); t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qt,qx,qy,qz,qxx,qyy,qzz,s] = taylorGreenVortex3d(x(:),y(:),z(:),L,t);

    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns3d(q,qt,qx,qy,qz,qxx,qyy,qzz,s,TOL)), return; end

    % Update figure
    plotNavierStokesEqns3d(x,y,z,q,t);
    % plotVorticityNavierStokesEqns3d(x,y,z,q,t)
    drawnow;
end

%% Compute Taylor-Green vortex solution (Antuono's Model)
figure(5); t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qt,qx,qy,qz,qxx,qyy,qzz,s] = taylorGreenVortex3d_AntuonoModel(x(:),y(:),z(:),L,t);

    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns3d(q,qt,qx,qy,qz,qxx,qyy,qzz,s,TOL)), return; end

    % Update figure
    % plotNavierStokesEqns3d(x,y,z,q,t);
    plotVorticityNavierStokesEqns3d(x,y,z,q,t)
    drawnow;
end
